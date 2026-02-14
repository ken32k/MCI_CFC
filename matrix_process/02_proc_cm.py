"""Compute SC-derived communication metrics for subjects."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Sequence

import numpy as np
from joblib import Parallel, delayed
from netneurotools import metrics as nmtr
from scipy.sparse import csr_matrix


@dataclass(frozen=True)
class PipelineConfig:
    """Typed container for frequently reused configuration values."""

    home_dir: Path = Path("/public/home/baishw/WMH_MCI")
    mica_subdir: str = "mica_out/micapipe_v0.2.0"
    data_subdir: str = "data"
    atlas_nodes: int = 200
    tract_labels: Sequence[str] = ("wb", "wmh", "int")
    threads: int = 12
    gammas: Sequence[float] = (0.25, 0.5, 1.0, 2.0, 4.0)
    markov_times: Sequence[float] = tuple(float(t) for t in range(1, 11))

    @property
    def atlas_name(self) -> str:
        return f"schaefer-{self.atlas_nodes}"

    @property
    def mica_out_dir(self) -> Path:
        return self.home_dir / self.mica_subdir

    @property
    def data_out_dir(self) -> Path:
        return self.home_dir / self.data_subdir


def sc_filter(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Remove nodes that have no SC connections to avoid singular metrics."""

    non_zero_mask = np.any(matrix != 0.0, axis=0)
    filtered = matrix[non_zero_mask][:, non_zero_mask]
    return filtered, non_zero_mask


def insert_nans(matrix: np.ndarray, non_zero_mask: np.ndarray, matrix_size: int) -> np.ndarray:
    """Restore a reduced matrix to its original shape, padding missing entries with NaNs."""

    restored = np.full((matrix_size, matrix_size), np.nan)
    restored[np.ix_(non_zero_mask, non_zero_mask)] = matrix
    return restored


def ensure_output_folders(base_dir: Path, subfolders: Iterable[str]) -> None:
    """Create metric-specific directories on demand to avoid race conditions."""

    for folder in subfolders:
        (base_dir / folder).mkdir(parents=True, exist_ok=True)


def purge_existing_outputs(metrics_dir: Path) -> None:
    """Delete stale metric files to guarantee a clean run."""

    if not metrics_dir.exists():
        return

    for file in metrics_dir.rglob("*.npy"):
        file.unlink()


def load_structural_connectomes(sub_id: str, config: PipelineConfig) -> Dict[str, np.ndarray]:
    """Load all SC variations for a subject, skipping missing entries."""

    connectomes: Dict[str, np.ndarray] = {}
    for tract_label in config.tract_labels:
        sc_path = config.data_out_dir / "scs" / f"{sub_id}_{tract_label}_sc-wei_{config.atlas_name}.npy"
        if not sc_path.exists():
            print(f"[Warn] Missing SC file for {sub_id} ({tract_label}).", flush=True)
            return {}
        connectomes[tract_label] = np.load(sc_path)
    return connectomes


def save_metric(array: np.ndarray, destination: Path) -> None:
    """Persist metric arrays ensuring the parent directories exist."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    np.save(destination, array)


def build_cost_matrix(sc_matrix: np.ndarray, gamma: float) -> np.ndarray:
    """Generate the cost matrix used across PL, SI, and PT metrics."""

    with np.errstate(divide="ignore"):
        return csr_matrix(sc_matrix).power(-gamma).toarray()


def process_single_connectome(
    sub_id: str,
    tract_label: str,
    sc_matrix: np.ndarray,
    config: PipelineConfig,
) -> None:
    """Compute and persist all cm for a single subject/tract pair."""
    # SC filter 
    filtered_sc, active_nodes = sc_filter(sc_matrix)
    matrix_size = config.atlas_nodes

    # CO
    try:
        communicability = nmtr.communicability_wei(filtered_sc)
        save_metric(
            np.array([insert_nans(communicability, active_nodes, matrix_size)]),
            config.data_out_dir
            / "cm"
            / "co-wei"
            / f"{sub_id}_{tract_label}_co-wei_{config.atlas_name}.npy",
        )
    except Exception as exc:  # pylint: disable=broad-except
        print(
            f"[Error] Communicability failed for {sub_id} ({tract_label}): {exc}",
            flush=True,
        )
        
    # --------------------------------------------------------------
    # FG
    flow_graphs = []
    for markov_time in config.markov_times:
        try:
            flow_graph = nmtr.flow_graph(filtered_sc, t=markov_time)
        except Exception as exc:  # pylint: disable=broad-except
            flow_graph = np.zeros_like(filtered_sc)
            print(
                f"[Error] Flow-graph failed for {sub_id} ({tract_label}, t={markov_time}): {exc}",
                flush=True,
            )
        flow_graphs.append(insert_nans(flow_graph, active_nodes, matrix_size))

    save_metric(
        np.asarray(flow_graphs),
        config.data_out_dir
        / "cm"
        / "fg-wei"
        / f"{sub_id}_{tract_label}_fg-wei_{config.atlas_name}.npy",
    )

    # --------------------------------------------------------------
    # PT, PL, SI
    path_lengths, path_transitivity, search_information = [], [], []
    for gamma in config.gammas:
        cost_matrix = build_cost_matrix(filtered_sc, gamma)

        try:
            pl_matrix = nmtr.distance_wei_floyd(cost_matrix)[0]
        except Exception:  # pylint: disable=broad-except
            pl_matrix = np.zeros_like(cost_matrix)
        path_lengths.append(insert_nans(pl_matrix, active_nodes, matrix_size))

        try:
            pt_matrix = nmtr.path_transitivity(cost_matrix)
            pt_matrix = (pt_matrix + pt_matrix.T) / 2
        except Exception:  # pylint: disable=broad-except
            pt_matrix = np.zeros_like(cost_matrix)
        path_transitivity.append(insert_nans(pt_matrix, active_nodes, matrix_size))

        try:
            si_matrix = nmtr.search_information(filtered_sc, cost_matrix)
            si_matrix = (si_matrix + si_matrix.T) / 2
        except Exception:  # pylint: disable=broad-except
            si_matrix = np.zeros_like(cost_matrix)
        search_information.append(insert_nans(si_matrix, active_nodes, matrix_size))

    # Save PL, PT, SI
    save_metric(
        np.asarray(path_lengths),
        config.data_out_dir
        / "cm"
        / "pl-wei"
        / f"{sub_id}_{tract_label}_pl-wei_{config.atlas_name}.npy",
    )
    save_metric(
        np.asarray(path_transitivity),
        config.data_out_dir
        / "cm"
        / "pt-wei"
        / f"{sub_id}_{tract_label}_pt-wei_{config.atlas_name}.npy",
    )
    save_metric(
        np.asarray(search_information),
        config.data_out_dir
        / "cm"
        / "si-wei"
        / f"{sub_id}_{tract_label}_si-wei_{config.atlas_name}.npy",
    )

def par_proc_mat(sub_id: str, config: PipelineConfig) -> None:
    """Compute SC-derived communication metrics for a single subject."""
    print(f"[Info] Processing subject {sub_id}", flush=True)

    connectomes = load_structural_connectomes(sub_id, config)
    if not connectomes:
        print(f"[Error] {sub_id}: Missing SC files. Skipping.", flush=True)
        return

    for tract_label, sc_matrix in connectomes.items():
        process_single_connectome(sub_id, tract_label, sc_matrix, config)


def main() -> None:
    """Main pipeline execution function."""
    config = PipelineConfig()
    print(f"Run {__file__} with atlas {config.atlas_name}", flush=True)

    metrics_dir = config.data_out_dir / "cm"
    ensure_output_folders(
        metrics_dir,
        ("co-wei", "pl-wei", "pt-wei", "si-wei", "fg-wei"),
    )
    purge_existing_outputs(metrics_dir)

    if not config.mica_out_dir.exists():
        raise FileNotFoundError(f"MICA directory missing: {config.mica_out_dir}")

    subjects = sorted(
        entry.name for entry in config.mica_out_dir.iterdir() if entry.is_dir()
    )

    Parallel(n_jobs=config.threads)(
        delayed(par_proc_mat)(subject, config) for subject in subjects
    )


if __name__ == "__main__":
    main()
