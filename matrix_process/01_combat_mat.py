"""Harmonize structural and functional connectivity matrices with neuroCombat."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from neuroCombat import neuroCombat


@dataclass(frozen=True)
class CombatConfig:
    """Container for paths and processing hyper-parameters."""

    home_dir: Path = Path("/public/home/baishw/WMH_MCI")
    atlas_nodes: int = 200
    threshold: float = 0.25 # Proportional threshold for SC matrices
    tract_labels: Sequence[str] = ("wb", "wmh", "int")
    reference_batch: str = "1"

    @property
    def atlas_name(self) -> str:
        return f"schaefer-{self.atlas_nodes}"

    @property
    def data_dir(self) -> Path:
        return self.home_dir / "data"

    @property
    def sc_raw_dir(self) -> Path:
        return self.data_dir / "scs_raw"

    @property
    def sc_dir(self) -> Path:
        return self.data_dir / "scs"

    @property
    def fc_raw_dir(self) -> Path:
        return self.data_dir / "fcs_raw"

    @property
    def fc_dir(self) -> Path:
        return self.data_dir / "fcs"

    @property
    def combat_dir(self) -> Path:
        return self.data_dir / "combat"

    @property
    def subject_info_path(self) -> Path:
        return self.data_dir / "sub_list_beforePSM_WMHstat.csv"


def threshold_proportional(matrix: np.ndarray, proportion: float, copy: bool = True) -> np.ndarray:
    """Preserve the strongest proportion of weights following Brain Connectivity Toolbox semantics."""

    if not 0 <= proportion <= 1:
        raise ValueError("Threshold must be in range [0, 1].")

    working = matrix.copy() if copy else matrix
    n_nodes = working.shape[0]
    is_symmetric = np.allclose(working, working.T)

    if is_symmetric:
        working[np.tril_indices(n_nodes)] = 0
        divisor = 2
    else:
        divisor = 1

    indices = np.where(working)
    sorted_idx = np.argsort(working[indices])[::-1]
    keep = int(round((n_nodes * n_nodes - n_nodes) * proportion / divisor))
    working[(indices[0][sorted_idx][keep:], indices[1][sorted_idx][keep:])] = 0

    if is_symmetric:
        working[:, :] = working + working.T

    return working


def upper_triangle_vector(matrix: np.ndarray, log_scale: bool = False) -> np.ndarray:
    """Extract the strictly upper-triangular entries as a vector with optional log transform."""

    tri_indices = np.triu_indices(matrix.shape[0], k=1)
    values = matrix[tri_indices]
    if log_scale:
        values = np.log(np.clip(values, a_min=1e-12, a_max=None))
    return values


def vector_to_symmetric_matrix(vector: np.ndarray, nodes: int, *, exp_transform: bool = False, zero_diagonal: bool = False) -> np.ndarray:
    """Reconstruct a symmetric matrix from a vectorized upper triangle."""

    data = np.exp(vector) if exp_transform else vector
    matrix = np.zeros((nodes, nodes))
    tri_indices = np.triu_indices(nodes, k=1)
    matrix[tri_indices] = data
    matrix += matrix.T
    if zero_diagonal:
        np.fill_diagonal(matrix, 0)
    return matrix


def build_covariate_table(subject_info: pd.DataFrame, batches: Sequence[int]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": batches,
            "age": subject_info.age.to_list(),
            "sex": subject_info.sex.to_list(),
            "logwmh": subject_info.logwmh.to_list(),
        }
    )


def available_subject_ids(sc_raw_dir: Path) -> set[str]:
    return {path.name.split("_")[0] for path in sc_raw_dir.glob("*_sc-wei_*.npy")}


def derive_center_labels(subject_ids: Sequence[str]) -> list[int]:
    labels = []
    for sub_id in subject_ids:
        try:
            numeric_id = int(sub_id.split("-")[1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"Subject ID '{sub_id}' does not contain a numeric center code") from exc
        labels.append(1 if numeric_id <= 630 else 2)
    return labels


def run_neurocombat(table: pd.DataFrame, covariates: pd.DataFrame, reference_batch: str, mean_only: bool) -> np.ndarray:
    result = neuroCombat(
        dat=table,
        covars=covariates,
        batch_col="batch",
        categorical_cols=["sex"],
        ref_batch=reference_batch,
        mean_only=mean_only,
    )
    return result["data"]


def process_structural_connectomes(subjects: Sequence[str], subject_info: pd.DataFrame, center_labels: Sequence[int], config: CombatConfig) -> None:
    """Process and harmonize structural connectivity matrices for all subjects."""
    print("[Info] Processing SC", flush=True)
    covariates = build_covariate_table(subject_info, center_labels)

    for tract in config.tract_labels:
        print(f"[Info] Processing SC tract {tract}", flush=True)
        vectors = []
        for subject in subjects:
            sc_path = config.sc_raw_dir / f"{subject}_{tract}_sc-wei_{config.atlas_name}.npy"
            if not sc_path.exists():
                raise FileNotFoundError(f"Missing SC file: {sc_path}")
            raw_matrix = np.load(sc_path)
            vectors.append(upper_triangle_vector(raw_matrix, log_scale=True))

        multicenter_df = pd.DataFrame(np.column_stack(vectors), columns=subjects)
        pre_path = config.combat_dir / f"sc_pre_{tract}.csv"
        multicenter_df.to_csv(pre_path, index=False)

        harmonized = run_neurocombat(multicenter_df, covariates, config.reference_batch, mean_only=True)
        post_path = config.combat_dir / f"sc_data_combat_{tract}.csv"
        pd.DataFrame(harmonized).to_csv(post_path, header=False, index=False)

        for col_idx, subject in enumerate(subjects):
            harmonized_vec = harmonized[:, col_idx]
            reconstructed = vector_to_symmetric_matrix(harmonized_vec, config.atlas_nodes, exp_transform=True)

            mask_path = config.sc_raw_dir / f"{subject}_{tract}_sc-wei_{config.atlas_name}.npy"
            mask = np.load(mask_path)
            mask[mask > 0] = 1
            reconstructed *= mask

            thresholded = threshold_proportional(reconstructed, config.threshold)
            output_path = config.sc_dir / f"{subject}_{tract}_sc-wei_{config.atlas_name}.npy"
            np.save(output_path, thresholded)
            density = np.count_nonzero(thresholded) / (config.atlas_nodes**2) * 100
            print(f"[Info] {subject} ({tract}) density {density:.2f}%", flush=True)


def process_functional_connectomes(subjects: Sequence[str], subject_info: pd.DataFrame, center_labels: Sequence[int], config: CombatConfig) -> None:
    """Process and harmonize functional connectivity matrices for all subjects."""
    print("[Info] Processing FC", flush=True)
    covariates = build_covariate_table(subject_info, center_labels)

    vectors = []
    for subject in subjects:
        fc_path = config.fc_raw_dir / f"{subject}_fc-wb_{config.atlas_name}.npy"
        if not fc_path.exists():
            raise FileNotFoundError(f"Missing FC file: {fc_path}")
        raw_matrix = np.load(fc_path)
        vectors.append(upper_triangle_vector(raw_matrix, log_scale=False))

    multicenter_df = pd.DataFrame(np.column_stack(vectors), columns=subjects)
    multicenter_df.to_csv(config.combat_dir / "fc_pre.csv", index=False)

    harmonized = run_neurocombat(multicenter_df, covariates, config.reference_batch, mean_only=False)
    pd.DataFrame(harmonized).to_csv(config.combat_dir / "fc_data_combat.csv", header=False, index=False)

    for col_idx, subject in enumerate(subjects):
        harmonized_vec = harmonized[:, col_idx]
        fc_matrix = vector_to_symmetric_matrix(harmonized_vec, config.atlas_nodes, zero_diagonal=True)
        np.save(config.fc_dir / f"{subject}_fc-wb_{config.atlas_name}.npy", fc_matrix)


def main() -> None:
    """Main function to harmonize SC and FC matrices using neuroCombat."""
    config = CombatConfig()
    print(f"Run {__file__} for atlas {config.atlas_name}", flush=True)

    config.combat_dir.mkdir(parents=True, exist_ok=True)
    config.sc_dir.mkdir(parents=True, exist_ok=True)
    config.fc_dir.mkdir(parents=True, exist_ok=True)

    subject_info = pd.read_csv(config.subject_info_path)
    available_subjects = available_subject_ids(config.sc_raw_dir)
    subject_info = subject_info[subject_info.pid.isin(available_subjects)].copy()
    if subject_info.empty:
        raise RuntimeError("No subjects available after filtering for SC data.")

    subjects = subject_info.pid.to_list()
    center_labels = derive_center_labels(subjects)
    print(f"[Info] {len(subjects)} subjects retained for harmonization", flush=True)

    process_structural_connectomes(subjects, subject_info, center_labels, config)
    process_functional_connectomes(subjects, subject_info, center_labels, config)


if __name__ == "__main__":
    main()
