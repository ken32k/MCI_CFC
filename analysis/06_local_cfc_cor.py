"""Network-level coupling between communication metrics (CM) and FC."""

from dataclasses import dataclass, field
import logging
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import ana_utils


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(Path(__file__).stem)


@dataclass(frozen=True)
class NetCorrConfig:
    proj_home: Path = Path(ana_utils.PROJ_HOME)
    atlas: str = ana_utils.ATLAS
    tract_labels: tuple = tuple(ana_utils.TRACT_LABS)
    cm_defs: tuple = tuple(ana_utils.CM_ARR[1:])
    yeo_networks: Dict[str, np.ndarray] = field(
        default_factory=lambda: dict(ana_utils.YEO7_DICT)
    )
    n_jobs: int = ana_utils.N_JOBS

    @property
    def output_dir(self) -> Path:
        return self.proj_home / "results" / "corr_net_lm_cfc"

    @property
    def sub_list_path(self) -> Path:
        return self.proj_home / "data" / "sub_list.csv"


CONFIG = NetCorrConfig()


def _load_subject_ids(config: NetCorrConfig) -> List[str]:
    try:
        return pd.read_csv(config.sub_list_path, encoding="gbk")["pid"].to_list()
    except Exception:
        return pd.read_csv(config.sub_list_path, encoding="utf-8")["pid"].to_list()


def _cm_slices(cm_mat: np.ndarray) -> List[np.ndarray]:
    if cm_mat.ndim == 2:
        return [cm_mat]
    if cm_mat.ndim == 3:
        return [cm_mat[i] for i in range(cm_mat.shape[0])]
    raise ValueError(f"Unsupported CM dimensions: {cm_mat.shape}")


def _clear_old_files(directory: Path) -> None:
    for file in directory.glob("*"):
        if file.is_file():
            try:
                file.unlink()
            except Exception:
                logger.exception("Failed to delete %s", file)


def par_net_cfc_corr(cm_lab: str, sub_id: str, tract_labs: str) -> List[Dict[str, Any]]:
    """Compute network-level CFC correlation for a single subject."""
    proj = CONFIG.proj_home
    cm_path = proj / "data" / "cm" / cm_lab / f"{sub_id}_{tract_labs}_{cm_lab}_{CONFIG.atlas}.npy"
    fc_path = proj / "data" / "fcs" / f"{sub_id}_fc-wb_{CONFIG.atlas}.npy"

    if not cm_path.exists():
        logger.warning("Missing CM for %s at %s", sub_id, cm_path)
        return []
    if not fc_path.exists():
        logger.warning("Missing FC for %s at %s", sub_id, fc_path)
        return []

    try:
        cm_mat = np.load(cm_path)
    except Exception:
        logger.exception("Failed to load CM %s", cm_path)
        return []

    try:
        fc_mat = np.square(np.load(fc_path))
    except Exception:
        logger.exception("Failed to load FC %s", fc_path)
        return []

    sub_results: List[Dict[str, Any]] = []

    try:
        cm_iter = _cm_slices(cm_mat)
    except ValueError:
        logger.warning("Skipping %s due to invalid CM shape %s", sub_id, cm_mat.shape)
        return []

    for idx, cm_slice in enumerate(cm_iter):
        cm_view = np.array(cm_slice, copy=True)
        np.fill_diagonal(cm_view, np.nan)
        for net_name, indices in CONFIG.yeo_networks.items():
            try:
                r_values, p_values = ana_utils.sub_mat_correlatrion(
                    fc_mat[indices, :],
                    cm_view[indices, :],
                )
                rsq_val = float(np.nanmean(np.square(r_values)))
                p_val = float(np.nanmean(p_values))
            except Exception:
                logger.exception(
                    "Correlation failed (sub=%s, cm=%s, idx=%s, net=%s)",
                    sub_id,
                    cm_lab,
                    idx,
                    net_name,
                )
                rsq_val, p_val = np.nan, np.nan

            sub_results.append(
                {
                    "PID": sub_id,
                    "Group": ana_utils.SUB_ID_DIAG_DICT.get(sub_id, ""),
                    "CM": cm_lab,
                    "IDX": idx,
                    "RSQ": rsq_val,
                    "PVAL": p_val,
                    "NET": net_name,
                }
            )

    return sub_results


def get_cfc_corr(cm_lab: str, tract_labs: str) -> List[List[Dict[str, Any]]]:
    """Compute network-level CFC correlation for all subjects."""
    sub_ids = _load_subject_ids(CONFIG)
    return Parallel(n_jobs=CONFIG.n_jobs)(
        delayed(par_net_cfc_corr)(cm_lab, sub_id, tract_labs) for sub_id in sub_ids
    )


def main() -> None:
    logger.info("Run %s", __file__)

    output_dir = CONFIG.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    _clear_old_files(output_dir)

    logger.info("Yeo networks: %s", list(CONFIG.yeo_networks.keys()))

    merged_frames: List[pd.DataFrame] = []

    for cm_lab, _, _ in CONFIG.cm_defs:
        logger.info("Processing CM: %s", cm_lab)
        for tract_label in CONFIG.tract_labels:
            cfc_results = get_cfc_corr(cm_lab, tract_label)
            flattened = [item for sublist in cfc_results for item in sublist if item]

            if not flattened:
                logger.info("No results for cm=%s tract=%s", cm_lab, tract_label)
                continue

            cfc_results_tbl = pd.json_normalize(flattened)
            out_path = output_dir / f"{cm_lab}-cfc_{tract_label}.csv"
            cfc_results_tbl.to_csv(out_path, index=False, encoding="utf-8")
            merged_frames.append(cfc_results_tbl)
            logger.info("Saved %s", out_path)

    if merged_frames:
        merged_tbl = pd.concat(merged_frames, axis=0, ignore_index=True)
        merged_path = output_dir / "merged-cfc.csv"
        merged_tbl.to_csv(merged_path, index=False, encoding="utf-8")
        logger.info("Merged results saved to %s", merged_path)
    else:
        logger.warning("No merged results generated")

    try:
        import plot_net_cfc  # noqa: F401

        logger.info("plot_net_cfc module imported for downstream plotting")
    except ImportError:
        logger.warning("plot_net_cfc module not found; skipping plotting step")


if __name__ == "__main__":
    main()