"""
Global univariate coupling: correlate communication metrics (CM) with whole-brain FC.

Improvements in this refactor:
- safer Path handling (`pathlib.Path`)
- structured logging instead of prints
- robust file existence checks and graceful skipping
- support for static (2D) and dynamic (3D) CM arrays
- `if __name__ == '__main__'` guard so module can be imported
"""

from pathlib import Path
import logging
from typing import List, Dict, Any

import ana_utils
import numpy as np
import pandas as pd
from joblib import Parallel, delayed


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(Path(__file__).stem)


def par_global_cfc_corr(cm_lab: str, sub_id: str, tract_labs: str) -> List[Dict[str, Any]]:
    """Compute per-subject correlation(s) between CM (may be dynamic) and FC.

    Returns a list of result dicts (one per dynamic index, or one if static).
    Missing files are skipped and produce an empty list.
    """
    proj = Path(ana_utils.PROJ_HOME)
    cm_path = proj / "data" / "CMs" / cm_lab / f"{sub_id}_{tract_labs}_{cm_lab}_{ana_utils.ATLAS}.npy"
    fc_path = proj / "data" / "FCs" / f"{sub_id}_fc-wb_{ana_utils.ATLAS}.npy"

    if not cm_path.exists():
        logger.warning("Missing CM for %s: %s", sub_id, cm_path)
        return []
    if not fc_path.exists():
        logger.warning("Missing FC for %s: %s", sub_id, fc_path)
        return []

    try:
        cm_mat = np.load(cm_path)
    except Exception as e:
        logger.exception("Failed to load CM for %s: %s", sub_id, e)
        return []

    try:
        fc_mat = np.load(fc_path)
    except Exception as e:
        logger.exception("Failed to load FC for %s: %s", sub_id, e)
        return []

    sub_cfc_result: List[Dict[str, Any]] = []

    # support static CM (2D) or dynamic CM (3D: time x parcels x parcels)
    if cm_mat.ndim == 2:
        cm_iters = [cm_mat]
    elif cm_mat.ndim == 3:
        cm_iters = [cm_mat[i] for i in range(cm_mat.shape[0])]
    else:
        logger.warning("Unexpected CM shape for %s: %s", sub_id, cm_mat.shape)
        return []

    for j, cm_sub_mat in enumerate(cm_iters):
        try:
            r_values, p_values = ana_utils.sub_mat_correlatrion(fc_mat, cm_sub_mat)
            rsq_values = np.power(r_values, 2)
            rsq_val, p_val = float(np.nanmean(rsq_values)), float(np.nanmean(p_values))
        except Exception:
            logger.exception("Correlation failed for %s (cm=%s, idx=%s)", sub_id, cm_lab, j)
            rsq_val, p_val = np.nan, np.nan

        sub_cfc_result.append(
            {
                "pid": sub_id,
                "Group": ana_utils.SUB_ID_DIAG_DICT.get(sub_id, ""),
                "cm": cm_lab,
                "idx": j,
                "rsq": rsq_val,
                "pval": p_val,
            }
        )

    return sub_cfc_result


def get_cfc_corr(cm_lab: str, tract_labs: str) -> List[List[Dict[str, Any]]]:
    """Run correlations across subjects (parallelized) and return nested result lists."""

    sub_list_path = Path(ana_utils.PROJ_HOME) / "data" / "sub_list.csv"
    # keep original encoding used in project, fallback to utf-8
    try:
        sub_id_list = pd.read_csv(sub_list_path, encoding="gbk")["pid"].to_list()
    except Exception:
        sub_id_list = pd.read_csv(sub_list_path, encoding="utf-8")["pid"].to_list()

    # Main parallel loop
    results = Parallel(n_jobs=ana_utils.N_JOBS)(
        delayed(par_global_cfc_corr)(cm_lab, sub_id, tract_labs) for sub_id in sub_id_list
    )

    return results


def _clear_dir_files(dirpath: Path) -> None:
    for p in dirpath.iterdir():
        if p.is_file():
            try:
                p.unlink()
            except Exception:
                logger.exception("Failed to delete file: %s", p)


def main() -> None:
    logger.info("Run %s", __file__)

    output_dir = Path(ana_utils.PROJ_HOME) / "results" / "corr_glob_lm_cfc"
    output_dir.mkdir(parents=True, exist_ok=True)
    _clear_dir_files(output_dir)

    cm_arr = ana_utils.CM_ARR[1:]

    for cm in cm_arr:
        cm_lab = cm[0]
        for TRACT_LABS in ana_utils.TRACT_LABS:
            logger.info("processing cm=%s tract=%s", cm_lab, TRACT_LABS)
            cfc_results = get_cfc_corr(cm_lab, TRACT_LABS)

            # flatten nested lists and skip empty results
            flattened_list: List[Dict[str, Any]] = [item for sublist in cfc_results for item in sublist if item]

            if not flattened_list:
                logger.info("No results for cm=%s tract=%s (skipping)", cm_lab, TRACT_LABS)
                continue

            cfc_results_tbl = pd.json_normalize(flattened_list)

            out_path = output_dir / f"{cm_lab}-cfc_{TRACT_LABS}.csv"
            try:
                cfc_results_tbl.to_csv(out_path, index=False, encoding="utf-8")
                logger.info("Wrote results to %s", out_path)
            except Exception:
                logger.exception("Failed to write CSV to %s", out_path)


if __name__ == "__main__":
    main()

