"""Select subjects with complete CM/FC data and export metadata before propensity score matching (PSM).

Filtering before PSM ensures that only subjects with complete data and within specified demographic
criteria are included in subsequent analyses. The script outputs a CSV file containing relevant
subject information for further processing.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

import ana_utils

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class SubjectListConfig:
    """Container for frequently reused paths and filtering parameters."""

    project_home: Path = Path(ana_utils.PROJ_HOME)
    atlas_name: str = ana_utils.ATLAS
    cm_definitions: Sequence[Sequence] = tuple(ana_utils.CM_ARR[1:])
    age_range: tuple[int, int] = (50, 80)
    n_nodes: int = 200
    tract_label: str = "wb"
    

    @property
    def data_dir(self) -> Path:
        return self.project_home / "data"

    @property
    def mica_dir(self) -> Path:
        return self.project_home / "mica_out" / "micapipe_v0.2.0"

    @property
    def bids_dir(self) -> Path:
        return self.project_home / "mica_bids"

    @property
    def cm_dir(self) -> Path:
        return self.data_dir / "cm"

    @property
    def fc_dir(self) -> Path:
        return self.data_dir / "fcs_raw"

    @property
    def subject_info(self) -> pd.DataFrame:
        return ana_utils.SUB_INFO_TBL.copy()

    @property
    def output_csv(self) -> Path:
        return self.data_dir / "sub_list_beforePSM.csv"


def matrix_is_invalid(matrix: np.ndarray) -> bool:
    """Return True when a matrix contains NaN/Inf values or is entirely zeros."""

    if not np.isfinite(matrix).all():
        return True
    if not np.any(matrix):
        return True
    return False


def load_cm_matrix(subject_id: str, cm_label: str, index: int, config: SubjectListConfig) -> Optional[np.ndarray]:
    path = config.cm_dir / cm_label / f"{subject_id}_{config.tract_label}_{cm_label}_{config.atlas_name}.npy"
    if not path.exists():
        return None

    matrices = np.load(path)
    if index >= matrices.shape[0]:
        return None

    matrix = matrices[index].copy()
    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(matrix, 0)

    if matrix_is_invalid(matrix):
        LOGGER.debug("[Skip] %s %s invalid matrix", subject_id, cm_label)
        return None
    return matrix


def has_fcs(subject_id: str, config: SubjectListConfig) -> bool:
    return (config.fc_dir / f"{subject_id}_fc-wb_{config.atlas_name}.npy").exists()


def count_available_subjects(subjects: Iterable[str], directory: Path, relative_pattern: str) -> int:
    """Count subjects for which a specific relative file exists under `directory/sub`.

    relative_pattern should include ``{sub}`` placeholder that will be replaced with
    the subject identifier (e.g., "anat/{sub}_FLAIR.nii.gz").
    """

    count = 0
    for subject_id in subjects:
        relative_path = relative_pattern.format(sub=subject_id)
        if (directory / subject_id / relative_path).exists():
            count += 1
    return count


def subject_has_all_cm(subject_id: str, config: SubjectListConfig) -> tuple[bool, list[str]]:
    """Check if a subject has all required connectivity matrices."""
    missing: list[str] = []
    for cm_label, index in config.cm_definitions:
        matrix = load_cm_matrix(subject_id, cm_label, index, config)
        if matrix is None:
            missing.append(cm_label)
            
    # Return Boolean indicating if all CMs are present, and list of missing CMs
    return (len(missing) == 0, missing)


def age_within_range(subject_id: str, config: SubjectListConfig) -> bool:
    """Age filter for a subject."""
    age = ana_utils.SUB_ID_AGE_DICT.get(subject_id)
    if age is None:
        return False
    lower, upper = config.age_range
    return lower <= age <= upper


def filter_subjects(config: SubjectListConfig) -> list[str]:
    """Filter subjects based on data completeness and demographic criteria."""
    selected: list[str] = []
    for subject in ana_utils.SUB_ID_LIST:
        # Check for required connectivity matrices
        has_all, missing = subject_has_all_cm(subject, config)
        if not has_all:
            LOGGER.info("[Exclude] %s missing CMs: %s", subject, ",".join(missing))
            continue
        if not has_fcs(subject, config):
            LOGGER.info("[Exclude] %s missing FC matrix", subject)
            continue
        if not age_within_range(subject, config):
            LOGGER.info("[Exclude] %s age outside range", subject)
            continue
        selected.append(subject)
        LOGGER.info("[Include] %s", subject)
    return selected


def build_subject_table(subject_ids: Iterable[str], config: SubjectListConfig) -> pd.DataFrame:
    """Build a subject metadata table for the selected subjects."""
    table = pd.DataFrame({"pid": subject_ids})
    table["Group"] = table["pid"].map(ana_utils.SUB_ID_DIAG_DICT)
    table["age"] = table["pid"].map(ana_utils.SUB_ID_AGE_DICT)
    table["sex"] = table["pid"].map(ana_utils.SUB_ID_GEN_DICT)
    table = table.merge(config.subject_info[["pid", "logwmh"]], on="pid", how="left")
    return table


def main() -> None:
    """Main function to generate subject list before PSM."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    LOGGER.info("Run %s", __file__)

    config = SubjectListConfig()
    config.output_csv.unlink(missing_ok=True)

    bids_count = count_available_subjects(
        ana_utils.SUB_ID_LIST,
        config.bids_dir,
        "anat/{sub}_FLAIR.nii.gz",
    )
    mica_count = count_available_subjects(
        ana_utils.SUB_ID_LIST,
        config.mica_dir,
        f"dwi/connectomes-wmh-2M/{{sub}}_space-dwi_atlas-schaefer-{config.n_nodes}_desc-iFOD2-2M-SIFT2_full-connectome.shape.gii",
    )
    LOGGER.info("Subjects with FLAIR: %d", bids_count)
    LOGGER.info("Subjects with connectomes: %d", mica_count)

    filtered_subjects = filter_subjects(config)
    LOGGER.info("Selected %d subjects", len(filtered_subjects))

    subject_table = build_subject_table(filtered_subjects, config)
    subject_table.to_csv(config.output_csv, index=False)
    LOGGER.info("Subject table saved to %s", config.output_csv)


if __name__ == "__main__":
    main()
