"""Compute global CFC linear regression metrics for each subject and tract."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import ana_utils

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class GlobalCFCConfig:
    project_home: Path = Path(ana_utils.PROJ_HOME)
    atlas: str = ana_utils.ATLAS
    cm_definitions: Sequence[Sequence] = tuple(ana_utils.CM_ARR[1:])
    tract_labels: Sequence[str] = tuple(ana_utils.TRACT_LABS)
    subject_list_name: str = "sub_list.csv"
    output_dirname: str = "results/1_global_CFC"
    n_jobs: int = -1

    @property
    def data_dir(self) -> Path:
        return self.project_home / "data"

    @property
    def cm_dir(self) -> Path:
        return self.data_dir / "CMs" # Communication measures

    @property
    def fc_dir(self) -> Path:
        return self.data_dir / "FCs" # Functional connectivities

    @property
    def subject_list_path(self) -> Path:
        return self.data_dir / self.subject_list_name

    @property
    def output_dir(self) -> Path:
        return self.project_home / self.output_dirname


def load_subject_ids(config: GlobalCFCConfig) -> list[str]:
    """Load subject IDs from the subject list CSV file."""
    if not config.subject_list_path.exists():
        raise FileNotFoundError(f"Subject list missing: {config.subject_list_path}")
    return pd.read_csv(config.subject_list_path)["pid"].tolist()


def load_cm_vector(subject_id: str, tract_label: str, cm_label: str, cm_index: int, config: GlobalCFCConfig) -> Optional[np.ndarray]:
    """Load the CM vector for a subject."""
    path = config.cm_dir / cm_label / f"{subject_id}_{tract_label}_{cm_label}_{config.atlas}.npy"
    if not path.exists():
        LOGGER.warning("Missing CM file %s", path)
        return None
    matrices = np.load(path)
    if cm_index >= matrices.shape[0]:
        LOGGER.warning("Index %d out of bounds for %s", cm_index, path)
        return None
    return matrices[cm_index][ana_utils.TRIU_IDX]


def load_fc_vector(subject_id: str, config: GlobalCFCConfig) -> Optional[np.ndarray]:
    """Load the FC vector for a subject."""
    path = config.fc_dir / f"{subject_id}_fc-wb_{config.atlas}.npy"
    if not path.exists():
        LOGGER.warning("Missing FC file %s", path)
        return None
    return np.load(path)[ana_utils.TRIU_IDX]


def compute_subject_metrics(subject_id: str, tract_label: str, config: GlobalCFCConfig) -> Optional[dict]:
    """Compute CFC linear regression metrics for a subject and tract."""
    # Load CM vectors
    cm_vectors = []
    for cm_label, cm_index in config.cm_definitions:
        cm_vec = load_cm_vector(subject_id, tract_label, cm_label, cm_index, config)
        if cm_vec is None:
            LOGGER.info("Skipping %s: missing %s", subject_id, cm_label)
            return None
        cm_vectors.append(cm_vec)

    # Load FC vector
    fc_vector = load_fc_vector(subject_id, config)
    if fc_vector is None:
        LOGGER.info("Skipping %s: missing FC", subject_id)
        return None

    # Linear regression
    try:
        dominance, rsq = ana_utils.sub_cfc_lr(fc_vector, np.array(cm_vectors))
    except Exception as exc:  # pylint: disable=broad-except
        LOGGER.warning("LR model failed for %s (%s): %s", subject_id, tract_label, exc)
        return None

    # Prepare result
    result = {
        "PID": subject_id,
        "Group": ana_utils.SUB_ID_DIAG_DICT.get(subject_id),
        "RSQ": rsq,
    }
    # Add dominance values
    for idx, (cm_label, _) in enumerate(config.cm_definitions):
        result[cm_label] = dominance[idx]
    return result


def run_for_tract(tract_label: str, subjects: Iterable[str], config: GlobalCFCConfig) -> list[dict]:
    """Run global CFC linear regression for all subjects for a given tract."""
    LOGGER.info("Processing tract %s", tract_label)
    processor = delayed(compute_subject_metrics)
    results = Parallel(n_jobs=config.n_jobs)(
        processor(subject_id, tract_label, config) for subject_id in subjects
    )
    return [res for res in results if res]


def main() -> None:
    """Main function to compute global CFC linear regression metrics."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config = GlobalCFCConfig()
    config.output_dir.mkdir(parents=True, exist_ok=True)

    subject_ids = load_subject_ids(config)
    for tract_label in config.tract_labels:
        tract_results = run_for_tract(tract_label, subject_ids, config)
        if not tract_results:
            LOGGER.warning("No valid results for tract %s", tract_label)
            continue
        output_path = config.output_dir / f"_lr-global_cfc-{tract_label}.csv"
        pd.json_normalize(tract_results).to_csv(output_path, index=False)
        LOGGER.info("Saved %s", output_path)


if __name__ == "__main__":
    main()