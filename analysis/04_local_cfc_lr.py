"""Compute network-level CFC linear regression metrics and visualizations."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import ana_utils

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class NetCFCConfig:
    project_home: Path = Path(ana_utils.PROJ_HOME)
    atlas: str = ana_utils.ATLAS
    cm_definitions: Sequence[Sequence] = tuple(ana_utils.CM_ARR[1:])
    tract_labels: Sequence[str] = tuple(ana_utils.TRACT_LABS)
    networks: dict = ana_utils.YEO7_DICT
    subject_list_name: str = "sub_list.csv"
    output_dirname: str = "results/lr_net_lm_cfc"
    n_jobs: int = -1

    @property
    def data_dir(self) -> Path:
        return self.project_home / "data"

    @property
    def cm_dir(self) -> Path:
        return self.data_dir / "CM"  # Communication measures

    @property
    def fc_dir(self) -> Path:
        return self.data_dir / "FCs"  # Functional connectivities

    @property
    def subject_list_path(self) -> Path:
        return self.data_dir / self.subject_list_name

    @property
    def output_dir(self) -> Path:
        return self.project_home / self.output_dirname


def load_subject_ids(config: NetCFCConfig) -> list[str]:
    """Load subject IDs from the subject list CSV file."""
    if not config.subject_list_path.exists():
        raise FileNotFoundError(f"Missing subject list: {config.subject_list_path}")
    return pd.read_csv(config.subject_list_path)["pid"].tolist()


def load_cm_block(
    subject_id: str,
    tract_label: str,
    cm_label: str,
    cm_index: int,
    network_mask: np.ndarray,
    config: NetCFCConfig,
) -> np.ndarray | None:
    """Load the CM block for a subject and network."""
    path = (
        config.cm_dir
        / cm_label
        / f"{subject_id}_{tract_label}_{cm_label}_{config.atlas}.npy"
    )
    if not path.exists():
        LOGGER.warning("Missing CM file %s", path)
        return None
    matrices = np.load(path)
    if cm_index >= matrices.shape[0]:
        LOGGER.warning("Index %d out of range for %s", cm_index, path)
        return None
    return matrices[cm_index][network_mask, :]


def load_fc_block(
    subject_id: str, network_mask: np.ndarray, config: NetCFCConfig
) -> np.ndarray | None:
    """Load the FC block for a subject and network."""
    path = config.fc_dir / f"{subject_id}_fc-wb_{config.atlas}.npy"
    if not path.exists():
        LOGGER.warning("Missing FC file %s", path)
        return None
    return np.load(path)[network_mask, :]


def compute_subject_network(
    subject_id: str, tract_label: str, network_name: str, config: NetCFCConfig
) -> dict | None:
    """Compute network-level CFC linear regression metrics for a subject and network."""
    # Get network mask
    mask = config.networks[network_name]

    # Load CM blocks
    cm_blocks = []
    for cm_label, cm_index in config.cm_definitions:
        block = load_cm_block(subject_id, tract_label, cm_label, cm_index, mask, config)
        if block is None:
            LOGGER.info(
                "Skipping %s (%s) missing %s", subject_id, network_name, cm_label
            )
            return None
        cm_blocks.append(block)

    # Load FC block
    fc_block = load_fc_block(subject_id, mask, config)
    if fc_block is None:
        LOGGER.info("Skipping %s (%s) missing FC", subject_id, network_name)
        return None
    # Perform linear regression
    try:
        dominance, rsq = ana_utils.sub_cfc_lr(fc_block, np.array(cm_blocks))
    except Exception as exc:  # pylint: disable=broad-except
        LOGGER.warning(
            "Regression failed for %s (%s, %s): %s",
            subject_id,
            tract_label,
            network_name,
            exc,
        )
        return None
    # Prepare result
    result = {
        "PID": subject_id,
        "Group": ana_utils.SUB_ID_DIAG_DICT.get(subject_id),
        "RSQ": rsq,
        "NETWORK": network_name,
    }
    # Add dominance values
    for idx, (cm_label, _) in enumerate(config.cm_definitions):
        result[cm_label] = dominance[idx]
    return result


def compute_subject(
    subject_id: str, tract_label: str, config: NetCFCConfig
) -> list[dict]:
    """Compute network-level CFC metrics for a subject across all networks."""
    entries = []
    for network_name in config.networks:
        entry = compute_subject_network(subject_id, tract_label, network_name, config)
        if entry:
            entries.append(entry)
    return entries


def run_for_tract(
    tract_label: str, subjects: Iterable[str], config: NetCFCConfig
) -> list[dict]:
    """Run network-level CFC linear regression for all subjects for a given tract."""
    LOGGER.info("Processing tract %s", tract_label)
    worker = delayed(compute_subject)
    results = Parallel(n_jobs=config.n_jobs)(
        worker(subject_id, tract_label, config) for subject_id in subjects
    )
    return [item for sublist in results for item in sublist]


def clear_directory(directory: Path) -> None:
    """Remove all files in the specified directory."""
    if not directory.exists():
        return
    for file_path in directory.glob("*"):
        if file_path.is_file():
            file_path.unlink()


def save_component_tables(
    cfc_table: pd.DataFrame, tract_label: str, config: NetCFCConfig
) -> None:
    """Save separate CSV files for each communication measure's dominance values."""
    for cm_label, cm_index in config.cm_definitions:
        cm_dom_tbl = cfc_table[["PID", "Group", "CM", "NETWORK"]].copy()
        cm_dom_tbl["IDX"] = cm_index
        cm_dom_tbl["RSQ"] = cfc_table[cm_label]
        cm_dom_tbl.to_csv(
            config.output_dir / f"{cm_label}-cfc_{tract_label}.csv", index=False
        )


def main() -> None:
    """Main function to compute network-level CFC linear regression metrics."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    LOGGER.info("Run %s", __file__)

    config = NetCFCConfig()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    clear_directory(config.output_dir)

    subjects = load_subject_ids(config)
    for tract_label in config.tract_labels:
        tract_results = run_for_tract(tract_label, subjects, config)
        if not tract_results:
            LOGGER.warning("No results for tract %s", tract_label)
            continue
        table = pd.json_normalize(tract_results)
        table.to_csv(config.output_dir / f"lrcms-cfc_{tract_label}.csv", index=False)
        save_component_tables(table, tract_label, config)
        LOGGER.info("Saved results for tract %s", tract_label)


if __name__ == "__main__":
    main()
