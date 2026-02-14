"""Summarize FC/SC matrix sums and generate QC plots."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ana_utils

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class MatrixSummaryConfig:
    project_home: Path = Path(ana_utils.PROJ_HOME)
    atlas: str = ana_utils.ATLAS
    tract_labels: tuple = tuple(ana_utils.TRACT_LABS)
    subject_ids: tuple = tuple(ana_utils.SUB_ID_LIST)
    output_dirname: str = "results/mat_stat"

    @property
    def data_dir(self) -> Path:
        return self.project_home / "data"

    @property
    def fc_dir(self) -> Path:
        return self.data_dir / "fcs"

    @property
    def cm_dir(self) -> Path:
        return self.data_dir / "cm" / "sc-wei"

    @property
    def output_dir(self) -> Path:
        return self.project_home / self.output_dirname


def clear_directory(directory: Path) -> None:
    """Remove all files in the specified directory."""
    if not directory.exists():
        return
    for file_path in directory.glob("*"):
        if file_path.is_file():
            file_path.unlink()


def load_matrix_sum(path: Path) -> float:
    """Load a connectivity matrix and return the sum of its elements."""
    if not path.exists():
        raise FileNotFoundError(f"Missing matrix file: {path}")
    matrix = np.load(path)
    if matrix.ndim == 3:
        matrix = matrix[0]
    return float(matrix.sum())


def collect_matrix_summaries(config: MatrixSummaryConfig) -> pd.DataFrame:
    """Collect sum statistics for FC and SC matrices across subjects."""
    records: List[Dict] = []
    for subject_id in config.subject_ids:
        group = ana_utils.SUB_ID_DIAG_DICT.get(subject_id)
        fc_path = config.fc_dir / f"{subject_id}_fc-wb_{config.atlas}.npy"
        fc_sum = load_matrix_sum(fc_path)
        records.append(
            {
                "pid": subject_id,
                "Group": group,
                "cm": "fc",
                "tg_lab": "wb",
                "value": fc_sum,
            }
        )
        # SC matrices
        for tract_label in config.tract_labels:
            cm_path = config.cm_dir / f"{subject_id}_{tract_label}_sc-wei_{config.atlas}.npy"
            cm_sum = load_matrix_sum(cm_path)
            records.append(
                {
                    "pid": subject_id,
                    "Group": group,
                    "cm": "sc-wei",
                    "tg_lab": tract_label,
                    "value": cm_sum,
                }
            )
    return pd.DataFrame.from_records(records)


def plot_summary(data: pd.DataFrame, title_suffix: str, output_path: Path) -> None:
    """Generate and save a swarm plot for the given matrix summary data."""
    if data.empty:
        LOGGER.warning("No data available for %s", title_suffix)
        return
    fig, ax = plt.subplots(figsize=(4, 4))
    ylim = data["value"].max() * 1.5 if data["value"].any() else 1.0
    ana_utils.plot_ax_swarmplot(data, ax=ax, i=0, y="value", ylim=ylim)
    ax.set_title(title_suffix)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Main function to summarize matrices and generate QC plots."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    LOGGER.info("Run %s", __file__)

    config = MatrixSummaryConfig()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    clear_directory(config.output_dir)

    summary_table = collect_matrix_summaries(config)
    csv_path = config.output_dir / "sub_mat_sum.csv"
    summary_table.to_csv(csv_path, index=False)
    LOGGER.info("Saved matrix summary to %s", csv_path)

    plot_summary(
        summary_table[summary_table.cm == "fc"],
        "FC",
        config.output_dir / "mat-fc.png",
    )

    for tract_label in config.tract_labels:
        subset = summary_table[(summary_table.cm == "sc-wei") & (summary_table.tg_lab == tract_label)]
        plot_summary(subset, f"SC {tract_label}", config.output_dir / f"mat-sc-{tract_label}.png")


if __name__ == "__main__":
    main()
