"""Aggregate LST lesion statistics and join with the pre-PSM subject list."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class WMHSummaryConfig:
    home_dir: Path = Path("/public/home/baishw/WMH_MCI")
    subject_list_name: str = "sub_list_beforePSM.csv"
    output_name: str = "sub_list_beforePSM_WMHstat.csv"

    @property
    def bids_dir(self) -> Path:
        return self.home_dir / "mica_bids"

    @property
    def subject_list_path(self) -> Path:
        return self.home_dir / "data" / self.subject_list_name

    @property
    def output_path(self) -> Path:
        return self.home_dir / "data" / self.output_name


def csdv_home_from_env() -> Path:
    env_path = os.getenv("MCI_CFC")
    if not env_path:
        raise EnvironmentError("MCI_CFC environment variable is not set.")
    return Path(env_path)


def load_subject_list(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Subject list not found: {path}")
    return pd.read_csv(path)


def lesion_stats_rows(bids_dir: Path) -> Iterable[pd.DataFrame]:
    for subject_dir in bids_dir.glob("sub-*"):
        stats_path = subject_dir / "wmh" / "lesion_stats.csv"
        if not stats_path.exists():
            continue
        data = pd.read_csv(stats_path)
        data["pid"] = subject_dir.name
        yield data


def collect_lesion_stats(bids_dir: Path) -> pd.DataFrame:
    frames = list(lesion_stats_rows(bids_dir))
    if not frames:
        raise RuntimeError("No lesion_stats.csv files were found.")
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    home_dir = csdv_home_from_env()
    config = WMHSummaryConfig(home_dir=home_dir)

    subject_list = load_subject_list(config.subject_list_path)
    lesion_stats = collect_lesion_stats(config.bids_dir)

    merged = subject_list.merge(lesion_stats, on="pid", how="left")
    merged = merged[merged["Lesion_Volume"] > 0].copy()
    merged["logwmh"] = np.log10(merged["Lesion_Volume"])

    config.output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(config.output_path, index=False)
    print(f"[Info] WMH stats saved to {config.output_path}")


if __name__ == "__main__":
    main()
