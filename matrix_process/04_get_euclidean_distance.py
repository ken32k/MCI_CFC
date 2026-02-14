"""Compute parcel-level Euclidean distance matrices from micapipe surfaces."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import nibabel as nib
import nilearn.surface
import numpy as np
from scipy.io import savemat
from scipy.spatial import distance, distance_matrix
from tqdm import tqdm


@dataclass(frozen=True)
class EuclideanConfig:
    """Centralize frequently reused paths and options."""

    home_dir: Path = Path("/public/home/baishw/WMH_MCI")
    atlas_name: str = "schaefer-200"
    output_subdir: str = "data/euc_distance"
    micapipe_subdir: str = "mica_out/micapipe_v0.2.0"
    freesurfer_subdir: str = "mica_out/freesurfer"
    metric_name: str = "euc"

    @property
    def output_dir(self) -> Path:
        return self.home_dir / self.output_subdir

    @property
    def micapipe_dir(self) -> Path:
        return self.home_dir / self.micapipe_subdir

    @property
    def freesurfer_dir(self) -> Path:
        return self.home_dir / self.freesurfer_subdir


def surface_coordinates(subject_id: str, hemisphere: str, config: EuclideanConfig) -> np.ndarray:
    """Load surface mesh coordinates for a given subject and hemisphere."""
    surf_path = (
        config.micapipe_dir
        / subject_id
        / "surf"
        / f"{subject_id}_hemi-{hemisphere}_surf-fsnative_label-midthickness.surf.gii"
    )
    if not surf_path.exists():
        raise FileNotFoundError(f"Missing surface mesh: {surf_path}")
    mesh = nilearn.surface.load_surf_mesh(str(surf_path))
    return mesh.coordinates


def read_parcellation_labels(subject_id: str, config: EuclideanConfig) -> np.ndarray:
    """Read parcellation labels from FreeSurfer annotation files for a given subject."""
    lh_path = config.freesurfer_dir / subject_id / "label" / f"lh.{config.atlas_name}_mics.annot"
    rh_path = config.freesurfer_dir / subject_id / "label" / f"rh.{config.atlas_name}_mics.annot"
    if not lh_path.exists() or not rh_path.exists():
        raise FileNotFoundError(f"Missing annotation for {subject_id}")

    labels_lh, ctab_lh, _ = nib.freesurfer.io.read_annot(str(lh_path), orig_ids=True)
    labels_rh, ctab_rh, _ = nib.freesurfer.io.read_annot(str(rh_path), orig_ids=True)

    parcel_ids = np.zeros(len(labels_lh) + len(labels_rh))
    for idx, label in enumerate(labels_lh):
        parcel_ids[idx] = np.where(ctab_lh[:, 4] == label)[0][0]
    offset = len(ctab_lh)
    for idx, label in enumerate(labels_rh):
        parcel_ids[idx + len(labels_lh)] = np.where(ctab_rh[:, 4] == label)[0][0] + offset
    return parcel_ids


def parcel_representatives(vertices: np.ndarray, parcels: np.ndarray, invalid_labels: Sequence[int]) -> np.ndarray:
    """Compute representative coordinates for each parcel using the geometric median approach."""
    valid_labels = [label for label in np.unique(parcels) if label not in invalid_labels]
    centers = np.zeros((len(valid_labels), 3))
    for idx, label in enumerate(valid_labels):
        parcel_vertices = vertices[parcels == label]
        if parcel_vertices.size == 0:
            raise ValueError(f"Parcel {label} contains no vertices")
        if parcel_vertices.shape[0] == 1:
            centers[idx] = parcel_vertices[0]
            continue
        condensed = distance.pdist(parcel_vertices, metric="euclidean")
        square = distance.squareform(condensed)
        min_index = np.argmin(square.sum(axis=1))
        centers[idx] = parcel_vertices[min_index]
    return centers


def compute_subject_distances(subject_id: str, config: EuclideanConfig) -> None:
    """Compute and save the Euclidean distance matrix for a given subject."""
    vertices_lh = surface_coordinates(subject_id, "L", config)
    vertices_rh = surface_coordinates(subject_id, "R", config)
    vertices = np.vstack((vertices_lh, vertices_rh))

    parcel_ids = read_parcellation_labels(subject_id, config)
    centers = parcel_representatives(vertices, parcel_ids, invalid_labels=(0, 200))

    distances = distance_matrix(centers, centers)
    output_path = config.output_dir / f"{subject_id}_{config.metric_name}_{config.atlas_name}.mat"
    savemat(output_path, {"ed": distances})


def subject_list(config: EuclideanConfig) -> list[str]:
    """Retrieve the list of subject IDs from the micapipe output directory."""
    if not config.micapipe_dir.exists():
        raise FileNotFoundError(f"Micapipe directory missing: {config.micapipe_dir}")
    return sorted(entry.name for entry in config.micapipe_dir.iterdir() if entry.is_dir())


def main() -> None:
    """Main function to compute Euclidean distance matrices for all subjects."""
    config = EuclideanConfig()
    config.output_dir.mkdir(parents=True, exist_ok=True)

    for subject in tqdm(subject_list(config)):
        if not subject.startswith("sub-"):
            continue
        output_file = config.output_dir / f"{subject}_{config.metric_name}_{config.atlas_name}.mat"
        if output_file.exists():
            print(f"[Info] {subject} already processed", flush=True)
            continue
        try:
            compute_subject_distances(subject, config)
            print(f"[Info] {subject} distances saved", flush=True)
        except Exception as exc:  # pylint: disable=broad-except
            print(f"[Error] {subject}: {exc}", flush=True)


if __name__ == "__main__":
    main()
