#!/bin/bash
set -euo pipefail

# Basic metadata for reference
# WMH segmentation using LST on bias-corrected T1w/FLAIR volumes

subject_id=${1:?"Subject ID (e.g., 001) is required"}

: "${MCI_CFC_DIR:?MCI_CFC_DIR must be defined}"
export PROJ_HOME="${MCI_CFC_DIR}"

if ! command -v lst >/dev/null; then
    # shellcheck disable=SC1091
    source ~/apps/lst/bin/activate
fi
export PATH="${PATH}:~/apps"

bids_dir="${PROJ_HOME}/mica_bids"
subject_dir="${bids_dir}/sub-${subject_id}"
anat_dir="${subject_dir}/anat"
wmh_dir="${subject_dir}/wmh"

t1_path="${anat_dir}/sub-${subject_id}_T1w.nii.gz"
t1_n4_path="${anat_dir}/sub-${subject_id}_T1w_n4.nii.gz"
flair_path="${anat_dir}/sub-${subject_id}_FLAIR.nii.gz"
flair_n4_path="${anat_dir}/sub-${subject_id}_FLAIR_n4.nii.gz"

if compgen -G "${wmh_dir}"/*lst*.nii.gz >/dev/null; then
    echo "[Info] Skip subject ${subject_id}: WMH output already exists"
    exit 0
fi

echo "[Info] Processing subject ${subject_id}"
mkdir -p "${wmh_dir}"

N4BiasFieldCorrection -d 3 -i "${flair_path}" -o "${flair_n4_path}"
N4BiasFieldCorrection -d 3 -i "${t1_path}" -o "${t1_n4_path}"

lst \
    --t1 "${t1_n4_path}" \
    --flair "${flair_n4_path}" \
    --output "${wmh_dir}" \
    --device cpu \
    >/dev/null

rm -f "${t1_n4_path}" "${flair_n4_path}"
echo "[Info] Subject ${subject_id} completed"
