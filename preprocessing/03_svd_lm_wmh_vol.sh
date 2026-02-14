#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-07-01 13:54:48
# @LastEditors: Ken32g ken32k@163.com
# @LastEditTime: 2024-12-07 10:52:30
# @FilePath: //svd_lm_wmh_vol.sh
# @Description: Get the mean values of DTI derived metrics and T1FLAIR ratio for mediation analysis
###

export PROJ_HOME=$MCI_CFC_DIR
# Make directory
mkdir -p ${PROJ_HOME}/results/lms_stat/

# Get subject list
sub_list_path=${PROJ_HOME}/data/sub_list.csv

while IFS= read -r line; do
    array_value=$(echo "$line" | awk -F',' '{print substr($3, 1, 7)}')
    sub_list+=("$array_value")
done < <(tail -n +2 "$sub_list_path")

# Get mean value of DTI and T1FLAIR
for sub_dir in ${sub_list[@]}; do

    sub_id=$(basename "${sub_dir}")
    # echo $sub_id

    # WMH DWI metrics
    dwi_dir=${PROJ_HOME}/mica_out/micapipe_v0.2.0/${sub_id}/dwi
    dwi_dti="${proc_dwi}/${sub_id}_space-dwi_model-DTI.nii.gz"
    ad=${dwi_dir}/${sub_id}_space-dwi_model-DTI_map-AD.nii.gz
    fa=${dwi_dir}/${sub_id}_space-dwi_model-DTI_map-FA.nii.gz
    rd=${dwi_dir}/${sub_id}_space-dwi_model-DTI_map-RD.nii.gz
    adc=${dwi_dir}/${sub_id}_space-dwi_model-DTI_map-ADC.nii.gz
    fivett=${dwi_dir}/${sub_id}_space-dwi_desc-5tt.nii.gz

    # Check 5tt file exists
    if [ ! -f "$fivett" ]; then
        echo "[Error] ... ${sub_dir}: 5tt file not found, pass."
        continue
    fi

    # Check dwi metrics file exists
    if [ ! -f "$rd" ]; then
        if [ -f "$dwi_dti" ]; then
            tensor2metric -nthreads 1 -fa "$fa" -ad "$ad" -rd "$rd" -adc "$adc" "$dwi_dti" -force
        else
            echo "[Error] ... ${sub_dir}: DTI file not found, pass."
            continue
        fi
    fi

    # NAWM
    wmparc_mgz=${PROJ_HOME}/mica_out/freesurfer/${sub_id}/mri/wmparc.mgz
    if [ ! -f "$wmparc_mgz" ]; then
        echo "[Error] ... ${sub_dir}: wmparc file not found, pass."
        continue
    fi
    dst_nifti=${PROJ_HOME}/mica_out/micapipe_v0.2.0/${sub_id}/anat/${sub_id}_cerebral_WM.nii.gz
    wmh_nifti=${PROJ_HOME}/mica_out/micapipe_v0.2.0/${sub_id}/anat/${sub_id}_space-5tt_wmh.nii.gz
    nawm_bin=${PROJ_HOME}/mica_out/micapipe_v0.2.0/${sub_id}/anat/${sub_id}_cerebral_WM_bin.nii.gz
    if [ ! -f "$nawm_bin" ]; then
        mrconvert $wmparc_mgz $dst_nifti &>/dev/null
        fslmaths $dst_nifti -thr 3000 -bin -sub $wmh_nifti $nawm_bin &>/dev/null
    fi

    fslstats ${nawm_bin} -V | awk -v X=${sub_id} '{print X",",$2}' >>${PROJ_HOME}/results/lms_stat/nawm_volume_raw.csv

    # # T1FLAIR
    # tfr_raw=${PROJ_HOME}/mica_bids/T1T2FLAIRMTR_ratio/${sub_id}/ses-${sub_id}/${sub_id}_ratio-T1FLAIR_calib-nonlin.nii.gz
    # tfr_anat=${PROJ_HOME}/mica_out/micapipe_v0.2.0/${sub_id}/anat/${sub_id}_ratio-T1FLAIR_calib-nonlin.nii.gz
    # tfr_mask=${PROJ_HOME}/mica_out/micapipe_v0.2.0/${sub_id}/anat/${sub_id}_ratio-T1FLAIR_mask.nii.gz
    # if [ -f $tfr_raw ]; then
    #     cp $tfr_raw $tfr_anat &>/dev/null
    #     antsApplyTransforms -d 3 -i $wmh_nifti -r $tfr_anat -o $tfr_mask -n NearestNeighbor -v 0
    # else
    #     echo [Error] ... No T1FLAIR file exist
    # fi

    # Stat DWI metric
    cd ${dwi_dir}
    fslsplit ${fivett} 5tt-split-
    for five_lab in 2 4; do
        wmh_mask="${dwi_dir}/5tt-split-000${five_lab}.nii.gz"
        csv_file="${PROJ_HOME}/results/lms_stat/${five_lab}_wmh_dwi_metric.csv"

        fslstats "${ad}" -k "${wmh_mask}" -M | awk -v X="ad" -v Y="${sub_id}" '{printf "%s,%s,%.12f\n", Y, X, $1}' >>"$csv_file"
        fslstats "${fa}" -k "${wmh_mask}" -M | awk -v X="fa" -v Y="${sub_id}" '{printf "%s,%s,%.12f\n", Y, X, $1}' >>"$csv_file"
        fslstats "${rd}" -k "${wmh_mask}" -M | awk -v X="rd" -v Y="${sub_id}" '{printf "%s,%s,%.12f\n", Y, X, $1}' >>"$csv_file"
        fslstats "${adc}" -k "${wmh_mask}" -M | awk -v X="adc" -v Y="${sub_id}" '{printf "%s,%s,%.12f\n", Y, X, $1}' >>"$csv_file"
        if [ -f $tfr_anat ]; then
            fslstats "${tfr_anat}" -k "${tfr_mask}" -M | awk -v X="tfr" -v Y="${sub_id}" '{printf "%s,%s,%.12f\n", Y, X, $1}' >>"$csv_file"
            fslstats "${tfr_anat}" -k "${tfr_mask}" -R | awk -v X="tfrmin" -v Y="${sub_id}" '{printf "%s,%s,%.12f\n", Y, X, $1}' >>"$csv_file"
        fi
    done

    # rm tmp files
    rm ${dwi_dir}/5tt-split-000*.nii.gz &>/dev/null
    echo "[Info] ... ${sub_dir}: Done."
done
