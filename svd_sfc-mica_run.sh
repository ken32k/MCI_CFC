#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-07-01 14:15:30
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-12-12 12:38:54
# @FilePath: /csvd-sfc/svd_sfc-mica_run.sh
# @Description: Structure DTI and fMRI preprocessing with MICA
###

# for hpc
# set -e

# Check the number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <patient_id> <number_of_threads> <micapipe_path>"
    echo "Example: $0 001 4 /path/to/micapipe"
    exit 1
fi

# Arguments
sub=$1           # patient id, e.g., 001, 101
nthr=$2          # number of threads
micapipe_path=$3 # path to MICAPIPE, 1 for single-shell, 2 for multi-shell

# Set micapipe version
if [ -z "$micapipe_path" ]; then
    echo "micapipe_path is not set."
    exit 1
fi

if [ "$micapipe_path" -eq 2 ]; then
    # Multishell version
    MICAPIPE="$HOME/micapipe/micapipe-master"
elif [ "$micapipe_path" -eq 1 ]; then
    # Single shell version
    MICAPIPE="$HOME/micapipe/micapipe-single-shell"
else
    echo "Invalid value for micapipe_path: $micapipe_path"
    exit 1 # 处理无效值
fi
# Source
source $HOME/.bashrc
# export PROJ_HOME=/public/home/baishw/WMH_MCI
export PROJ_HOME=$CSVD_DIR

# Set PATH and MICAPIPE
PATH="$PATH:${MICAPIPE}:${MICAPIPE}/functions"
export PATH
export MICAPIPE

export FREESURFER_HOME=/public/software/apps/freesurfer_infant/freesurfer7.3.2/freesurfer/7.3.2-1
source /public/software/apps/freesurfer_infant/freesurfer7.3.2/freesurfer/7.3.2-1/SetUpFreeSurfer.sh
export SUBJECTS_DIR=${PROJ_HOME}/csvd_mica_out/fastsurfer
export PROC=container_micapipe-v0.2.2

# set subject path
bids=${PROJ_HOME}/csvd_mica_bids
out=${PROJ_HOME}/csvd_mica_out
tmp=${PROJ_HOME}/csvd_mica_temp

# subject path
out_sub_dir=${out}/micapipe_v0.2.0/sub-${sub}
fs_out_dir=${out}/freesurfer/sub-${sub}

# fs_licence
fslc=${HOME}/micapipe/FreeSurfer_License_71515.txt
# tract number
tract_num='2M'

atlas='schaefer-200'
# atlas2='schaefer-200'

# start subject
echo "################################"

start_time=$(date +%s)
formatted_start_time=$(date -d @$start_time +"%Y-%m-%d %H:%M:%S")
lopuu=$(date -d @$start_time +"%Y-%m-%d")
echo "${formatted_start_time}---subject: ${sub}"
printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# Check bids files
bidst1="${bids}/sub-${sub}/anat/sub-${sub}_T1w.nii.gz"
bidsflair="${bids}/sub-${sub}/anat/sub-${sub}_FLAIR.nii.gz"

# bidswmh="${bids}/sub-${sub}/anat/sub-${sub}_WMH_mask.nii.gz"
bidswmh="${bids}/sub-${sub}/wmh/space-flair_seg-lst.nii.gz"
bidsdwi="${bids}/sub-${sub}/dwi/sub-${sub}_dir-AP_dwi.nii.gz"
bidsbold="${bids}/sub-${sub}/func/sub-${sub}_dir-AP_rest_bold.nii.gz"

# if any file is lost, exit
if [ -f $bidst1 ] && [ -f $bidsflair ] && [ -f $bidswmh ] && [ -f $bidsdwi ]; then
    echo [Info] ... bids files checked
else
    echo [Error] ... bids files loss
    exit
fi

# ---------------------------------------------------------
# MICAPIPE-proc_struct_surf
log_proc_struct=${out_sub_dir}/logs/proc_structural*
if [ ! -f "${log_proc_struct}" ]; then
    echo [Info] ... MICA-proc_struct_surf -no-output
    micapipe -sub ${sub} -out ${out} -bids ${bids} -tmpDir ${tmp} -threads ${nthr} -fs_licence ${fslc} \
        -proc_structural -mf 45 -freesurfer \
        -proc_surf &>/dev/null
    # &>/dev/null
    echo [Info] ... MICA-proc_struct_surf -no-output: done!
else
    echo [Warning] ... Pass MICA-proc_struct: files exist
fi
printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# MICAPIPE-post_struct
parc_post_struct="${out_sub_dir}/parc/sub-${sub}_space-nativepro_T1w_atlas-${atlas}.nii.gz"
# parc_post_struct2="${out_sub_dir}/parc/sub-${sub}_space-nativepro_T1w_atlas-${atlas2}.nii.gz"

if [ -f "${parc_post_struct}" ]; then
    echo [Warning] ... Pass MICA-post_struct: files exist
else
    echo [Info] ... MICA-post_struct -no-output
    micapipe_cleanup -sub ${sub} -out ${out} -bids ${bids} -post_structural -fs_licence ${fslc} &>/dev/null ###### Sat Dec 9 14:42:27 CST 2023
    micapipe -sub ${sub} -out ${out} -bids ${bids} -tmpDir ${tmp} -threads ${nthr} -fs_licence ${fslc} \
        -post_structural -atlas ${atlas} \
        &>/dev/null
    echo [Info] ... MICA-post_struct -no-output: done!
fi
printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# MICAPIPE-proc_fun
func_gii_dir="${out_sub_dir}/func/desc-se_dir-AP_rest_bold/surf"

if [ "$(ls -A ${func_gii_dir})" ]; then
    echo [Warning] ... Pass MICA-post_struct: files exist
else
    echo [Info] ... MICA-proc_fun -cleanup
    micapipe_cleanup -sub ${sub} -out ${out} -bids ${bids} -fs_licence ${fslc} -proc_func -acqStr desc-se_dir-AP_rest_bold ###### Sat Dec 9 14:42:27 CST 2023
    wait

    echo [Info] ... MICA-proc_fun -no-output
    micapipe -sub ${sub} -out ${out} -bids ${bids} -tmpDir ${tmp} -threads ${nthr} -fs_licence ${fslc} \
        -proc_func -mainScanStr "dir-AP_rest_bold" -NSR -GSR -dropTR &>/dev/null
    # -func_pe "$bidsbold"

    echo [Info] ... MICA-proc_fun -no-output: done!
fi
#### Wed Sep 18 15:10:00 CST 2024

# ---------------------------------------------------------
# MICAPIPE-proc_fun only change NSR and GSR

python /public/home/baishw/micapipe/micapipe/functions/03_FC.py \
    sub-${sub} \
    ${out}/micapipe_v0.2.0/sub-${sub}/func/desc-se_dir-AP_rest_bold \
    /public/home/baishw/micapipe/micapipe/parcellations/ \
    /public/home/baishw/micapipe/micapipe/parcellations \
    ${out}/micapipe_v0.2.0/sub-${sub}/parc 1 1 _space-func_desc-se FALSE 0

# ---------------------------------------------------------
# replace 5tt file with 5tt freesurfer
mica_5tt_old=${out_sub_dir}/anat/sub-${sub}_space-nativepro_T1w_5tt_mica.nii.gz
t1_5tt=${out_sub_dir}/anat/sub-${sub}_space-nativepro_T1w_5tt.nii.gz
t1_5tt_mif=${out_sub_dir}/anat/5tt-freesurfer.mif
if [ ! -f "${mica_5tt_old}" ]; then
    echo [Info] ... Run 5ttgen-freesurfer
    5ttgen freesurfer ${fs_out_dir}/mri/aseg.mgz ${t1_5tt_mif} -nocrop -force &>/dev/null
    mv ${t1_5tt} ${mica_5tt_old} &>/dev/null
    mrconvert ${t1_5tt_mif} ${t1_5tt} &>/dev/null
    echo [Info] ... Run 5ttgen-freesurfer: done!
else
    echo [Warning] ... Pass replace 5tt file: freesurfer files exist
fi

printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# WMH registration
wmh_to_dwi=${out_sub_dir}/anat/sub-${sub}_space-dwi_wmh.nii.gz
wmh_to_5tt=${out_sub_dir}/anat/sub-${sub}_space-5tt_wmh.nii.gz
if [ ! -f "${wmh_to_5tt}" ]; then
    echo [Info] ... Run wmh registration

    rbidswmh=${out_sub_dir}/anat/sub-${sub}_space-flair_wmh.nii.gz
    # required files
    nativa_t1=${out_sub_dir}/anat/sub-${sub}_space-nativepro_T1w.nii.gz
    fs_t1=${out_sub_dir}/anat/sub-${sub}_space-fsnative_T1w.nii.gz

    # affine DF files
    flair_to_native_affine=${out_sub_dir}/xfm/sub-${sub}_from-flair_to-nativepro_wmh-reg-affine_
    flair_to_fs_affine=${out_sub_dir}/xfm/sub-${sub}_from-flair_to-fsnative_wmh-reg-affine_
    dwi_to_native_affine=${out_sub_dir}/xfm/sub-${sub}_space-dwi_from-dwi_to-nativepro_mode-image_desc-affine_0GenericAffine.mat

    # Check if the bidswmh has been resampled to the flair space
    shape_wmh=$(fslhd $bidswmh | grep "^dim1" | awk '{print $2, $3, $4}')
    shape_flair=$(fslhd $bidsflair | grep "^dim1" | awk '{print $2, $3, $4}')
    if [ "$shape_wmh" != "$shape_flair" ]; then
        echo "Shapes are different, resampling bidswmh to match bidsflair"
        # flirt -in $bidswmh -ref $bidsflair -out ${rbidswmh} -interp nearestneighbour
        applywarp --ref=$bidsflair --in=$bidswmh --out=${rbidswmh} --interp=nn
    else
        echo "Shapes are the same, no resampling needed"
    fi

    # generate DF registration to dwi
    bash antsRegistrationSyN.sh -d 3 -f "$nativa_t1" -m "$bidsflair" -o "$flair_to_native_affine" -t a -n ${nthr} -p d #&> /dev/null
    bash antsRegistrationSyN.sh -d 3 -f "$fs_t1" -m "$bidsflair" -o "$flair_to_fs_affine" -t a -n ${nthr} -p d         #&> /dev/null

    # registration to t1 5tt
    antsApplyTransforms -d 3 -i ${bidswmh} -r ${fs_t1} \
        -t ${flair_to_fs_affine}0GenericAffine.mat \
        -o ${wmh_to_5tt} -n NearestNeighbor -v -u -e 3 #&> /dev/null

    # 5ttedit
    5ttedit ${t1_5tt_mif} ${t1_5tt_mif} -path ${wmh_to_5tt} -force
    mrconvert ${t1_5tt_mif} ${t1_5tt} -force &>/dev/null
    echo [Info] ... Run wmh registration: done!
else
    echo [Warning] ... Pass wmh registration: files exist
fi
printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# MICAPIPE dwi
if [ ! -f "${out_sub_dir}/dwi/sub-${sub}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" ]; then
    echo [Info] ... run MICA-DWI
    micapipe_cleanup -sub ${sub} -out ${out} -bids ${bids} -proc_dwi ###### Sat Dec 9 14:44:39 CST 2023
    micapipe -sub ${sub} -out ${out} -bids ${bids} -tmpDir ${tmp} -threads ${nthr} -fs_licence ${fslc} \
        -proc_dwi
    echo [Info] ... run MICA-DWI: done!
else echo [Warning] ... Pass MICA-DWI: DWI folder exist; fi

# Deal with failure of DWI processing
if [ ! -f "${out_sub_dir}/dwi/sub-${sub}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" ]; then exit; fi
printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# create two types of tcks
dwi_space_prex=${out_sub_dir}/dwi/sub-${sub}_space-dwi

# tcks
dwi_tck=${dwi_space_prex}_desc-iFOD2-${tract_num}_tractography.tck
dwi_tck_wb=${dwi_space_prex}-wb_desc-iFOD2-${tract_num}_tractography.tck
dwi_tck_wb_small=${dwi_space_prex}-wb_smaller_tractography.tck
dwi_tck_wmh=${dwi_space_prex}-wmh_desc-iFOD2-${tract_num}_tractography.tck
dwi_tck_wmh_small=${out_sub_dir}/dwi/sub-${sub}_vis_wmh-smaller_tractography.tck
dwi_tck_int=${dwi_space_prex}-int_desc-iFOD2-${tract_num}_tractography.tck
dwi_tck_int_small=${out_sub_dir}/dwi/sub-${sub}_vis_int-smaller_tractography.tck

# tckweights
dwi_tck_wei=${dwi_space_prex}_desc-iFOD2-${tract_num}_tractography_weights.txt
dwi_tck_wei_wb=${dwi_space_prex}-wb_desc-iFOD2-${tract_num}_tractography_weights.txt
dwi_tck_wei_wmh=${dwi_space_prex}-wmh_desc-iFOD2-${tract_num}_tractography_weights.txt
dwi_tck_wei_int=${dwi_space_prex}-int_desc-iFOD2-${tract_num}_tractography_weights.txt

# ---------------------------------------------------------
# MICAPIPE SC
wb_connectome="${out_sub_dir}/dwi/connectomes-wmh-${tract_num}/sub-${sub}_space-dwi_atlas-${atlas}_desc-iFOD2-${tract_num}-SIFT2_full-connectome.shape.gii"
if [ ! -f $wb_connectome ]; then
    echo [Info] ... run MICA-SC
    # always run cleanup before generating SC
    rm -rf "${out_sub_dir}/dwi/connectomes-int-${tract_num}"
    rm -rf "${out_sub_dir}/dwi/connectomes-wmh-${tract_num}"
    rm -rf "${out_sub_dir}/dwi/connectomes-wb-${tract_num}"
    
    micapipe_cleanup -sub ${sub} -out ${out} -bids ${bids} -SC -tracts ${tract_num} -fs_licence ${fslc} &>/dev/null

    # micapipe SC
    micapipe -sub ${sub} -out ${out} -bids ${bids} -tmpDir ${tmp} -threads ${nthr} -fs_licence ${fslc} \
        -SC -tracts ${tract_num} -keep_tck &>/dev/null

    mv ${out_sub_dir}/QC/sub-${sub}_module-SC-${tract_num}.json ${out_sub_dir}/QC/sub-${sub}_module-SC-wb-${tract_num}.json
    mv ${out_sub_dir}/logs/SC-${tract_num}*.txt ${out_sub_dir}/logs/proc_SC-wb-${tract_num}-${lopuu}.txt
    echo [Info] ... run MICA-SC: done!
else
    echo [Warning] ... Pass MICA-SC: connectomes folder exist
fi
if [ ! "$(ls -A ${out_sub_dir}/dwi/connectomes)" ]; then exit; fi
printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# rename connectome folder for wb
if [ ! -d "${out_sub_dir}/dwi/connectomes-wb-${tract_num}" ] || [ ! -d "$dwi_tck_wei_wb" ]; then
    echo [Info] ... Renaming folder and files of wb tck
    mv ${out_sub_dir}/dwi/connectomes ${out_sub_dir}/dwi/connectomes-wb-${tract_num}
    # rename tck weight
    mv $dwi_tck_wei $dwi_tck_wei_wb
    # rename tck json
    mv ${out_sub_dir}/dwi/sub-${sub}_space-dwi_desc-iFOD2-${tract_num}_tractography.json \
        ${out_sub_dir}/dwi/sub-${sub}_space-dwi-wb_desc-iFOD2-${tract_num}_tractography.json
    # rename tck tdi
    mv ${out_sub_dir}/dwi/sub-${sub}_space-dwi_desc-iFOD2-${tract_num}_tdi.nii.gz \
        ${out_sub_dir}/dwi/sub-${sub}_space-dwi-wb_desc-iFOD2-${tract_num}_tdi.nii.gz
    # rename tck
    mv $dwi_tck $dwi_tck_wb

else
    echo [Warning] ... Pass renaming connectome folder: folder exist
fi
printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# tckedit tck
fslsplit ${out_sub_dir}/dwi/sub-${sub}_space-dwi_desc-5tt.nii.gz ${out_sub_dir}/dwi/sub-${sub}_space-dwi_desc-
wmh_to_dwi=${out_sub_dir}/dwi/sub-${sub}_space-dwi_desc-0004.nii.gz

if [ ! -f "${dwi_tck_wmh}" ]; then
    echo [Info] ... Create WMH connected and intact tractography
    echo [Info] ...... Create WMH connected tractography
    tckedit $dwi_tck_wb $dwi_tck_wmh -include $wmh_to_dwi -tck_weights_in $dwi_tck_wei_wb -tck_weights_out $dwi_tck_wei_wmh -force

    echo [Info] ...... Create intact tractography
    tckedit $dwi_tck_wb $dwi_tck_int -exclude $wmh_to_dwi -tck_weights_in $dwi_tck_wei_wb -tck_weights_out $dwi_tck_wei_int -force
else
    echo [Warning] ... Pass WMH connected tractography: files exist
fi
printf -- '-%.0s' {1..20}
echo

# ---------------------------------------------------------
# Create WMH connected SC

if [ ! -d "${out_sub_dir}/dwi/connectomes-wmh-${tract_num}" ] || [ ! -d "${out_sub_dir}/dwi/connectomes-int-${tract_num}" ]; then
    if [ ! -f $dwi_tck_int ] && [ ! -f $dwi_tck_wmh ]; then
        exit
    else
        echo [Info] ... WMH SC

        echo [Info] ... Create WMH SC
        if [ ! -d "${out_sub_dir}/dwi/connectomes-wmh-${tract_num}" ]; then
            # cleanup the existing SC folder
            # note additional args will cause err, which is not the same as docker
            # must set the -tracts arg, otherwise the default is 40M can will cause err

            micapipe_cleanup -sub ${sub} -out ${out} -bids ${bids} -SC -tracts ${tract_num} -fs_licence ${fslc} &>/dev/null
            mv $dwi_tck_wei_wmh $dwi_tck_wei
            micapipe -sub ${sub} -out ${out} -bids ${bids} -tmpDir ${tmp} -threads ${nthr} -fs_licence ${fslc} \
                -SC -tracts ${tract_num} -tck $dwi_tck_wmh -keep_tck False \
                &>/dev/null
            mv ${out_sub_dir}/dwi/connectomes ${out_sub_dir}/dwi/connectomes-wmh-${tract_num}
            mv $dwi_tck_wei $dwi_tck_wei_wmh
            mv ${out_sub_dir}/QC/sub-${sub}_module-SC-${tract_num}.json ${out_sub_dir}/QC/sub-${sub}_module-SC-wmh-${tract_num}.json
            mv ${out_sub_dir}/logs/SC-${tract_num}*.txt ${out_sub_dir}/logs/proc_SC-wmh-${tract_num}-${lopuu}.txt
        fi

        echo [Info] ... Create intact SC
        if [ ! -d "${out_sub_dir}/dwi/connectomes-int-${tract_num}" ]; then
            micapipe_cleanup -sub ${sub} -out ${out} -bids ${bids} -SC -tracts ${tract_num} -fs_licence ${fslc} &>/dev/null
            mv $dwi_tck_wei_int $dwi_tck_wei
            micapipe -sub ${sub} -out ${out} -bids ${bids} -tmpDir ${tmp} -threads ${nthr} -fs_licence ${fslc} \
                -SC -tracts ${tract_num} -tck $dwi_tck_int -keep_tck False \
                &>/dev/null
            mv ${out_sub_dir}/dwi/connectomes ${out_sub_dir}/dwi/connectomes-int-${tract_num}
            mv $dwi_tck_wei $dwi_tck_wei_int
            mv ${out_sub_dir}/QC/sub-${sub}_module-SC-${tract_num}.json ${out_sub_dir}/QC/sub-${sub}_module-SC-int-${tract_num}.json
            mv ${out_sub_dir}/logs/SC-${tract_num}*.txt ${out_sub_dir}/logs/proc_SC-int-${tract_num}-${lopuu}.txt
        fi
    fi
    printf -- '-%.0s' {1..20}
    echo

else
    echo [Warning] ... Pass sub SC: files exist
fi

# ---------------------------------------------------------
# #  Create smaller tractography
# # Create smaller tractography for wmh tck

# if [ ! -f $dwi_tck_wmh_small ]; then
#     echo [Info] ... Create smaller tractography
#     tckedit -number 100000 $dwi_tck_wmh $dwi_tck_wmh_small &>/dev/null
# fi

# # Create smaller tractography for intact tck
# if [ ! -f $dwi_tck_int_small ]; then
#     tckedit -number 100000 $dwi_tck_int $dwi_tck_int_small &>/dev/null
# fi

# ----------------------------------------------------------------
# Remove tck files
# Some subjects can only get either the wmh or int?
if [ -f ${dwi_tck_wmh} ] || [ -f ${dwi_tck_int} ]; then
    echo [Info] ... Remove TCK files
    rm -f $dwi_tck_wmh
    rm -f $dwi_tck_int
    rm -f $dwi_tck_wb
    rm ${out_sub_dir}/dwi/sub-${sub}_space-dwi_desc-preproc_dwi.mif
    rm ${out_sub_dir}/anat/5tt-freesurfer.mif
    rm -rf ${out_sub_dir}/dwi/eddy
else
    echo [Warning] ... Pass remove: TCK files not exist
fi

# rm ${out_sub_dir}/func/desc-se_dir-AP_rest_bold/surf/sub-${sub}_hemi*
# rm ${out_sub_dir}/func/desc-se_dir-AP_rest_bold/surf/sub-${sub}_surf-fsLR-5k_desc-FC.shape.gii
# rm ${out_sub_dir}/func/desc-se_dir-AP_rest_bold/surf/sub-${sub}_surf-fsLR-32k_desc-timeseries_clean.shape.gii

echo "################################"
echo "${sub}: End at ${formatted_start_time}"
echo "################################"
