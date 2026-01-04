#!/bin/bash
###
 # @Author: Ken32g ken32k@163.com
 # @Date: 2024-04-10 09:12:26
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-09-17 21:41:53
 # @FilePath: /csvd-sfc/svd_file-dcm2nii.sh
 # @Description: Dcm2nii
### 

sub=$1
echo CSVD_BATCH_DCM2NII
printf -- '-%.0s' {1..20}; echo
echo start on $(date)
echo proc sub-$sub
printf -- '-%.0s' {1..20}; echo

# Paths
PROJ_HOME=/public/home/baishw/WMH_MCI
flair_wmh_dir="/public_bme/share/RENJI/RENJI/flair_wmh"
dcm_subject_dir=/public_bme/share/RENJI/RENJI/VCI/sub${sub}*
dcm_temp_dir=${PROJ_HOME}/dcm_temp/sub-${sub}
bids=${PROJ_HOME}/csvd_mica_bids
out_dir=${PROJ_HOME}/csvd_mica_bids/sub-${sub}


# Reg for the MRI sequences
# seq_dirs=("*DTI*" "*MPR*" "*lair*" "*BrainWave*")
seq_dirs=("*BrainWave*")
# Target nifti files
# seq_tars=("/dwi/sub-${sub}_dir-AP_dwi.nii.gz" "/anat/sub-${sub}_T1w.nii.gz" "/anat/sub-${sub}_FLAIR.nii.gz" "/func/sub-${sub}_dir-AP_rest_bold.nii.gz")
seq_tars=("/func/sub-${sub}_dir-AP_rest_bold.nii.gz")

if [ ! -d ${dcm_subject_dir} ]; then
    echo  "[Info] ... dcm folder not exits, exit!"
    exit
fi

# Create dir if not exist
if [ ! -d ${out_dir} ]; then
    echo  "[Info] ... create output folder, create dst folder."
    mkdir -p ${out_dir}
fi

# Copy seq dicom folders to the temp folder, if target file not exists
for i in ${!seq_tars[@]};do
    echo ${seq_tars[$i]}-----------------------
    if [ -f ${out_dir}${seq_tars[$i]} ];then
        echo [Info] ... ${seq_tars[$i]} target exits!
    else
        echo [Warning] ... ${seq_tars[$i]} target not exists!
        
        # Create tmp folder
        if [ ! -d $dcm_temp_dir ];then
            echo "[Info] ... create tmp folder."
            mkdir $dcm_temp_dir
        fi
        seq_name=${seq_dirs[$i]}
        
        # Some subjects have secondary folder
        dcm_sub_seq=${dcm_subject_dir}/${seq_name}
        if [ ! "$(ls -A ${dcm_sub_seq})" ]; then
            dcm_sub_seq=${dcm_subject_dir}/*/${seq_name}
        fi
        # Check if multiple scans exist
        seq_count=$(echo $dcm_sub_seq | grep -o "sub${sub}" | wc -l)
        
        if [ $seq_count -gt 1 ]; then
            echo "[Info] ... multiple folders exit."
            max_count=0
            max_seq=""
            for seq in ${dcm_sub_seq[@]};do
                count=$(ls -1 "${seq}" | wc -l)
                if [ $count -gt $max_count ]; then
                    max_count=$count
                    max_seq=$seq
                fi
            done
            seq_dir=$max_seq
        else
            echo "[Info] ... Set seq_dir to $dcm_sub_seq."
            seq_dir=$dcm_sub_seq
        fi
        seq_dir_r=$(echo $seq_dir)
        
        # Copy the folder
        if [ "$(ls -A "$seq_dir_r")" ]; then
            echo "[Info] ... cp dcm files from ${seq_dir_r} to ${dcm_temp_dir}"
            cp -r "$seq_dir_r" "$dcm_temp_dir"
            
        else
            echo ${seq_dir_r}
            echo "[Warning] ... Source folder ${seq_dir_r} not exist or empty."
        fi
    fi
    echo "[Info] ... Copying file done."
    echo -----------------------
done

echo "[Info] ... Perform dcm2niix for all dcm sequences in the temp folder"
dcm2niix -f %d -z o -o ${out_dir} ${dcm_temp_dir}/*/*.dcm

cd ${out_dir}

mkdir -p anat
mkdir -p dwi
mkdir -p func

# mv *MPRAG*.nii* anat/sub-${sub}_T1w.nii.gz
# mv *MPRAG*.json anat/sub-${sub}_T1w.json

# mv *FLAIR*nii* anat/sub-${sub}_FLAIR.nii.gz
# mv *FLAIR*json anat/sub-${sub}_FLAIR.json

# mv *Flair*nii* anat/sub-${sub}_FLAIR.nii.gz
# mv *Flair*json anat/sub-${sub}_FLAIR.json

# mv *DTI*.nii* dwi/sub-${sub}_dir-AP_dwi.nii.gz
# mv *DTI*.json dwi/sub-${sub}_dir-AP_dwi.json
# mv *.bval dwi/sub-${sub}_dir-AP_dwi.bval
# mv *.bvec dwi/sub-${sub}_dir-AP_dwi.bvec

mv *BrainWave*.nii* func/sub-${sub}_dir-AP_rest_bold.nii.gz
mv *BrainWave*.json func/sub-${sub}_dir-AP_rest_bold.json

#  Remove temp folder
echo "[Info] ... rm temp folder."
# rm -r ${dcm_temp_dir}

# match WMH
if [ ! -f "${out_dir}/anat/sub-${sub}_WMH.nii.gz" ]; then
    flair_sub_dir=$(find ${flair_wmh_dir} -type d -name "sub${sub}*" | head -n 1)
    
    # Check if there is a corresponding file in the flair_wmh directory
    if [ -f "${flair_sub_dir}/flair_wmh_mask_0.nii.gz" ]; then
        # Copy the file to the subfolder in the bids_subject directory
        cp "${flair_sub_dir}/flair_wmh_mask_0.nii.gz" "${out_dir}/anat/sub-${sub}_WMH.nii.gz"
        echo "Succeed: ${sub} done!"
    else
        # If there is no corresponding file, print an error message
        echo "Error: ${sub} not found in ${flair_sub_dir}/"
    fi
else
    echo WMH file exists.
fi


# # The latter script 
# # Match FLAIR from wmh folder
# if [ ! -f "${out_dir}/anat/sub-${sub}_FLAIR.nii.gz" ]; then
#     flair_sub_dir=$(find ${flair_wmh_dir} -type d -name "sub${sub}*" | head -n 1)
#     if [ -f "${flair_sub_dir}/flair_0.nii.gz" ]; then
#         # Copy the file to the subfolder in the bids_subject directory
#         cp "${flair_sub_dir}/flair_0.nii.gz" "${out_dir}/anat/sub-${sub}_FLAIR.nii.gz"
#         echo "Succeed: ${sub} done!"
#     else
#         # If there is no corresponding file, print an error message
#         echo "Error: ${sub} not found in ${flair_sub_dir}/"
#         # tree $flair_sub_dir
#     fi

#     # Check if there is a corresponding file in the flair_wmh directory
#     if [ -f "${flair_sub_dir}/flair_0.nii" ]; then
#         # Copy the file to the subfolder in the bids_subject directory
#         cp "${flair_sub_dir}/flair_0.nii" "${out_dir}/anat/sub-${sub}_FLAIR.nii"
#         gzip "${out_dir}/anat/sub-${sub}_FLAIR.nii.gz"
#         echo "Succeed: ${sub} done!"
#     else
#         # If there is no corresponding file, print an error message
#         echo "Error: ${sub} not found in ${flair_sub_dir}/"
#         # tree $flair_sub_dir
#     fi
# else
#     echo WMH file exists.
# fi

# alterT1source="/public_bme/share/RENJI/RENJI/T1"
# if [ ! -f "${out_dir}/anat/sub-${sub}_T1w.nii.gz" ]; then
#     altt1=$(find ${alterT1source}/sub${sub}* -name "*.dcm")
#     echo $altt1
#     if [ ! -d $dcm_temp_dir ];then
#             echo create tmp folder.
#             mkdir $dcm_temp_dir
#     fi
#     cp $altt1 $dcm_temp_dir
#     dcm2niix -f %d -z o -o ${out_dir} ${dcm_temp_dir}
        
#     cd ${out_dir}
#     mv *MPRAG*.nii* anat/sub-${sub}_T1w.nii.gz
#     mv *MPRAG*.json anat/sub-${sub}_T1w.json
#     rm -rf ${dcm_temp_dir}
# else
#     echo WMH file exists.
# fi

# cd ${out_dir}
# mv *mpr*.nii* anat/sub-${sub}_T1w.nii.gz
# mv *mpr*.json anat/sub-${sub}_T1w.json


# if [ ! -f "${out_dir}/dwi/sub-${sub}_dir-AP_dwi.nii.gz" ]; then
#     cd ${out_dir}
#     mv *DTI*20.nii* dwi/sub-${sub}_dir-AP_dwi.nii.gz
#     mv *DTI*20.json dwi/sub-${sub}_dir-AP_dwi.json
#     mv *20.bval dwi/sub-${sub}_dir-AP_dwi.bval
#     mv *20.bvec dwi/sub-${sub}_dir-AP_dwi.bvec
# fi

# Finish
echo end on $(date)