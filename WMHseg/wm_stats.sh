###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-08-29 08:08:37
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-09-26 22:16:30
# @FilePath: /csvd-sfc/WMHseg/wmh_stats.sh
# @Description: This script calculates the WMH lesion statistics
###
export PROJ_HOME=$CSVD_DIR
wmh_vol_stat=${PROJ_HOME}/data/sub_wmh_vol_lst.csv
echo pid, wmh_vol >$wmh_vol_stat

for sub_dir in ${PROJ_HOME}/csvd_mica_bids/*; do
    (
        # cd ${sub_dir}/anat
        # sub_id=$(basename $sub_dir)
        # bids_flair=${sub_id}_FLAIR.nii.gz
        # wmh_seg=${sub_id}_WMH_seg.nii.gz
        # wmh_mask=${sub_id}_WMH_mask.nii.gz
        cd ${sub_dir}/wmh
        sub_id=$(basename $sub_dir)
        wmh_mask=space-flair_seg-lst.nii.gz
        # Lesion volume statistics
        # fslmaths $wmh_seg -thr 77 -uthr 77 -bin $wmh_mask &>/dev/null
        fslstats $wmh_mask -V | awk -v X=${sub_id} '{print X",",$2}' >>$wmh_vol_stat

        # Visualization
        # mkdir -p ${PROJ_HOME}/view
        # slicer $bids_flair -a ${PROJ_HOME}/view/lst_${sub_id}_FLAIR.png
        # slicer $wmh_mask -a ${PROJ_HOME}/view/lst_${sub_id}_WMH.png
    ) &
    if (($(wc -w <<<$(jobs -p)) % 12 == 0)); then wait; fi
done
