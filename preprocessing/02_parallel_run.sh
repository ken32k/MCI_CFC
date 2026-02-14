#/bin/bash
# Parallel run mci_cfc-mica_run.sh for all subjects in the MCI-CFC dataset. This script will run 3 parallel processes at a time to process the subjects. The output logs will be saved in the logs directory.

export WMH_MCI=/public/home/baishw/WMH_MCI
bids=${WMH_MCI}/mica_bids

# Parallel
for SUB_DIR in "$bids"/*; do
    (
        # sub_pid=$(basename "$SUB_DIR")
        sub_pid="${SUB_DIR: -3}"
        echo "[Info] $(date +'%Y-%m-%d %H:%M:%S') - Processing $sub_pid"
        bash $WMH_MCI/script/preprocessing/mci_cfc-mica_run.sh $sub_pid 12 2 >$WMH_MCI/logs/${sub_pid}.log &>/dev/null
    ) &
    if (($(wc -w <<<$(jobs -p)) % 3 == 0)); then
        echo $(wc -w <<<$(jobs -p))
        echo full wait
        wait
    fi
done
