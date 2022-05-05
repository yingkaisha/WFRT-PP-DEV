#!/bin/bash

# The main GEFS pre-processing routine
# ---------------------------------------------------
#     Check if today's GEFS files exists.
#     Download the ensemble mean.
#     (Planned) subset to given variables and domain.

echo
echo "[`date +"%F %T %Z"`] Starting ${0} ..."

current_time=$(date -u +%Y%m%d)
#current_time='20220406'
# ---------- UBC/NCAR Switch ---------- #
# UBC
#download_status="/home/ibcs/GEFS0p25/${current_time}/download.status"
# NCAR
download_status="/glade/scratch/ksha/DATA/GEFS/${current_time}_members/download.status"
# ------------------------------------- #

run=true

while $run; do

    cd $PWD
    python3 ./GEFS/GEFS_download_members.py $current_time
    # ------------------------------------- #
    
    log_info=$(cat $download_status)
    flag_success=$[10#${log_info:0:2}]
    
    if [ $flag_success == "0" ]; then
        echo "No new files, waiting 30 min ..."
        sleep 1800
    else
        run=false
    fi
done

echo "[`date +"%F %T %Z"`] ${0} is done."
exit 0
