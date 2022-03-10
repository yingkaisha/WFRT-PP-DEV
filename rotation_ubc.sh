#!/bin/bash

# The main GEFS pre-processing routine
# ---------------------------------------------------
#     Check if today's GEFS files exists.
#     Download the ensemble mean.
#     (Planned) subset to given variables and domain.

echo
echo "[`date +"%F %T %Z"`] Starting ${0} ..."

current_time=$(date -u +%Y%m%d)

# <---------- !!!! status file
# Double check --> "status_dir_namelist" @ namelist_ubc.py
nowcast_status="/glade/scratch/ksha/DRIVE/${current_time}/nowcast_kyle.status"
# ------------------------------------- #

run=true

while $run; do

    # main nonDL
    python3 main_nonDL.py $current_time
    
    log_info=$(cat $nowcast_status)
    flag_success=$[10#${log_info:0:2}]
    
    if [ $flag_success == "0" ]; then
        echo "Raw forecasts do not exist, waiting 30 min ..."
        sleep 1800
    else
        run=false
    fi
done

echo "[`date +"%F %T %Z"`] ${0} is done."
exit 0


