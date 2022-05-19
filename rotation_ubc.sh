#!/bin/bash

# The main GEFS pre-processing routine
# ---------------------------------------------------
#     Check if today's GEFS files exists.
#     Download the ensemble mean.
#     (Planned) subset to given variables and domain.

module load oper/NowCastingML

echo
echo "[`date +"%F %T %Z"`] Starting ${0} ..."

current_time=$(date -u +%Y%m%d)

# <---------- !!!! status file
# Double check --> "status_dir_namelist" @ namelist_ubc.py
nowcast_status="/oper_data/NowCastingML/example_output/${current_time}/nowcast_kyle.status"
download_status="/scratch/NowCastingML/GEFS_downloads/${current_time}_members/download.status"
# ------------------------------------- #

run=true

while $run; do
    
    cd $PWD
    bash ./GEFS/main.sh
    bash ./GEFS/main_members.sh
    
    download_info=$(cat $download_status)
    flag_download=$[10#${download_info:0:2}]
    
    #echo $flag_download_member
    if [ $flag_download == "0" ]; then
        echo "GEFS member download failed, re-try in 30 min ..."
        sleep 1800
    else
        #  main nonDL
        echo "Analog ensemble routine starts"
        python3 main_nonDL.py $current_time
        
        log_info=$(cat $nowcast_status)
        flag_success=$[10#${log_info:0:2}]
        
        if [ $flag_success == "0" ]; then
            echo "Raw forecasts do not exist or incomplete, re-try in 30 min ..."
            sleep 1800
        else
            echo "Deep learning routine starts"
            python3 main_DL.py $current_time
            echo "Plot routine starts"
            python3 main_plot.py $current_time
         
            rm -r /scratch/NowCastingML/GEFS_downloads/${current_time}/
            rm -r /scratch/NowCastingML/GEFS_downloads/${current_time}_members/

            echo "Delete raw GEFS to free space"
            run=false
        fi
    fi
done

echo "[`date +"%F %T %Z"`] ${0} is done."
exit 0
