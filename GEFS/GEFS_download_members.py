from __future__ import absolute_import
from namelist import *

import os
import argparse
import subprocess
import urllib.request
from datetime import datetime

print("[{}] Starting GEFS_Download.py ...".format(datetime.strftime(datetime.now(), '%F %T %z')))

parser = argparse.ArgumentParser()
parser.add_argument('date_str', help='date_str. e.g., 20100101')
args = vars(parser.parse_args())

dt_fmt = args['date_str']
if dt_fmt == 'auto':
    # Get the current UTC datetime
    dt_fmt = datetime.strftime(datetime.utcnow(), '%Y%m%d')

print("Downloading GEFS for {}".format(dt_fmt))

# set up download directory
target_dir = target_path_members.format(dt_fmt)

if os.path.isdir(target_dir):
    print("{} already exists".format(target_dir))
else:
    print("Creating {}".format(target_dir))
    os.mkdir(target_dir)

# downloading status
status = 0

for lead_ in fcst_leads:
    
    for member_ in range(N_member):
        
        if member_ == 0:
            filename_ = filename_memberc.format(member_, lead_)
        else:
            filename_ = filename_memberp.format(member_, lead_)
        
        print("Downloading {}".format(filename_))

        try:
            urllib.request.urlretrieve(url_fmt.format(dt_fmt, filename_), os.path.join(target_dir, filename_));
            status = 1
        except:
            # if file access is not successful, terminate all downloads
            print("{} is not available".format(filename_))
            status = 0
            break;
        
with open("{}/download.status".format(target_dir), "w") as log_io:
    log_io.write(str(status))

print("[{}] GEFS_Download.py is done.".format(datetime.strftime(datetime.now(), '%F %T %Z')))
