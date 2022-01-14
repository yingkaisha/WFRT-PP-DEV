
import sys
import os.path
from glob import glob
from datetime import datetime, timedelta

# data tools
import h5py
import zarr
import numpy as np

# custom tools
sys.path.insert(0, '/glade/u/home/ksha/PUBLISH/WFRT-PP-DEV/')

from namelist_casper import * 

# --------------------------------------------------------------
# Import indices

# supplemental locations
SL_xy_dict = {}

with h5py.File(save_dir+'SL20_d4_unique.hdf', 'r') as h5io:
    IxIy_unique = h5io['unique_inds'][...]
    for i in range(12):
        temp = h5io['mon_{}_inds'.format(i)][...]
        temp = temp.astype(int)
        SL_xy_dict['{}'.format(i)] = temp
        
IxIy_unique = IxIy_unique.astype(int)
SL_xy = tuple(SL_xy_dict.values())

# --------------------------------------------------------------

for year in range(2000, 2020):
    for lead in range(N_fcst):
        
        print("Save hdf to zarr: year{}, lead{}".format(year, lead))

        with h5py.File(ERA_dir+'ERA5_GEFS-fcst_{}.hdf'.format(year), 'r') as h5io:
            era_temp = h5io['era_fcst'][:, lead, ...][..., IxIy_unique[:, 0], IxIy_unique[:, 1]]

        name_output_era = 'BASE_ERA5_year{}_lead{}.zarr'.format(year, lead)

        zarr.save(BASE_dir+name_output_era, era_temp)


        