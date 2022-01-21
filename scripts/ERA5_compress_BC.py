
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
with h5py.File(save_dir+'BC_domain_info.hdf', 'r') as h5io:
    lat_bc = h5io['bc_lat'][...] # lats of the BC domain
    lon_bc = h5io['bc_lon'][...] # lons of the BC domain
    land_mask_bc = h5io['land_mask_bc'][...] # selecting OCEAN grids from the BC domain

ocean_mask_bc = np.logical_not(land_mask_bc) # selecting LAND grids from the BC domain

# --------------------------------------------------------------

for year in range(2000, 2020):
    for lead in range(N_fcst):
        
        print("Save hdf to zarr: year{}, lead{}".format(year, lead))

        with h5py.File(ERA_dir+'ERA5_GEFS-fcst_{}.hdf'.format(year), 'r') as h5io:
            era_temp = h5io['era_fcst'][:, lead, ...][:, bc_inds[0]:bc_inds[1], bc_inds[2]:bc_inds[3]][:, ocean_mask_bc]

        name_output_era = 'BC_ERA5_year{}_lead{}.zarr'.format(year, lead)

        zarr.save(BASE_dir+name_output_era, era_temp)


        