
import os
import sys
import time
import argparse

import h5py
import zarr
import pygrib
import numpy as np
from datetime import datetime

from tensorflow import keras
from tensorflow.keras.layers import Input
from tensorflow.keras.models import Model

base_dir = os.getcwd()
lib_dir = base_dir + '/libs/'

sys.path.insert(0, base_dir)
sys.path.insert(0, lib_dir)

import utils
import DL_lib as DL
import nonDL_lib as nDL

# !!!! <---- change to your namelist
from namelist_casper import * 
#

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('date_str', help='date_str. e.g., 20100101')
args = vars(parser.parse_args())

# Datetime info
dt_fmt = args['date_str']

if dt_fmt == 'auto':
    dt_utc_now = datetime.utcnow()
    dt_fmt_string = datetime.strftime(dt_utc_now, '%Y%m%d')
else:
    dt_fmt_string = dt_fmt
    dt_utc_now = datetime.strptime(dt_fmt_string, '%Y%m%d')

# ========== Check nonDL outputs ========== #

print('Checking file existence')

name_output = filename_output_namelist.format(dt_fmt_string)
output_dir = output_dir_namelist.format(dt_fmt_string)

flag_exist = os.path.isfile(output_dir+name_output)
    
if flag_exist:
    with h5py.File(output_dir+name_output, 'r') as h5io:
        gefs_apcp = h5io['gefs_apcp'][...]
else:
    sys.exit('The main program is terminated because of missing files')
    
# ========== Import domain information ========== #

start_time = time.time()

print('Import domain information ...')
with h5py.File(path_domain_namelist, 'r') as h5io:
    lat_bc = h5io['bc_lat'][...] # lats of the BC domain
    lon_bc = h5io['bc_lon'][...] # lons of the BC domain
    etopo_025 = h5io['etopo_bc'][...]
    land_mask_bc = h5io['land_mask_bc'][...] # selecting OCEAN grids from the BC domain
grid_shape = land_mask_bc.shape

# normalize elevation
etopo_025[etopo_025<0] = 0
max_ = np.nanmax(etopo_025)
min_ = np.nanmin(etopo_025)
etopo_025 = (etopo_025-min_)/(max_-min_)

# ERA5 climatology (one of the CNN inputs)

era_3h_clim = zarr.load(path_era5_clim_namelist)

# select current month
mon_ind = dt_utc_now.month-1
temp_clim = era_3h_clim[mon_ind, ...]

# normalize clim
temp_clim = np.log(temp_clim+1)
temp_clim[..., land_mask_bc] = 0.0

print('\t ... done. {} sec'.format((time.time() - start_time)))

# ========== CNN preparation ========== #

# Normalize precip
gefs_dress = gefs_apcp
gefs_dress[..., land_mask_bc] = 0.0 # nan --> 0.0

temp_precip = gefs_dress
temp_precip[temp_precip<0] = 0
temp_precip = np.log(temp_precip+1)

# CNN hyperparam
input_tensor = Input((None, None, 3))
filter_num_down = [80, 160, 320, 640]
filter_num_skip = [80, 80, 80,]
filter_num_aggregate = 320
filter_num_sup = 80
stack_num_down = 2
stack_num_up = 1
activation = 'GELU'
batch_norm = True
pool = False
unpool = False
name = 'denoise'

# CNN config
X_decoder = DL.denoise_base(input_tensor, filter_num_down, filter_num_skip, filter_num_aggregate, 
                            stack_num_down=stack_num_down, stack_num_up=stack_num_up, activation=activation, 
                            batch_norm=batch_norm, pool=pool, unpool=unpool, name=name)

OUT_stack = DL.denoise_sup_head(X_decoder, filter_num_sup, activation=activation, 
                                batch_norm=batch_norm, pool=pool, unpool=unpool, name=name)

model = Model([input_tensor,], OUT_stack)

# CNN weights

W = DL.dummy_loader(path_CNN_namelist)
model.set_weights(W)

# Allocate CNN outputs

CNN_output = np.empty(temp_precip.shape)
N_seq = temp_precip.shape[0]
N_leads = temp_precip.shape[1]
single_seq = np.empty((N_seq,)+grid_shape+(3,)) # Three predictors: GEFS, clim, elev

# ========== CNN prediction ========== #

start_time = time.time()
print('CNN prediction starts')

for lead in range(N_leads):
    print('\tLead time index = {}'.format(lead))
    single_seq[..., 0] = temp_precip[:, lead, ...]
    single_seq[..., 1] = temp_clim[None, ...]
    single_seq[..., 2] = etopo_025[None, ...]
    
    single_out = model.predict([single_seq,])
    single_out = single_out[-1][..., 0]
    single_out = np.exp(single_out)-1
    
    CNN_output[:, lead, ...] = single_out
    
print('\t ... done. {} sec'.format((time.time() - start_time)))

# ========== save ========== #

CNN_output = nDL.cnn_precip_fix(CNN_output)
CNN_output[..., land_mask_bc] = np.nan

name_output = filename_CNN_output_namelist.format(dt_fmt_string)
utils.save_hdf5((CNN_output,), ['gefs_apcp',], output_dir, name_output)

# if the script runs to here, write status=1
status=1
with open(output_dir+"nowcast_kyle.status", "w") as log_io:
    log_io.write(str(status))


