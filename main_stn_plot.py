
import os
import sys
import time
import argparse

import h5py
import pygrib
import numpy as np
from datetime import datetime

import cartopy.feature as cfeature

base_dir = os.getcwd()
lib_dir = base_dir + '/libs/'

sys.path.insert(0, base_dir)
sys.path.insert(0, lib_dir)

import utils
import plot_lib as plib

from namelist_plot import *


# !!!! <---- change to your namelist
from namelist_ubc import * 
#
# ========== Datetime informtion ========== #

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

    
    
with h5py.File(path_domain_namelist, 'r') as h5io:
    lat_bc = h5io['bc_lat'][...] # lats of the BC domain
    lon_bc = h5io['bc_lon'][...] # lons of the BC domain
    lon_4km = h5io['lon_4km'][...]
    lat_4km = h5io['lat_4km'][...]
    land_mask_bc = h5io['land_mask_bc'][...] # selecting OCEAN grids from the BC domain
    land_mask_bc_4km = h5io['land_mask_bc_4km'][...]
# ocean_mask_bc = np.logical_not(land_mask_bc) # selecting LAND grids from the BC domain
grid_shape = lon_bc.shape


output_dir = output_dir_namelist.format(dt_fmt_string)
name_output = filename_CNN_output_namelist.format(dt_fmt_string)

with h5py.File(output_dir+name_output, 'r') as h5io:
    CNN_output = h5io['gefs_apcp'][...]
    
CNN_output[..., land_mask_bc] = np.nan
    

stn_names = list(STN_LOCs.keys())
stn_lons = []
stn_lats = []

for stn in stn_names:
    stn_lons.append(STN_LOCs[stn][0])
    stn_lats.append(STN_LOCs[stn][1])
    
stn_lons = np.array(stn_lons)
stn_lats = np.array(stn_lats)



indx, indy = utils.grid_search(lon_bc, lat_bc, stn_lons, stn_lats)
L_stn = len(stn_names)

Qs = [0.1, 0.25, 0.5, 0.75, 0.9]
Qs_str = ['P10', 'P25', 'P50', 'P75', 'P90']
COLORS = [plib.xcolor('light steel blue'), plib.xcolor('corn flower blue'), plib.xcolor('slate blue')]

LEADs_3H_ind = np.arange(0, N_leads_namelist, dtype=np.int)
LEADs_3H_hrs = np.arange(9.0, 24*7+3, 3)[:N_leads_namelist]



for n in range(L_stn):
    stn_name = stn_names[n]
    ix = indx[n]
    iy = indy[n]
    STN_CNN = CNN_output[:, :, ix, iy]
    
    DATA = {}
    
    for i, q in enumerate(Qs):
        q_str = Qs_str[i]
        DATA['CNN_3_{}'.format(q_str)] = np.quantile(STN_CNN, q, axis=0)
        DATA['CNN_3_max'] =  np.nanmax(STN_CNN, axis=0)
        DATA['CNN_3_min'] =  np.nanmin(STN_CNN, axis=0)
    
    # -------------------- 24 h -------------------- #
    accum_window = 8 # 8x3h = 1 day
    output_freq = 8 # 2x3h = 6h per output
    skip_start = 0 # skip the first 3hr, start from 12hr instead of 9hr
    
    STN_CNN_24, inds_start, inds_end = utils.accum_slide_window_stn(STN_CNN, accum_window, output_freq, skip_start)
    
    LEADs_24H_hrs = LEADs_3H_hrs[inds_end]
    
    for i, q in enumerate(Qs):
        q_str = Qs_str[i]
        DATA['CNN_24_{}'.format(q_str)] = np.quantile(STN_CNN_24, q, axis=0)
        DATA['CNN_24_max'] =  np.nanmax(STN_CNN_24, axis=0)
        DATA['CNN_24_min'] =  np.nanmin(STN_CNN_24, axis=0)
        
    # -------------------- 72 h -------------------- #
    accum_window = 24 # 24x3h = 3 day
    output_freq = 2 # 2x3h = 6h per output
    skip_start = 1 # skip the first 3hr, start from 12hr instead of 9hr 

    STN_CNN_72, inds_start, inds_end = utils.accum_slide_window_stn(STN_CNN, accum_window, output_freq, skip_start)
    
    LEADs_72H_hrs = LEADs_3H_hrs[inds_end]
    
    for i, q in enumerate(Qs):
        q_str = Qs_str[i]
        DATA['CNN_72_{}'.format(q_str)] = np.quantile(STN_CNN_72, q, axis=0)
        DATA['CNN_72_max'] =  np.nanmax(STN_CNN_72, axis=0)
        DATA['CNN_72_min'] =  np.nanmin(STN_CNN_72, axis=0)
    
    
    STN_CNN_accum = np.cumsum(STN_CNN, axis=1)
    LEADS_accum_d0 = np.copy(LEADs_3H_hrs)
    
    for i, q in enumerate(Qs):
        q_str = Qs_str[i]
        DATA['CNN_d0_{}'.format(q_str)] = np.quantile(STN_CNN_accum, q, axis=0)        
        DATA['CNN_d0_max'] =  np.nanmax(STN_CNN_accum, axis=0)
        DATA['CNN_d0_min'] =  np.nanmin(STN_CNN_accum, axis=0)
        
    
    STN_CNN_accum = np.cumsum(STN_CNN[:, 16:], axis=1)
    LEADS_accum_d2 = np.copy(LEADs_3H_hrs[16:])
    
    for i, q in enumerate(Qs):
        q_str = Qs_str[i]
        DATA['CNN_d2_{}'.format(q_str)] = np.quantile(STN_CNN_accum, q, axis=0)        
        DATA['CNN_d2_max'] =  np.nanmax(STN_CNN_accum, axis=0)
        DATA['CNN_d2_min'] =  np.nanmin(STN_CNN_accum, axis=0)
       
    
    STN_CNN_accum = np.cumsum(STN_CNN[:, 24:], axis=1)
    LEADS_accum_d3 = np.copy(LEADs_3H_hrs[24:])
    
    for i, q in enumerate(Qs):
        q_str = Qs_str[i]
        DATA['CNN_d3_{}'.format(q_str)] = np.quantile(STN_CNN_accum, q, axis=0)        
        DATA['CNN_d3_max'] =  np.nanmax(STN_CNN_accum, axis=0)
        DATA['CNN_d3_min'] =  np.nanmin(STN_CNN_accum, axis=0)
    
    STN_CNN_accum = np.cumsum(STN_CNN[:, 40:], axis=1)
    LEADS_accum_d5 = np.copy(LEADs_3H_hrs[40:])
    
    for i, q in enumerate(Qs):
        q_str = Qs_str[i]
        DATA['CNN_d5_{}'.format(q_str)] = np.quantile(STN_CNN_accum, q, axis=0)        
        DATA['CNN_d5_max'] =  np.nanmax(STN_CNN_accum, axis=0)
        DATA['CNN_d5_min'] =  np.nanmin(STN_CNN_accum, axis=0)    
    
    LEADs = (LEADs_3H_hrs, LEADs_24H_hrs, LEADS_accum_d0, LEADS_accum_d2, LEADS_accum_d3, LEADS_accum_d5, LEADs_72H_hrs,)
    accums = [3, 24, 'd0', 'd2', 'd3', 'd5', 72]
    accum_strs = ['3 hourly', 'daily', 'starting on day 0', 'starting on day 2', 'starting on day 3', 'starting on day 5', 'rolling 3-day'] 
    
    plib.plot_stn(DATA, LEADs, accums, accum_strs, dt_utc_now, stn_name=stn_name, COLORS=COLORS, font_text=font_text, fig_keys=fig_keys, png_stn_name=png_stn_name)
    
    
    
    