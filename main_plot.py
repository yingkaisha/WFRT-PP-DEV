
'''
* Post-processed precipitation as an unit of mm per 3 hours 
'''

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

# -------------------- Geo data -------------------- #

with h5py.File(path_domain_namelist, 'r') as h5io:
    lat_bc = h5io['bc_lat'][...] # lats of the BC domain
    lon_bc = h5io['bc_lon'][...] # lons of the BC domain
    lon_4km = h5io['lon_4km'][...]
    lat_4km = h5io['lat_4km'][...]
    land_mask_bc = h5io['land_mask_bc'][...] # selecting OCEAN grids from the BC domain
    land_mask_bc_4km = h5io['land_mask_bc_4km'][...]
# ocean_mask_bc = np.logical_not(land_mask_bc) # selecting LAND grids from the BC domain
grid_shape = lon_bc.shape

# US states and CAN-US boundary
PROVINCE = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale=scale_param,
    facecolor='none')

name_list = ['United States of America']
geom_US = plib.get_country_geom(name_list, scale_param)

# -------------------- GEFS outputs -------------------- #

output_dir = output_dir_namelist.format(dt_fmt_string)
name_output = filename_CNN_output_namelist.format(dt_fmt_string)

with h5py.File(output_dir+name_output, 'r') as h5io:
    CNN_output = h5io['gefs_apcp'][...]
    
CNN_output[..., land_mask_bc] = np.nan

# -------------------- Raw GEFS members -------------------- #

GEFS_raw = np.empty((ensemble_mumber, p_raw_gefs_namelist+1, N_leads_namelist)+grid_shape)

# GEFS file path creation
GEFS_dir_base = path_gefs_member_namelist.format(dt_fmt_string)

for i, lead_ in enumerate(FCSTs_namelist):

    print("Process raw GEFS members on lead_time[{}]".format(i))
    
    for member_ in range(ensemble_mumberp_raw_gefs_namelist+1):
        
        if member_ == 0:
            filename_ = filename_memberc_namelist.format(member_, int(lead_))
        else:
            filename_ = filename_memberp_namelist.format(member_, int(lead_))
    
        GEFS_dir_full = GEFS_dir_base+filename_

        with pygrib.open(GEFS_dir_full) as grb_io:
            grb_reader_apcp = grb_io.select(name='Total Precipitation')[0]
            apcp, _, _ = grb_reader_apcp.data(lat1=48.25, lat2=60.00, lon1=-141.0+360, lon2=-113.25+360)
        
        apcp = np.flipud(apcp)
        GEFS_raw[member_, i, ...] = apcp
        
GEFS_raw[..., land_mask_bc] = np.nan
DATA = {}

# ==================== 3 hourly results ==================== #
        
LEADs_3H_ind = np.arange(0, N_leads_namelist, dtype=np.int)
LEADs_3H_hrs = np.arange(9.0, 24*7+3, 3)[:N_leads_namelist]

DATA['CNN_3_P10'] = np.quantile(CNN_output, 0.1, axis=0)
DATA['CNN_3_P50'] = np.quantile(CNN_output, 0.5, axis=0)
DATA['CNN_3_P90'] = np.quantile(CNN_output, 0.9, axis=0)

DATA['GEFS_3_P10'] = np.quantile(GEFS_raw, 0.1, axis=0)
DATA['GEFS_3_P50'] = np.quantile(GEFS_raw, 0.5, axis=0)
DATA['GEFS_3_P90'] = np.quantile(GEFS_raw, 0.9, axis=0)

# ==================== daily results ==================== #

accum_window = 8 # 8x3h = 1 day
output_freq = 2 # 2x3h = 6h per output
skip_start = 1 # skip the first 3hr, start from 12hr instead of 9hr 

CNN_accum_24, inds_start, inds_end = utils.accum_slide_window(CNN_output, accum_window, output_freq, skip_start)
CNN_accum_24[..., land_mask_bc] = np.nan

LEADs_24H_hrs = LEADs_3H_hrs[inds_end]

GEFS_accum_24, inds_start, inds_end = utils.accum_slide_window(GEFS_raw, accum_window, output_freq, skip_start)
GEFS_accum_24[..., land_mask_bc] = np.nan

DATA['CNN_24_P10'] = np.quantile(CNN_accum_24, 0.1, axis=0)
DATA['CNN_24_P50'] = np.quantile(CNN_accum_24, 0.5, axis=0)
DATA['CNN_24_P90'] = np.quantile(CNN_accum_24, 0.9, axis=0)

DATA['GEFS_24_P10'] = np.quantile(GEFS_accum_24, 0.1, axis=0)
DATA['GEFS_24_P50'] = np.quantile(GEFS_accum_24, 0.5, axis=0)
DATA['GEFS_24_P90'] = np.quantile(GEFS_accum_24, 0.9, axis=0)

# ==================== 3 days results ==================== #

accum_window = 24 # 24x3h = 3 day
output_freq = 2 # 2x3h = 6h per output
skip_start = 1 # skip the first 3hr, start from 12hr instead of 9hr 

CNN_accum_72, inds_start, inds_end = utils.accum_slide_window(CNN_output, accum_window, output_freq, skip_start)
CNN_accum_72[..., land_mask_bc] = np.nan

LEADs_72H_hrs = LEADs_3H_hrs[inds_end]

GEFS_accum_72, inds_start, inds_end = utils.accum_slide_window(GEFS_raw, accum_window, output_freq, skip_start)
GEFS_accum_72[..., land_mask_bc] = np.nan

DATA['CNN_72_P10'] = np.quantile(CNN_accum_72, 0.1, axis=0)
DATA['CNN_72_P50'] = np.quantile(CNN_accum_72, 0.5, axis=0)
DATA['CNN_72_P90'] = np.quantile(CNN_accum_72, 0.9, axis=0)

DATA['GEFS_72_P10'] = np.quantile(GEFS_accum_72, 0.1, axis=0)
DATA['GEFS_72_P50'] = np.quantile(GEFS_accum_72, 0.5, axis=0)
DATA['GEFS_72_P90'] = np.quantile(GEFS_accum_72, 0.9, axis=0)

# ==================== Plot ==================== #

camp_precip, label_precip = plib.precip_cmap(return_label=True, accum_map=False)
lon = lon_bc
lat = lat_bc

ACCUMs = ['3', '24', '72',]
ACCUM_strs = ['3 hourly', '1 day', '3 days'] 
ACCUM_flags = [False, True, True,]
LEADs = (LEADs_3H_hrs, LEADs_24H_hrs, LEADs_72H_hrs,)

Ps = ['P10', 'P50', 'P90']
P_strs = ['10-th', '50-th', '90-th']

EDGEs = [edge_bc, edge_sw]
CENTERs = [center_lon_bc, center_lon_sw]

    
for j, accum in enumerate(ACCUMs):
    accum_hrs = int(accum)
    accum_str = ACCUM_strs[j]
    lead_hrs = LEADs[j]
    cmap_precip, label_precip = plib.precip_cmap(return_label=True, accum_map=ACCUM_flags[j])

    for k, P in enumerate(Ps):
        P_str = P_strs[k]

        data_pair = (DATA['GEFS_{}_{}'.format(accum, P)], DATA['CNN_{}_{}'.format(accum, P)])
        
        for l, edge in enumerate(EDGEs):
            center_lon = CENTERs[l]
            
            plib.precip_map(data_pair, lon, lat, lead_hrs, accum_hrs, dt_utc_now, 
                            camp_precip, label_precip, linewidth_map, P_str, accum_str, 
                            edge, center_lon, shape_watershed_dir, PROVINCE, geom_US, 
                            fig_keys, png_bch_name, scale_param='50m', font_text=14)
            
        
