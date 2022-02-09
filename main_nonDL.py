
'''
* Post-processed precipitation as an unit of mm per 3 hours 
'''

# sys tools
import sys
import time
import argparse
import os.path

# data tools
import h5py
import zarr
import pygrib
import numpy as np
from datetime import datetime

sys.path.insert(0, '/glade/u/home/ksha/PUBLISH/WFRT-PP-DEV/')
sys.path.insert(0, '/glade/u/home/ksha/PUBLISH/WFRT-PP-DEV/libs/')

import nonDL_lib as nDL
import utils
from namelist_casper import * 

# ========== Datetime informtion ========== #

# UTC time
dt_utc_now = datetime.utcnow()
dt_fmt_string = datetime.strftime(dt_utc_now, '%Y%m%d')

# Day of the year (starts from zero)
dt_day_of_year = dt_utc_now.timetuple().tm_yday - 1

# Month (starts frp, zero)
dt_month_from_zero = dt_utc_now.month-1

# Leap vs non-leap year flag
flag_leap_year = utils.leap_year_checker(dt_utc_now.year)

# ========== Check if files exist ========== #

flag_exist = True

print('Checking file existence')

for i, lead in enumerate(LEADs_namelist):
    
    lead_int_h = int(FCSTs_namelist[lead])
    
    # GEFS file path creation
    GEFS_dir_base = path_gefs_nrt_namelist.format(dt_fmt_string)
    GEFS_dir_full = GEFS_dir_base+filename_gefs_namelist+'.f{0:03d}'.format(lead_int_h)
    
    # check file existence
    flag_temp = os.path.isfile(GEFS_dir_full)
    
    # raise in terms of missing files
    if flag_temp is False:
        print('\tNot found: {}'.format(GEFS_dir_full))
    
        flag_exist = False

# exit the main program if any files missing
if flag_exist is False:
    
    print('The main program is terminated for missing files')
    sys.exit()

# ========== Import domain information ========== #

start_time = time.time()
print('Import domain information ...')
with h5py.File(path_domain_namelist, 'r') as h5io:
    lat_bc = h5io['bc_lat'][...] # lats of the BC domain
    lon_bc = h5io['bc_lon'][...] # lons of the BC domain
    land_mask_bc = h5io['land_mask_bc'][...] # selecting OCEAN grids from the BC domain

ocean_mask_bc = np.logical_not(land_mask_bc) # selecting LAND grids from the BC domain

bc_shape = land_mask_bc.shape
N_grids = np.sum(ocean_mask_bc)

# supplemental locations
SL_xy_dict = {}
with h5py.File(path_sl_namelist , 'r') as h5io:
    for i in range(12):
        temp = h5io['mon_{}_inds'.format(i)][...]
        temp = temp.astype(int)
        SL_xy_dict['{}'.format(i)] = temp

SL_xy = tuple(SL_xy_dict.values())
print('\t ... done. {} sec'.format((time.time() - start_time)))

# ========== AnEn post-processing ========== #

start_time = time.time()
print('AnEn post-processing ...')

AnEn_out = np.empty((N_leads_namelist, N_grids, ensemble_number_namelist))

for i, lead in enumerate(LEADs_namelist):
    lead_int_h = int(FCSTs_namelist[lead])
    
    print("\tLead time = {}".format(lead))
    
    # Import reforecast (Numpy arrays saved as zarr format)
    # ------------------------------------------------------ #
    
    APCP = ()
    PWAT = ()
    
    for year in year_anen_namelist:
        apcp_temp = zarr.load(path_gefs_apcp_namelist.format(year, lead))
        pwat_temp = zarr.load(path_gefs_pwat_namelist.format(year, lead))
        
        APCP += (apcp_temp,)
        PWAT += (pwat_temp,)
    
    # Import reanalysis (Numpy arrays saved as zarr format)
    # ----------------------------------------------------- #
    
    ERA5 = ()
    
    for year in year_anen_namelist:
        era_temp = zarr.load(path_era5_namelist.format(year, lead))
        
        ERA5 += (era_temp,)
    
    # Import today's GEFS 
    # ----------------------------------------------------- #
    
    # GEFS file path creation
    GEFS_dir_base = path_gefs_nrt_namelist.format(dt_fmt_string)
    GEFS_dir_full = GEFS_dir_base+filename_gefs_namelist+'.f{0:03d}'.format(lead_int_h)
    
    with pygrib.open(GEFS_dir_full) as grb_io:
        
        # Import APCP from today's forecast
        grb_reader_apcp = grb_io.select(name='Total Precipitation')[0]
        apcp, _, _ = grb_reader_apcp.data(lat1=48.25, lat2=60.00, lon1=-141.0+360, lon2=-113.25+360)
        apcp = np.flipud(apcp) # GEFS default: kg/m**-2 (or mm) per 3 hours
        
        # Import PWAT from today's forecast
        grb_reader_pwat = grb_io.select(name='Precipitable water')[0]
        pwat, _, _ = grb_reader_pwat.data(lat1=48.25, lat2=60.00, lon1=-141.0+360, lon2=-113.25+360)
        pwat = np.flipud(pwat) # GEFS default: kg/m**-2 (or mm) per 3 hours
    
    apcp_flat = apcp[ocean_mask_bc]
    pwat_flat = apcp[ocean_mask_bc]
    
    # AnEn search
    # ------------------------------------------------------ #
    
    AnEn = nDL.analog_search_SL_single_day(dt_day_of_year, year_anen_namelist, apcp_flat, pwat_flat, 
                                           APCP, PWAT, ERA5, ensemble_number_namelist, SL_xy, flag_leap_year)    
    AnEn_out[i, ...] = AnEn
    
print('\t ... done. {} sec'.format((time.time() - start_time)))

# ========== Import MDSS database ========== #

start_time = time.time()
print('Import MDSS database ...')

ERA5_mdss = ()

window_day = 30
N_days = window_day*2 + 1

# loop over years for MDSS training
for year in year_mdss_namelist:
    
    # separate leap year
    if utils.leap_year_checker(year):
        flag_pick = nDL.search_nearby_days(dt_day_of_year, window=30, leap_year=True)
    else:
        flag_pick = nDL.search_nearby_days(dt_day_of_year, window=30, leap_year=False)
        
    flag_pick = flag_pick == 1
    era_all_lead = np.empty((N_days, N_grids, N_leads_namelist))
    
    # loop over lead times
    for i, lead in enumerate(LEADs_namelist):
        era_temp = zarr.load(BASE_dir+'BC_ERA5_year{}_lead{}.zarr'.format(year, lead))
        era_all_lead[..., i] = era_temp[flag_pick, :]
    
    ERA5_mdss += (era_all_lead,)

ERA5_mdss = np.concatenate(ERA5_mdss)

print('\t ... done. {} sec'.format((time.time() - start_time)))

# ========== MDSS post-processing ========== #

start_time = time.time()
print('MDSS post-processing ...')

AnEn_out = np.transpose(AnEn_out, (2, 0, 1))
ERA5_mdss = np.transpose(ERA5_mdss, (0, 2, 1))

flag_pick = nDL.search_nearby_days(dt_day_of_year, window=30, leap_year=True)
flag_clean, count_trial = nDL.MDSS_main(ERA5_mdss, AnEn_out, factor=5, max_trial=50000)

if count_trial > ensemble_number_namelist:
    print('\t Warning: The MDSS does not coverge fully. Number of sequence: {}, expecting: {}'.format(count_trial, ensemble_number_namelist))
else:
    print('\t The MDSS has converged.')

print('Schaake Shuffle starts ...')

ERA5_pick = ERA5_mdss[flag_clean][:ensemble_number_namelist, ...]                                                                                       
AnEn_MDSS_out = nDL.schaake_shuffle(AnEn_out, ERA5_pick)

print('\t ... MDSS done. {} sec'.format((time.time() - start_time)))

# ========== Save output ========== #

print('Save output ...')

anen_grid = np.empty((ensemble_number_namelist, N_leads_namelist,)+bc_shape)
for en in range(ensemble_number_namelist):
    for l in LEADs_namelist:
        anen_grid[en, l, ocean_mask_bc] = AnEn_MDSS_out[en, l, :]
anen_grid[..., land_mask_bc] = np.nan


name_output = filename_output_namelist.format(dt_fmt_string)
utils.save_hdf5((anen_grid,), ['gefs_apcp',], output_dir_namelist, name_output)



