
# sys tools
import sys
import time
import argparse

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

# parser = argparse.ArgumentParser()
# parser.add_argument('year_fcst', help='year_fcst')
# parser.add_argument('part', help='part')
# args = vars(parser.parse_args())
# year_fcst = int(args['year_fcst'])
# part_ = int(args['part'])



dt_utc_now = datetime.utcnow()
dt_fmt_string = datetime.strftime(dt_utc_now, '%Y%m%d')
dt_day_of_year = dt_utc_now.timetuple().tm_yday
dt_month_from_zero = dt_utc_now.month-1
flag_leap_year = utils.leap_year_checker(dt_utc_now.year)

# ========== Parse namelist args ========== #

EN = ensemble_number_namelist

LEADs = LEADs_namelist
N_leads = N_leads_namelist

year_analog = year_anen_namelist
year_mdss = year_mdss_namelist


# importing domain information
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

AnEn_out = np.empty((N_leads, N_grids, EN))

for i, lead in enumerate(LEADs):
    lead_int_h = int(FCSTs_namelist[lead])
    print("Processing lead time = {}".format(lead))
    #print("Main program starts ...")
    # ------------------------------------------------- #
    
    # Import reforecast
    APCP = ()
    PWAT = ()
    
    for year in year_analog:
        apcp_temp = zarr.load(path_gefs_apcp_namelist.format(year, lead))
        pwat_temp = zarr.load(path_gefs_pwat_namelist.format(year, lead))
        
        APCP += (apcp_temp,)
        PWAT += (pwat_temp,)
        
    # ------------------------------------------------- #
    
    # Import reanalysis
    ERA5 = ()
    
    for year in year_analog:
        era_temp = zarr.load(path_era5_namelist.format(year, lead))
        
        ERA5 += (era_temp,)
        
    # ------------------------------------------------- #
    
    GEFS_dir_base = path_gefs_nrt_namelist.format(dt_fmt_string)
    GEFS_dir_full = GEFS_dir_base+filename_gefs_namelist+'.f{0:03d}'.format(lead_int_h)
    
    with pygrib.open(GEFS_dir_full) as grb_io:
        #
        grb_reader_apcp = grb_io.select(name='Total Precipitation')[0]
        apcp, _, _ = grb_reader_apcp.data(lat1=48.25, lat2=60.00, lon1=-141.0+360, lon2=-113.25+360)
        apcp = np.flipud(apcp) # GEFS default: kg/m**-2 (or mm) per 3 hours
        
        #
        grb_reader_pwat = grb_io.select(name='Precipitable water')[0]
        pwat, _, _ = grb_reader_pwat.data(lat1=48.25, lat2=60.00, lon1=-141.0+360, lon2=-113.25+360)
        pwat = np.flipud(pwat) # GEFS default: kg/m**-2 (or mm) per 3 hours
    
    apcp_flat = apcp[ocean_mask_bc]
    pwat_flat = apcp[ocean_mask_bc]
    
    start_time = time.time()
    AnEn = nDL.analog_search_SL_single_day(dt_day_of_year, year_analog, apcp_flat, pwat_flat, APCP, PWAT, ERA5, EN, SL_xy, flag_leap_year)
    #print("... Completed. Time = {} sec ".format((time.time() - start_time)))
    
    AnEn_out[i, ...] = AnEn
    
    



