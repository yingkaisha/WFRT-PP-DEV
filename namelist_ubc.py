'''
The namelist of AnEn-CNN post-processing pipeline @ WFRT-UBC
'''

import numpy as np

# ========== Modity & Check when migrate to a new envir ========== #
# ---------------------------------------------------------------- #

# !!!! <----- The path of near-real-time GEFS
# ????  '/nfs/kitsault/home/ibcs/GEFS0p25/{}'
path_gefs_nrt_namelist = '/glade/scratch/ksha/DATA/GEFS/{}/'


# !!!! <----- The path of reference data and geographical data
# ????  '/oper_data/NowCastingML/BASE_DATA/'
save_dir = '/glade/work/ksha/data/Keras/BIAS_publish/'


# !!!! <----- The path of reforecast and reanalysis
# ????  '/oper_data/NowCastingML/BASE_DATA/'
BASE_dir = '/glade/scratch/ksha/DRIVE/BASE_DATA/'


# !!!! <----- The path of output files. "{}" means separate on each day
# ????  '/oper_data/NowCastingML/example_output/{}'
output_dir_namelist = '/glade/scratch/ksha/DRIVE/{}/'


# !!!! <----- The path of log files
# ????  '/oper_data/NowCastingML/example_output/{}'
status_dir_namelist = '/glade/scratch/ksha/DRIVE/{}/'


# !!!! <----- Output filename, e.g., gefs_nonDL_20220310.hdf
filename_output_namelist = 'gefs_nonDL_{}.hdf' # AnEn and MDSS only
filename_CNN_output_namelist = 'gefs_CNN_{}.hdf' # AnEn+MDSS+CNN
# ---------------------------------------------------------------- #

# ========== Parameters ========== #

# Number of forecast lead times, starts at 9 hrs, can go up to 168 hrs (7 days)
N_leads_namelist = 54

# Years to identify analog days (the training period of AnEn)
year_anen_namelist = np.arange(2000, 2019)

# Years to identify MDSS dependence templates
year_mdss_namelist = np.arange(2000, 2019)

# The number of ensemble member output 
ensemble_number_namelist = 75

# Num of forecast lead time to lists
LEADs_namelist = np.arange(0, N_leads_namelist, dtype=np.int)
FCSTs_namelist = np.arange(9.0, 24*7+3, 3)[:N_leads_namelist]


# ---------- Touch with caution ---------- #

filename_gefs_namelist = 'geavg.t00z.pgrb2s.0p25'

path_domain_namelist = save_dir+'BC_domain_info.hdf'
path_sl_namelist = save_dir+'SL20_d4_unique.hdf'

path_gefs_apcp_namelist = BASE_dir+'BASE_APCP_year{}_lead{}.zarr'
path_gefs_pwat_namelist = BASE_dir+'BASE_PWAT_year{}_lead{}.zarr'
path_era5_namelist = BASE_dir+'BASE_ERA5_year{}_lead{}.zarr'
path_era5_clim_namelist = BASE_dir+'BASE_ERA5_clim.zarr'
path_CNN_namelist = BASE_dir+'AnEn_UNET3M_RAW_tune.hdf'


