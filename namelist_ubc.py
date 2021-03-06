'''
The namelist of GEFS post-processing pipeline @ WFRT-UBC
'''

import numpy as np

# ========== Modity & Check when migrate to a new envir ========== #
# ---------------------------------------------------------------- #

# The path of near-real-time GEFS
## download to /scratch
path_gefs_nrt_namelist = '/scratch/NowCastingML/GEFS_downloads/{}/'
path_gefs_member_namelist = '/scratch/NowCastingML/GEFS_downloads/{}_members/'

# Name convensions of GEFS individual members
filename_memberc_namelist = 'gec{:02d}.t00z.pgrb2s.0p25.f{:03d}'
filename_memberp_namelist = 'gep{:02d}.t00z.pgrb2s.0p25.f{:03d}'
## Download 29 perturbed members only.
## Member 30 sometimes missing @ NOAA server, drop
ensemble_mumberp_raw_gefs_namelist = 29


# The path of reference data and geographical data
save_dir = '/oper_data/NowCastingML/BASE_DATA/'


# The path of reforecast and reanalysis
BASE_dir = '/oper_data/NowCastingML/BASE_DATA/'


# The path of output files. "{}" means separate on each day
## !!!! <--- NOT output figures. 
## !!!! The path of output figures is @ namelist_plot.py
output_dir_namelist = '/oper_data/NowCastingML/example_output/{}/'


# The path of log files
status_dir_namelist = '/oper_data/NowCastingML/example_output/{}/'


# Output filename, e.g., gefs_nonDL_20220310.hdf
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
ensemble_number_namelist = 25

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

STN_LOCs = {'BCK': (-123.91388888888889, 48.50333333333333),
            'CMX': (-125.09444444444443, 49.643055555555556),
            'COQ': (-122.77777777777777, 49.355555555555554),
            'CRS': (-116.51666666666667, 49.1),
            'ELK': (-125.76416666666667, 49.87277777777778),
            'PYN': (-122.6375, 55.35),
            'STA': (-122.32638888888889, 49.5575),
            'WAH': (-121.61861111111111, 49.231944444444444),
            'WON': (-121.8, 56.733333333333334),
            'YRV': (-118.18333333333334, 50.96666666666667),
            'YVR': (-123.18333333333334, 49.18333333333333),
            'YXJ': (-120.73333333333333, 56.233333333333334),
            'ING': (-125.11, 56.72),
            'DLU': (-123.1672, 50.86)}
