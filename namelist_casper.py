'''
The namelist of AnEn-CNN post-processing pipeline @ casper.ucar.edu
'''

# ========== Python modules ========== #

import numpy as np

# ========== Parameters ========== #

# Number of forecast lead times, starts at 9 hrs, can go up to 168 hrs (7 days)
N_leads_namelist = 54

# Years to identify analog days (the training period of AnEn)
year_anen_namelist = np.arange(2000, 2019)

# Years to identify MDSS dependence templates
year_mdss_namelist = np.arange(2000, 2019)

# The number of ensemble member output 
ensemble_number_namelist = 25

# -------------------------------- #

LEADs_namelist = np.arange(0, N_leads_namelist, dtype=np.int)
FCSTs_namelist = np.arange(9.0, 24*7+3, 3)[:N_leads_namelist]

# ========== Data file locations ========== #

# The path of near-real-time GEFS
path_gefs_nrt_namelist = '/glade/scratch/ksha/DATA/GEFS/{}/'
filename_gefs_namelist = 'geavg.t00z.pgrb2s.0p25'

# The path of reference data and geographical data
save_dir = '/glade/work/ksha/data/Keras/BIAS_publish/'

# The path of reforecast and reanalysis
BASE_dir = '/glade/scratch/ksha/DRIVE/BASE_DATA/'

# ---------- Derived path ---------- #

path_domain_namelist = save_dir+'BC_domain_info.hdf'
path_sl_namelist = save_dir+'SL20_d4_unique.hdf'

path_gefs_apcp_namelist = BASE_dir+'BASE_APCP_year{}_lead{}.zarr'
path_gefs_pwat_namelist = BASE_dir+'BASE_PWAT_year{}_lead{}.zarr'
path_era5_namelist = BASE_dir+'BASE_ERA5_year{}_lead{}.zarr'
path_era5_clim_namelist = BASE_dir+'BASE_ERA5_clim.zarr'
path_CNN_namelist = BASE_dir+'AnEn_UNET3M_RAW_tune.hdf'

# ========== Output ========== #

output_dir_namelist = '/glade/scratch/ksha/DRIVE/{}/'
filename_output_namelist = 'gefs_nonDL_{}.hdf'
filename_CNN_output_namelist = 'gefs_CNN_{}.hdf'
# ----------------------------------------------- #


# # Other old path

# DATA_dir = '/glade/scratch/ksha/DATA/'
# BACKUP_dir = '/glade/scratch/ksha/BACKUP/'
# drive_dir = '/glade/scratch/ksha/DRIVE/'

# # ========== Data ========== #

# # NAEFS
# NAEFS_dir = DATA_dir+'NAEFS/'

# # GEFS reforecast
# REFCST_dir = DATA_dir+'REFCST/'
# REFCST_source_dir = drive_dir+'/REFCST/refcstv12/'
# GFS_APCP_single_dir = REFCST_source_dir+'apcp/'
# GFS_PWAT_single_dir = REFCST_source_dir+'pwat/'

# # ERA5 reanalysis
# ERA_dir = BACKUP_dir + 'ERA5/ERA5_PCT/'

# # PRISM
# PRISM_dir = BACKUP_dir + 'PRISM/'

# # ======= Neural network ======= #

# # Path of model check point


# # Path of batch files
# BATCH_dir = DATA_dir+'REFCST_PCT/'
# BATCH_dir2 = DATA_dir+'REFCST_PCT2/'

# # ========== AnEn-CNN ========== #
# # SL search domain indices (; GEFS and ERA5)
# ## [lat0, lat1, lon0, lon1]
# domain_inds = [120, 280, 120, 340]

# # BC domain indices
# ## [lat0, lat1, lon0, lon1]
# bc_inds = [73, 121, 36, 148]

# # Evaluation results


# # ========== Graphics ========== #

# # figure storage
# fig_dir = '/glade/u/home/ksha/figures/'

# # Matplotlib figure export settings
# fig_keys = {'dpi':250, 
#             'orientation':'portrait', 
#             'papertype':'a4',
#             'bbox_inches':'tight', 
#             'pad_inches':0.1, 
#             'transparent':False}

# # colors
# rgb_array = np.array([[0.85      , 0.85      , 0.85      , 1.        ],
#                       [0.66666667, 1.        , 1.        , 1.        ],
#                       [0.33333333, 0.62745098, 1.        , 1.        ],
#                       [0.11372549, 0.        , 1.        , 1.        ],
#                       [0.37647059, 0.81176471, 0.56862745, 1.        ],
#                       [0.10196078, 0.59607843, 0.31372549, 1.        ],
#                       [0.56862745, 0.81176471, 0.37647059, 1.        ],
#                       [0.85098039, 0.9372549 , 0.54509804, 1.        ],
#                       [1.        , 1.        , 0.4       , 1.        ],
#                       [1.        , 0.8       , 0.4       , 1.        ],
#                       [1.        , 0.53333333, 0.29803922, 1.        ],
#                       [1.        , 0.09803922, 0.09803922, 1.        ],
#                       [0.8       , 0.23921569, 0.23921569, 1.        ],
#                       [0.64705882, 0.19215686, 0.19215686, 1.        ],
#                       [0.55      , 0.        , 0.        , 1.        ]])
    
# blue   = rgb_array[3, :]  # blue
# cyan   = rgb_array[2, :]  # cyan
# lgreen = rgb_array[4, :]  # light green
# green  = rgb_array[5, :]  # dark green
# yellow = rgb_array[8, :]  # yellow
# orange = rgb_array[-6, :] # orange
# red    = rgb_array[-3, :] # red
