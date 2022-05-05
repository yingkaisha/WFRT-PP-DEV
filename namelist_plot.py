
fig_keys = {'dpi':250, 
            'orientation':'portrait', 
            'papertype':'a4',
            'bbox_inches':'tight', 
            'pad_inches':0.1, 
            'transparent':False}

scale_param = '50m'

font_text = 14
linewidth_map = 1.5

center_lon_bc = -120
edge_bc = [-141, -113.25, 48.25, 60]

center_lon_sw = -125
edge_sw = [-128.85, -121.15, 48, 51.75]

shape_watershed_dir = '/glade/u/home/ksha/figures/wshed_hires/MajorHydroWatershedsProject.shp'

png_bch_name = {}
png_bch_name['base_'] = 'gridbeta.pcp.{}.{}.{}F'
png_bch_name['accum_'] = 'qpf{}hb'
png_bch_name['dt_fmt_'] = '%Y%m%d%H'
png_bch_name['tail'] = '{:03d}.png'
png_bch_name['dir'] = '/glade/scratch/ksha/DRIVE/'
#png_bch_name['dir'] = '/www/results/gridbeta/YYMMDDHH/PNG/g2/'