import sys
import numpy as np
from datetime import datetime, timedelta

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.io.shapereader import natural_earth

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import transforms
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# # custom tools
# sys.path.insert(0, '/glade/u/home/ksha/PUBLISH/WFRT-PP-DEV/')
# sys.path.insert(0, '/glade/u/home/ksha/PUBLISH/WFRT-PP-DEV/libs/')

# import utils

def string_partial_format(fig, ax, x_start, y_start, ha, va, string_list, color_list, fontsize_list, fontweight_list):
    '''
    String partial formatting (experimental).
    
    handles = string_partial_format(fig, ax, 0., 0.5, 'left', 'bottom',
                                    string_list=['word ', 'word ', 'word'], 
                                    color_list=['r', 'g', 'b'], 
                                    fontsize_list=[12, 24, 48], 
                                    fontweight_list=['normal', 'bold', 'normal'])
    Input
    ----------
        fig: Matplotlib Figure instance. Must contain a `canvas` subclass. e.g., `fig.canvas.get_renderer()`
        ax: Matplotlib Axis instance.
        x_start: horizonal location of the text, scaled in [0, 1] 
        y_start: vertical location of the text, scale in [0, 1]
        ha: horizonal alignment of the text, expected to be either "left" or "right" ("center" may not work correctly).
        va: vertical alignment of the text
        string_list: a list substrings, each element can have a different format.
        color_list: a list of colors that matches `string_list`
        fontsize_list: a list of fontsizes that matches `string_list`
        fontweight_list: a list of fontweights that matches `string_list`
    
    Output
    ----------
        A list of Matplotlib.Text instance.
    
    * If `fig` is saved, then the `dpi` keyword must be fixed (becuase of canvas). 
      For example, if `fig=plt.figure(dpi=100)`, then `fig.savefig(dpi=100)`.
      
    '''
    L = len(string_list)
    Handles = []
    relative_loc = ax.transAxes
    renderer = fig.canvas.get_renderer()
    
    for i in range(L):
        handle_temp = ax.text(x_start, y_start, '{}'.format(string_list[i]), ha=ha, va=va,
                              color=color_list[i], fontsize=fontsize_list[i], 
                              fontweight=fontweight_list[i], transform=relative_loc)
        loc_shift = handle_temp.get_window_extent(renderer=renderer)
        relative_loc = transforms.offset_copy(handle_temp._transform, x=loc_shift.width, units='dots')
        Handles.append(handle_temp)
        
    return Handles


def precip_map(CNN_Px, lon, lat, lead_hrs, accum_hrs, dt_utc_now, 
               cmap_precip, label_precip, linewidth_map, Px_str, accum_str, tag, 
               edge, center_lon, shape_watershed_dir, PROVINCE, geom_US, fig_keys, scale_param='50m', font_text=14,):
    '''
    xxx
    '''
    vmin = label_precip[0]
    vmax = label_precip[-1]
    N_colors = len(label_precip)-1 # color = label-1

    dt_ini_str = datetime.strftime(dt_utc_now, '%Y-%h-%d %HZ')

    for lead, fcst_h in enumerate(lead_hrs):

        # --------------------------------------------- #
        # Calculate & format datetime
        dt_valid = dt_utc_now+timedelta(hours=fcst_h)
        dt_accum_start = dt_valid-timedelta(hours=accum_hrs)

        dt_valid_str = datetime.strftime(dt_valid, '%Y-%h-%d %HZ')
        dt_accum0_str = datetime.strftime(dt_accum_start, '%Y-%h-%d %HZ')
        dt_accum1_str = datetime.strftime(dt_valid, '%h-%d %HZ')

        h_str = ' ({} hrs)'.format(int(fcst_h))

        fig = plt.figure(figsize=(13, 11.5), dpi=fig_keys['dpi'])

        # --------------------------------------------- #
        # Map configuration
        proj_ = ccrs.NorthPolarStereo(central_longitude=center_lon)
        ax = fig.gca(projection=proj_)
        ax.set_extent(edge, ccrs.PlateCarree())

        ax.add_feature(cfeature.LAND.with_scale(scale_param), facecolor='none', zorder=1)
        ax.add_feature(cfeature.COASTLINE.with_scale(scale_param), edgecolor='k', linewidth=linewidth_map, zorder=6)
        ax.add_feature(cfeature.BORDERS.with_scale(scale_param), linestyle='--', linewidth=linewidth_map, zorder=6)
        ax.add_feature(PROVINCE, edgecolor='k', linestyle=':', linewidth=linewidth_map, zorder=5)

        ax.add_geometries(Reader(shape_watershed_dir).geometries(), ccrs.PlateCarree(),
                          facecolor='none', edgecolor='0.25', linewidth=linewidth_map, hatch='//', zorder=3) #3
        ax.add_geometries(geom_US, ccrs.PlateCarree(),
                          facecolor='w', edgecolor='none', linewidth=0, zorder=4)
        ax.spines['geo'].set_linewidth(2.5)
        ax.spines['geo'].set_zorder(9)

        # --------------------------------------------- #

        CS_ = ax.contourf(lon, lat, CNN_Px[lead, ...], cmap=cmap_precip, extend='max', transform=ccrs.PlateCarree(),  
                          levels=label_precip, norm=mcolors.BoundaryNorm(label_precip, ncolors=N_colors), zorder=2)

        CS_ = ax.contourf(lon[24:, :46], lat[24:, :46], CNN_Px[lead, 24:, :46], cmap=cmap_precip, extend='max', 
                          transform=ccrs.PlateCarree(), levels=label_precip, 
                          norm=mcolors.BoundaryNorm(label_precip, ncolors=N_colors), zorder=5)

        ax.contour(lon, lat, CNN_Px[lead, ...], levels=label_precip, colors=('0.5',), linewidths=(linewidth_map,), extend='max', 
                   transform=ccrs.PlateCarree(), norm=mcolors.BoundaryNorm(label_precip, ncolors=N_colors), zorder=2)

        # --------------------------------------------- #
        # Title and text boxes
        handle_title = []
        ax_t1 = fig.add_axes([0.14, 0.26, 0.285, 0.045], facecolor='w')
        [j.set_linewidth(0) for j in ax_t1.spines.values()]
        ax_t1.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                          labelleft=False, labeltop=False, labelright=False, labelbottom=False)

        handle_title += string_partial_format(fig, ax_t1, 0, 1, 'left', 'top', 
                                              ['The ', Px_str, ' percentile ', tag], 
                                              ['k',]*4, [font_text,]*4, ['normal', 'bold', 'normal', 'bold'])
        handle_title += string_partial_format(fig, ax_t1, 0, 0.55, 'left', 'top', 
                                              ['GEFS ', accum_str, ' total precip ensemble'], 
                                              ['k',]*3, [font_text,]*3, ['normal', 'bold', 'normal'])

        ax_t2 = fig.add_axes([0.145, 0.215, 0.285, 0.045], facecolor='w')
        [j.set_linewidth(0) for j in ax_t2.spines.values()]
        ax_t2.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                          labelleft=False, labeltop=False, labelright=False, labelbottom=False)

        handle_title += string_partial_format(fig, ax_t2, 0, 1, 'left', 'top', 
                                              ['Initi time: ', dt_ini_str], 
                                              ['k',]*2, [font_text,]*2, ['normal', 'normal'])
        handle_title += string_partial_format(fig, ax_t2, 0, 0.55, 'left', 'top', 
                                              ['Valid time: ', dt_valid_str, h_str], 
                                              ['k',]*3, [font_text,]*3, ['normal', 'normal', 'normal'])


        ax_t3 = fig.add_axes([0.145, 0.205, 0.355, 0.015], facecolor='w')
        [j.set_linewidth(0) for j in ax_t3.spines.values()]
        ax_t3.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                          labelleft=False, labeltop=False, labelright=False, labelbottom=False)

        handle_title += string_partial_format(fig, ax_t3, 0, 1.0, 'left', 'top', 
                                              ['Accum time: from ', dt_accum0_str, ' to ', dt_accum1_str], 
                                              ['k',]*4, [font_text,]*4, ['normal', 'normal', 'normal', 'normal'])

        for handle in handle_title:
            handle.set_bbox(dict(facecolor='w', edgecolor='none', pad=0.0, zorder=6))

        ax_base = fig.add_axes([0.925, 0.25, 0.075, 0.5], facecolor='none')
        [j.set_linewidth(0) for j in ax_base.spines.values()]
        ax_base.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                            labelleft=False, labeltop=False, labelright=False, labelbottom=False)
        cax = inset_axes(ax_base, height='100%', width='25%', borderpad=0, loc=2)
        CBar = plt.colorbar(CS_, orientation='vertical', extend='max', ticks=label_precip, cax=cax)
        CBar.ax.tick_params(axis='y', labelsize=font_text, direction='in', length=17.5)
        CBar.set_label('[mm]', fontsize=font_text)
        CBar.outline.set_linewidth(2.5)

def precip_cmap(return_label=True, accum_map=True):
    #
    color_over = np.array([135, 135, 140])
                          
    if accum_map:
        # colors
        rgb_array = np.array([[200, 200, 255], [160, 160, 255], [120, 120, 255], [80, 80, 255], [0, 0, 255], [0, 0, 200],
                              [0, 200, 0], [0, 255, 0], [50, 255, 50], [100, 255, 100], [150, 255, 150],
                              [255, 255, 150], [255, 255, 100], [255, 255, 50], [255, 255, 0], 
                              [200, 200, 0], [150, 150, 5], [100, 100, 0],
                              [240, 100, 0], [255, 135, 0], [255, 165, 0], [255, 195, 0], 
                              [255, 150, 150], [255, 100, 100], [255, 50, 50], [255, 0, 0], [200, 0, 0],])
        
        label = np.array([1, 5, 10, 15, 30, 40, 50, 60, 70, 80, 90, 100, 
                          125, 150, 175, 200, 225, 250, 300, 350, 400, 
                          450, 500, 600, 700, 800, 900, 1000])


    else:
        rgb_array = np.array([[200, 200, 255], [160, 160, 255], [120, 120, 255], [80, 80, 255], [0, 0, 255],
                              [0, 200, 0], [0, 255, 0], [50, 255, 50], [100, 255, 100], [150, 255, 150],
                              [255, 255, 150], [255, 255, 100], [255, 255, 50], [255, 255, 0],
                              [200, 200, 0], [150, 150, 0], [100, 100, 0],
                              [240, 100, 0], [255, 135, 0], [255, 165, 0], [255, 195, 0],
                              [255, 150, 150], [255, 100, 100], [255, 50, 50], [255, 0, 0], [200, 0, 0],])
        
        label = np.array([0.25, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
                          12.5, 15, 17.5, 20, 22.5, 25, 30, 35, 
                          40, 45, 50, 60, 70, 80, 90, 100])

    cmap_ = mcolors.ListedColormap(rgb_array/255.0, 'precip_cmap')
    cmap_.set_over(color_over/255.0)
    cmap_.set_under('w')
    if return_label:
        return cmap_, label
    else:
        return cmap_

def get_country_geom(name_list, res='10m'):
    country_shapes = Reader(natural_earth(resolution=res, category='cultural', name='admin_0_countries')).records()
    geoms = []
    for name in name_list:
        for shape_temp in country_shapes:
            if name == shape_temp.attributes['NAME_EN']:
                geoms.append(shape_temp.geometry)
    return geoms

def aspc_cal(edge):
    return (edge[1]-edge[0])/(edge[3]-edge[2])


def precip_map(data_pair, lon, lat, lead_hrs, accum_hrs, dt_utc_now, 
               cmap_precip, label_precip, linewidth_map, Px_str, accum_str,
               edge, center_lon, shape_watershed_dir, PROVINCE, geom_US, fig_keys, png_bch_name, scale_param='50m', font_text=14, ):
    '''
    xxx
    '''
    vmin = label_precip[0]
    vmax = label_precip[-1]
    N_colors = len(label_precip)-1 # color = label-1

    dt_ini_str = datetime.strftime(dt_utc_now, '%Y-%h-%d %HZ')
    
    TAGs = ['raw', 'post-processed']
    # use center_lon to separate swbc and bc 
    if center_lon == -125:
        region = 'swbc'
    else:
        region = 'bc'

    # use Px_str to separate pth
    if Px_str == '10-th':
        pth = png_bch_name['pth'][0]
    elif Px_str == '50-th':
        pth = png_bch_name['pth'][1]
    else:
        pth = png_bch_name['pth'][2]

    # format filename. titlename == filename
    accum_ = png_bch_name['accum_'].format(pth, str(accum_hrs))
    dt_fmt = datetime.strftime(dt_utc_now, png_bch_name['dt_fmt_'])

    TITLE_base = png_bch_name['base_'].format(accum_, region, dt_fmt)+png_bch_name['tail']
    
    for lead, fcst_h in enumerate(lead_hrs):
        
        TITLE_text = TITLE_base.format(int(fcst_h))
        handle_title = []
        # --------------------------------------------- #
        # Calculate & format datetime
        dt_valid = dt_utc_now+timedelta(hours=fcst_h)
        dt_accum_start = dt_valid-timedelta(hours=accum_hrs)

        dt_valid_str = datetime.strftime(dt_valid, '%Y-%h-%d %HZ')
        dt_accum0_str = datetime.strftime(dt_accum_start, '%Y-%h-%d %HZ')
        dt_accum1_str = datetime.strftime(dt_valid, '%h-%d %HZ')

        h_str = ' ({} hrs)'.format(int(fcst_h))

        fig = plt.figure(figsize=(1.25*15, 1.25*8), dpi=fig_keys['dpi'])
        gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
        
        proj_ = ccrs.NorthPolarStereo(central_longitude=center_lon)
        
        ax1 = plt.subplot(gs[0, 0], projection=proj_)
        ax2 = plt.subplot(gs[0, 1], projection=proj_)
        plt.subplots_adjust(0, 0, 1, 1, hspace=0, wspace=0)
        
        AX = [ax1, ax2]
        
        TEXT_poos1 = [[0.010, 0.248, 0.20, 0.0450], [0.510, 0.248, 0.20, 0.0450]]
        TEXT_poos2 = [[0.020, 0.198, 0.20, 0.0450], [0.520, 0.198, 0.20, 0.0450]]
        TEXT_poos3 = [[0.020, 0.178, 0.26, 0.0225], [0.520, 0.178, 0.26, 0.0225]]
        
        for i, ax in enumerate(AX):
            if i == 0:
                ax.set_title(TITLE_text, x=0.7, ha='left', fontsize=font_text+4)
            # --------------------------------------------- #
            # Map configuration
            ax.set_extent(edge, ccrs.PlateCarree())

            ax.add_feature(cfeature.LAND.with_scale(scale_param), facecolor='none', zorder=1)
            ax.add_feature(cfeature.COASTLINE.with_scale(scale_param), edgecolor='k', linewidth=linewidth_map, zorder=6)
            ax.add_feature(cfeature.BORDERS.with_scale(scale_param), linestyle='--', linewidth=linewidth_map, zorder=6)
            ax.add_feature(PROVINCE, edgecolor='k', linestyle=':', linewidth=linewidth_map, zorder=5)

            ax.add_geometries(Reader(shape_watershed_dir).geometries(), ccrs.PlateCarree(),
                              facecolor='none', edgecolor='0.25', linewidth=linewidth_map, hatch='//', zorder=3) #3
            ax.add_geometries(geom_US, ccrs.PlateCarree(),
                              facecolor='w', edgecolor='none', linewidth=0, zorder=4)
            ax.spines['geo'].set_linewidth(2.5)
            ax.spines['geo'].set_zorder(9)
            
            # --------------------------------------------- #
            
            data = data_pair[i]
            tag = TAGs[i]
            
            CS_ = ax.contourf(lon, lat, data[lead, ...], cmap=cmap_precip, extend='max', transform=ccrs.PlateCarree(),  
                              levels=label_precip, norm=mcolors.BoundaryNorm(label_precip, ncolors=N_colors), zorder=2)

            CS_ = ax.contourf(lon[24:, :46], lat[24:, :46], data[lead, 24:, :46], cmap=cmap_precip, extend='max', 
                              transform=ccrs.PlateCarree(), levels=label_precip, 
                              norm=mcolors.BoundaryNorm(label_precip, ncolors=N_colors), zorder=5)

            ax.contour(lon, lat, data[lead, ...], levels=label_precip, colors=('0.5',), linewidths=(linewidth_map,), extend='max', 
                       transform=ccrs.PlateCarree(), norm=mcolors.BoundaryNorm(label_precip, ncolors=N_colors), zorder=2)
            
    
            # --------------------------------------------- #
            # Title and text boxes

            ax_t1 = fig.add_axes(TEXT_poos1[i], facecolor='w')
            [j.set_linewidth(0.0) for j in ax_t1.spines.values()]
            ax_t1.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                              labelleft=False, labeltop=False, labelright=False, labelbottom=False)

            handle_title += string_partial_format(fig, ax_t1, 0, 1, 'left', 'top', 
                                                  ['The ', Px_str, ' percentile ', tag], 
                                                  ['k',]*4, [font_text,]*4, ['normal', 'bold', 'normal', 'bold'])
            handle_title += string_partial_format(fig, ax_t1, 0, 0.5, 'left', 'top', 
                                                  ['GEFS ', accum_str, ' total precip ensemble'], 
                                                  ['k',]*3, [font_text,]*3, ['normal', 'bold', 'normal'])

            ax_t2 = fig.add_axes(TEXT_poos2[i], facecolor='w')
            [j.set_linewidth(0.0) for j in ax_t2.spines.values()]
            ax_t2.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                              labelleft=False, labeltop=False, labelright=False, labelbottom=False)

            handle_title += string_partial_format(fig, ax_t2, 0, 1, 'left', 'top', 
                                                  ['Initi time: ', dt_ini_str], 
                                                  ['k',]*2, [font_text,]*2, ['normal', 'normal'])
            handle_title += string_partial_format(fig, ax_t2, 0, 0.525, 'left', 'top', 
                                                  ['Valid time: ', dt_valid_str, h_str], 
                                                  ['k',]*3, [font_text,]*3, ['normal', 'normal', 'normal'])
        
        
            ax_t3 = fig.add_axes(TEXT_poos3[i], facecolor='w')
            [j.set_linewidth(0.0) for j in ax_t3.spines.values()]
            ax_t3.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                              labelleft=False, labeltop=False, labelright=False, labelbottom=False)
        
            handle_title += string_partial_format(fig, ax_t3, 0, 1.0, 'left', 'top', 
                                                  ['Accum time: from ', dt_accum0_str, ' to ', dt_accum1_str], 
                                                  ['k',]*4, [font_text,]*4, ['normal', 'normal', 'normal', 'normal'])
        
        for handle in handle_title:
            handle.set_bbox(dict(facecolor='w', edgecolor='none', pad=0.0, zorder=6))
        
        ax_base = fig.add_axes([1.01, 0.2, 0.065, 0.6], facecolor='none')
        [j.set_linewidth(0) for j in ax_base.spines.values()]
        ax_base.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                            labelleft=False, labeltop=False, labelright=False, labelbottom=False)
        cax = inset_axes(ax_base, height='100%', width='25%', borderpad=0, loc=2)
        CBar = plt.colorbar(CS_, orientation='vertical', extend='max', ticks=label_precip, cax=cax)
        CBar.ax.tick_params(axis='y', labelsize=font_text, direction='in', length=21)
        CBar.set_label('[mm]', fontsize=font_text)
        CBar.outline.set_linewidth(2.5)
        output_dir = png_bch_name['dir']+TITLE_text
        print(output_dir)
        fig.savefig(output_dir, format='png', **fig_keys)

