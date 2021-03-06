import sys
import numpy as np
from datetime import datetime, timedelta

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.io.shapereader import natural_earth

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
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

        output_dir_now =  datetime.strftime(dt_utc_now, output_dir)
        print(output_dir_now)
        fig.savefig(output_dir_now, format='png', **fig_keys)


def ax_decorate_box(ax, font_text, left_flag=True, bottom_flag=True):
    ax.xaxis.set_tick_params(labelsize=font_text)
    ax.yaxis.set_tick_params(labelsize=font_text)
    [j.set_linewidth(2.5) for j in ax.spines.values()]
    ax.tick_params(axis="both", which="both", bottom=False, top=False, labelbottom=bottom_flag, left=False, right=False, labelleft=left_flag)
    return ax

def xcolor(key):
    xcolor = {
    "maroon":"#800000", "dark red":"#8B0000", "brown":"#A52A2A", "firebrick":"#B22222", "crimson":"#DC143C", "red":"#FF0000",
    "tomato":"#FF6347", "coral":"#FF7F50", "indian red":"#CD5C5C", "light coral":"#F08080", "dark salmon":"#E9967A", "salmon":"#FA8072",
    "light salmon":"#FFA07A", "orange red":"#FF4500", "dark orange":"#FF8C00", "orange":"#FFA500", "gold":"#FFD700", "dark golden rod":"#B8860B",
    "golden rod":"#DAA520", "pale golden rod":"#EEE8AA", "dark khaki":"#BDB76B", "khaki":"#F0E68C", "olive":"#808000", "yellow":"#FFFF00",
    "yellow green":"#9ACD32", "dark olive green":"#556B2F", "olive drab":"#6B8E23", "lawn green":"#7CFC00", "chart reuse":"#7FFF00", "green yellow":"#ADFF2F",
    "dark green":"#006400", "green":"#008000", "forest green":"#228B22", "lime":"#00FF00", "lime green":"#32CD32", "light green":"#90EE90",
    "pale green":"#98FB98", "dark sea green":"#8FBC8F", "medium spring green":"#00FA9A", "spring green":"#00FF7F", "sea green":"#2E8B57", "medium aqua marine":"#66CDAA",
    "medium sea green":"#3CB371", "light sea green":"#20B2AA", "dark slate gray":"#2F4F4F", "teal":"#008080", "dark cyan":"#008B8B", "aqua":"#00FFFF",
    "cyan":"#00FFFF", "light cyan":"#E0FFFF", "dark turquoise":"#00CED1", "turquoise":"#40E0D0", "medium turquoise":"#48D1CC", "pale turquoise":"#AFEEEE",
    "aqua marine":"#7FFFD4", "powder blue":"#B0E0E6", "cadet blue":"#5F9EA0", "steel blue":"#4682B4", "corn flower blue":"#6495ED", "deep sky blue":"#00BFFF",
    "dodger blue":"#1E90FF", "light blue":"#ADD8E6", "sky blue":"#87CEEB", "light sky blue":"#87CEFA", "midnight blue":"#191970",
    "navy":"#000080", "dark blue":"#00008B", "medium blue":"#0000CD", "blue":"#0000FF", "royal blue":"#4169E1", "blue violet":"#8A2BE2",
    "indigo":"#4B0082", "dark slate blue":"#483D8B", "slate blue":"#6A5ACD", "medium slate blue":"#7B68EE", "medium purple":"#9370DB", "dark magenta":"#8B008B",
    "dark violet":"#9400D3", "dark orchid":"#9932CC", "medium orchid":"#BA55D3", "purple":"#800080", "thistle":"#D8BFD8", "plum":"#DDA0DD",
    "violet":"#EE82EE", "magenta":"#FF00FF", "orchid":"#DA70D6", "medium violet red":"#C71585", "pale violet red":"#DB7093", "deep pink":"#FF1493",
    "hot pink":"#FF69B4","light pink":"#FFB6C1","pink":"#FFC0CB","antique white":"#FAEBD7","beige":"#F5F5DC","bisque":"#FFE4C4",
    "blanched almond":"#FFEBCD","wheat":"#F5DEB3","corn silk":"#FFF8DC","lemon chiffon":"#FFFACD","light golden rod yellow":"#FAFAD2","light yellow":"#FFFFE0",
    "saddle brown":"#8B4513","sienna":"#A0522D","chocolate":"#D2691E","peru":"#CD853F","sandy brown":"#F4A460","burly wood":"#DEB887",
    "tan":"#D2B48C","rosy brown":"#BC8F8F","moccasin":"#FFE4B5","navajo white":"#FFDEAD","peach puff":"#FFDAB9","misty rose":"#FFE4E1",
    "lavender blush":"#FFF0F5","linen":"#FAF0E6","old lace":"#FDF5E6","papaya whip":"#FFEFD5","sea shell":"#FFF5EE","mint cream":"#F5FFFA",
    "slate gray":"#708090","light slate gray":"#778899", "light steel blue":"#B0C4DE","lavender":"#E6E6FA","floral white":"#FFFAF0","alice blue":"#F0F8FF",
    "ghost white":"#F8F8FF","honeydew":"#F0FFF0","ivory":"#FFFFF0","azure":"#F0FFFF","snow":"#FFFAFA","black":"#000000",
    "dim gray":"#696969","gray":"#808080","dark gray":"#A9A9A9","silver":"#C0C0C0","light gray":"#D3D3D3","gainsboro":"#DCDCDC",
    "white smoke":"#F5F5F5","white":"#FFFFFF"}
    return xcolor[key]

def plot_stn(DATA, fcst_hrs, accums, accum_strs, dt_utc_now, stn_name, COLORS, font_text, fig_keys, png_stn_name):
#     PDT    
#     utc_converter = pytz.utc
#     dt_utc_now = utc_converter.localize(dt_utc_now)
#     dt_utc_now = dt_utc_now.astimezone(pytz.timezone("America/Los_Angeles"))
    
    dt_utc_now = dt_utc_now-timedelta(hours=8)
    dt_ini_str = datetime.strftime(dt_utc_now, '%Y-%h-%d %H00 PST (-8 UTC)')
    
    fig = plt.figure(figsize=(18, 42), dpi=70)#fig_keys['dpi']
    gs = gridspec.GridSpec(7, 1, height_ratios=[1, 1, 1, 1, 1, 1, 1])
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[2, 0])
    ax4 = plt.subplot(gs[3, 0])
    ax5 = plt.subplot(gs[4, 0])
    ax6 = plt.subplot(gs[5, 0])
    ax7 = plt.subplot(gs[6, 0])
    
    AX = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
    
    plt.subplots_adjust(0, 0, 1, 1, hspace=0.4, wspace=0)
    
    # ===== Text ===== #
    
    handle_title = []
    AX_text_order = []
    gap = 0.14285+0.00625 #0.176
    for i in range(7):
        ax_temp = fig.add_axes([0.01, 0.99-i*gap, 0.375, 0.006], facecolor='w')
        [j.set_linewidth(0.0) for j in ax_temp.spines.values()]
        
        ax_temp.tick_params(axis='both', left=False, top=False, right=False, bottom=False, 
                            labelleft=False, labeltop=False, labelright=False, labelbottom=False)
        AX_text_order.append(ax_temp)

    AX_text = AX_text_order
    # ===== Bars ===== #
    
    
    STEPs = [3, 24]
    # ---------- 3 hrs ---------- #
    i = 0
    ax = AX[i]
    accum = accums[i]
    accum_str = accum_strs[i]
    fcst_hr = fcst_hrs[i]
        
    ax = ax_decorate_box(ax, font_text, left_flag=True, bottom_flag=True)
    ax.grid(linestyle=':', color='k'); ax.xaxis.grid(False)
        
    dt_valid = []
    for hrs in fcst_hr:
        dt_valid.append(datetime.strftime(dt_utc_now+timedelta(hours=hrs), '%H00 %d %h'))
    
    for l, temp_h in enumerate(fcst_hr):
        # ((bottom left), width, height)
        ax.add_patch(patches.Rectangle((temp_h-STEPs[i], DATA['CNN_{}_min'.format(accum)][l]), STEPs[i], DATA['CNN_{}_max'.format(accum)][l], facecolor=COLORS[0], edgecolor='0.5', linewidth=1.5))
        ax.add_patch(patches.Rectangle((temp_h-STEPs[i], DATA['CNN_{}_P10'.format(accum)][l]), STEPs[i], DATA['CNN_{}_P90'.format(accum)][l], facecolor=COLORS[1], edgecolor='0.5', linewidth=1.5))
        ax.add_patch(patches.Rectangle((temp_h-STEPs[i], DATA['CNN_{}_P25'.format(accum)][l]), STEPs[i], DATA['CNN_{}_P75'.format(accum)][l], facecolor=COLORS[2], edgecolor='0.5', linewidth=1.5))
        ax.hlines(y=DATA['CNN_{}_P50'.format(accum)][l], xmin=temp_h, xmax=temp_h, color='k', linestyle='--', linewidth=2.5);
            
    ax.set_xticks(fcst_hr);
    ax.set_xticklabels(dt_valid, rotation='90');
    ax.set_xlim([fcst_hrs[0][0]-6, fcst_hrs[0][-1]+3]);
    
    # ---------- 24 hrs ---------- #
    i = 1
    ax = AX[i]
    accum = accums[i]
    accum_str = accum_strs[i]
    fcst_hr = fcst_hrs[i]
        
    ax = ax_decorate_box(ax, font_text, left_flag=True, bottom_flag=True)
    ax.grid(linestyle=':', color='k'); ax.xaxis.grid(False)
        
    dt_valid = []
    for hrs in fcst_hr:
        dt_valid.append(datetime.strftime(dt_utc_now+timedelta(hours=hrs-STEPs[i]), '%d %h'))
        
    for l, temp_h in enumerate(fcst_hr):
        # ((bottom left), width, height)
        ax.add_patch(patches.Rectangle((temp_h-STEPs[i], DATA['CNN_{}_min'.format(accum)][l]), STEPs[i], DATA['CNN_{}_max'.format(accum)][l], facecolor=COLORS[0], edgecolor='0.5', linewidth=1.5))
        ax.add_patch(patches.Rectangle((temp_h-STEPs[i], DATA['CNN_{}_P10'.format(accum)][l]), STEPs[i], DATA['CNN_{}_P90'.format(accum)][l], facecolor=COLORS[1], edgecolor='0.5', linewidth=1.5))
        ax.add_patch(patches.Rectangle((temp_h-STEPs[i], DATA['CNN_{}_P25'.format(accum)][l]), STEPs[i], DATA['CNN_{}_P75'.format(accum)][l], facecolor=COLORS[2], edgecolor='0.5', linewidth=1.5))
        ax.hlines(y=DATA['CNN_{}_P50'.format(accum)][l], xmin=temp_h-STEPs[i], xmax=temp_h, color='k', linestyle='--', linewidth=2.5);
            
    ax.set_xticks(fcst_hr-0.5*STEPs[i]);
    ax.set_xticklabels(dt_valid, rotation='90');
    ax.set_xlim([fcst_hrs[0][0]-6, fcst_hrs[0][-1]+3]);
        
    # ===== Lines ===== #
    for i in range(2, 7, 1):
        ax = AX[i]
        accum = accums[i]
        accum_str = accum_strs[i]
        fcst_hr = fcst_hrs[i]
        
        ax = ax_decorate_box(ax, font_text, left_flag=True, bottom_flag=True)
        ax.grid(linestyle=':', color='k'); ax.xaxis.grid(True)
        
        dt_valid = []
        for hrs in fcst_hr:
            dt_valid.append(datetime.strftime(dt_utc_now+timedelta(hours=hrs), '%H00 %d %h'))
        
        ax.fill_between(fcst_hr, DATA['CNN_{}_min'.format(accum)], DATA['CNN_{}_max'.format(accum)], color=COLORS[0])
        ax.fill_between(fcst_hr, DATA['CNN_{}_P10'.format(accum)], DATA['CNN_{}_P90'.format(accum)], color=COLORS[1])
        ax.fill_between(fcst_hr, DATA['CNN_{}_P25'.format(accum)], DATA['CNN_{}_P75'.format(accum)], color=COLORS[2])
        handle_line = ax.plot(fcst_hr, DATA['CNN_{}_P50'.format(accum)], 'k--', lw=2.5, label='50th');
        
        ax.set_xticks(fcst_hr);
        ax.set_xticklabels(dt_valid, rotation='90');
        ax.set_xlim([fcst_hrs[0][0]-6, fcst_hrs[0][-1]+3]);
    
    # ===== Ticks & Legend ===== #
    for i in range(7):
        ax = AX[i]
        accum = accums[i]
        accum_str = accum_strs[i]
        fcst_hr = fcst_hrs[i]
        ax.set_ylim([0, 0.5+np.nanmax(DATA['CNN_{}_max'.format(accum)])*1.25])
        #ax.set_xlabel('Valid time (accumulated hours backward)', fontsize=font_text);
        ax.set_ylabel('[mm]', fontsize=font_text);
        
    for i in [0, 1, -1]:
        
        handle_title += string_partial_format(fig, AX_text[i], 0, 1, 'left', 'top',
                                              ['Post-processed GEFS ', accum_strs[i], ' total precipitation ensemble at ', stn_name, '. Initialization time: ', dt_ini_str],
                                              ['k',]*6, [font_text+2,]*6, ['normal', 'bold', 'normal', 'bold', 'normal', 'normal'])
        
    for i in range(2, 6, 1):
        
        handle_title += string_partial_format(fig, AX_text[i], 0, 1, 'left', 'top',
                                              ['Post-processed GEFS accumulated total precipitation ensemble ', accum_strs[i], ', at ', stn_name, '. Initialization time: ', dt_ini_str],
                                              ['k',]*6, [font_text+2,]*6, ['normal', 'bold', 'normal', 'bold', 'normal', 'normal'])
        
    for handle in handle_title:
        handle.set_bbox(dict(facecolor='w', edgecolor='none', pad=0.0, zorder=6))
        
    handle_legneds = []
    handle_legneds.append(patches.Patch(facecolor=COLORS[0], edgecolor='k', linewidth=0, label='min - max'))
    handle_legneds.append(patches.Patch(facecolor=COLORS[1], edgecolor='k', linewidth=0, label='10th - 90th'))
    handle_legneds.append(patches.Patch(facecolor=COLORS[2], edgecolor='k', linewidth=0, label='25th - 75th'))


    
    ax_lg2 = fig.add_axes([1.0, 1.0-0.04, 0.11, 0.04])
    ax_lg2.set_axis_off()
    LG2 = ax_lg2.legend(handles=handle_legneds+handle_line, bbox_to_anchor=(0.0, 0.5), ncol=1, loc=6, prop={'size':font_text}, fancybox=False);
    LG2.get_frame().set_facecolor('w')
    LG2.get_frame().set_linewidth(2.0)
    LG2.get_frame().set_alpha(1.0)
    
    
    output_dir = png_stn_name['dir']+png_stn_name['base_'].format(stn_name)
    output_dir_now =  datetime.strftime(dt_utc_now, output_dir)
    
    print(output_dir_now)
    fig.savefig(output_dir_now, format='png', **fig_keys)
        
        