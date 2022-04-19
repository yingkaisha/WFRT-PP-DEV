
import h5py
import numpy as np

# Stats modules
from scipy.spatial import cKDTree
from scipy.interpolate import interp2d
from scipy.interpolate import NearestNDInterpolator
# from sklearn.metrics import mean_squared_error
# from sklearn.metrics import mean_absolute_error

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import transforms

from cartopy.io.shapereader import Reader, natural_earth

def leap_year_checker(y):
    if y % 400 == 0:
        return True
    if y % 100 == 0:
        return False
    if y % 4 == 0:
        return True
    else:
        return False

# saving function

def save_hdf5(p_group, labels, out_dir, filename='example.hdf'):
    '''
    Save data into a signle hdf5
        - p_group: datasets combined in one tuple;
        - labels: list of strings;
        - out_dir: output path;
        - filename: example.hdf;
    **label has initial 'x' means ENCODED strings
    '''    
    name = out_dir+filename
    hdf = h5py.File(name, 'w')
    for i, label in enumerate(labels):
        if label[0] != 'x':
            hdf.create_dataset(label, data=p_group[i])
        else:
            string = p_group[i]
            hdf.create_dataset(label, (len(string), 1), 'S10', string)
    hdf.close()
    print('Save to {}'.format(name))

# def KDTree_wraper(dist_lon, dist_lat):
#     '''
#     A warper of scipy.spatial.cKDTree
#     Tree = KDTree_wraper(dist_lon, dist_lat)
#     '''
#     return cKDTree(list(zip(dist_lon.ravel(), dist_lat.ravel())))

def grid_search(xgrid, ygrid, stn_lon, stn_lat):
    '''
    kdtree-based nearest gridpoint search
    output: indices_lon, indices_lat
    '''
    gridTree = cKDTree(list(zip(xgrid.ravel(), ygrid.ravel()))) #KDTree_wraper(xgrid, ygrid)
    grid_shape = xgrid.shape
    dist, indexes = gridTree.query(list(zip(stn_lon, stn_lat)))
    return np.unravel_index(indexes, grid_shape)

def accum_slide_window(data, accum_window, output_freq, skip_start):

    inds_start = []
    inds_end = []
    
    EN, N_lead, Nx, Ny = data.shape

    N_output = (N_lead - accum_window - skip_start) // output_freq + 1

    Accum_output = np.empty((EN, N_output, Nx, Ny))

    for n in range(N_output):

        ind_start = skip_start+n*output_freq
        ind_end = ind_start+accum_window
        Accum_output[:, n, ...] = np.nansum(data[:, ind_start:ind_end, ...], axis=1)

        inds_start.append(ind_start)
        inds_end.append(ind_end)
        
    return Accum_output, inds_start, inds_end


def fillzero(arr):
    '''
    replace NaNs with zeros
    '''
    flag = np.isnan(arr)
    arr[flag] = 0.0
    return arr

def fillnan(arr):
    '''
    fill NaNs with nearest neighbour grid point val
    The first grid point (left and bottom) cannot be NaNs
    output: grid
    '''
    mask = np.isnan(arr)
    idx = np.where(~mask, np.arange(mask.shape[1]), 0)
    np.maximum.accumulate(idx, axis=1, out=idx)
    out = arr[np.arange(idx.shape[0])[:, None], idx]
    return out

def fill_coast_interp(arr, flag=False):
    '''
    Fill ocean grid points with the nearest land val
    sequence: left > top > right > bottom
    '''
    out = np.copy(arr) # copy
    # left fill
    out = np.fliplr(fillnan(np.fliplr(out)))
    # top fill
    out = np.rot90(fillnan(np.rot90(out, k=1)), k=3)
    # right fill
    out = fillnan(out)
    # bottom fill
    out = np.rot90(fillnan(np.rot90(out, k=3)), k=1)
    if type(flag) == bool:
        return out
    else:
        out[flag] = np.nan
        return out

def fill_coast(arr, flag=False):
    '''
    Fill ocean grid points with the nearest land val
    sequence: left >  bottom > right > top
    '''
    out = np.copy(arr) # copy
    # left fill
    out = np.fliplr(fillnan(np.fliplr(out)))
    # bottom fill
    out = np.rot90(fillnan(np.rot90(out, k=-1)), k=-3)
    # right fill
    out = fillnan(out)
    # top fill
    out = np.rot90(fillnan(np.rot90(out, k=-3)), k=-1)
    
    if type(flag) == bool:
        return out
    else:
        out[flag] = np.nan
        return out

def interp2d_wraper(nav_lon, nav_lat, grid_z, out_lon, out_lat, method='linear'):
    '''
    wrapper of interp2d, works for 2-d grid to grid interp.
    method = 'linear' or 'cubic'
    output: grid
    '''
    if np.sum(np.isnan(grid_z)) > 0:
        grid_z = fill_coast_interp(grid_z, np.zeros(grid_z.shape).astype(bool))
        
    interp_obj = interp2d(nav_lon[0, :], nav_lat[:, 0], grid_z, kind=method)
    return interp_obj(out_lon[0, :], out_lat[:, 0])

def nearest_wraper(nav_lon, nav_lat, grid_z, out_lon, out_lat):
    '''
    wrapper of nearest neighbour
    '''
    f = NearestNDInterpolator((nav_lon.ravel(), nav_lat.ravel()), grid_z.ravel())
    out = f((out_lon.ravel(), out_lat.ravel()))
    return out.reshape(out_lon.shape)
    
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

    
    