
'''
This Python script contains functions used my the "main_nonDL.py" script.


'''

import numpy as np
import numba as nb

@nb.njit()
def shift_one(arr):
    """
    !!! interal function, do not touch.
    """
    arr[1:] = arr[:-1]
    return arr

@nb.njit()
def search_nearby_days(day, window=30, leap_year=False):
    """"
    Given a day as "the n-th day of a year", the program returns an array of neighboured 
    days based on a window length.
    
    Author:
    Yingkai Sha <yingkaisha@gmail.com>
    
    Parameters
    ----------
    day : int
        Date of the year. e.g., the 1st day of the year is day_ini=0
    window : int
        the length of day window
    
    Returns
    -------
    ind_date : array
        an array of days
        
    """
    if leap_year:
        N_days = 366
        ind_date = np.zeros((N_days,))
    else:
        N_days = 365
        ind_date = np.zeros((N_days,))
        
    ind_right = day+window+1
    ind_left = day-window
    
    if ind_left >= 0 and ind_right <= N_days:
        ind_date[day-window:day+window+1] = True
    elif ind_left < 0:
        ind_date[0:day+window+1] = True
        ind_date[day-window:] = True
    else:
        ind_diff = day+window+1-N_days
        ind_date[day-window:] = True
        ind_date[:ind_diff] = True
    return ind_date


@nb.njit()
def analog_search_SL_single_day(day_ini, year_analog, fcst_apcp, fcst_pwat, APCP, PWAT, ERA5, EN, SL_xy, flag_leap_year):
    """
    The program performs univariate Analog Ensemble (AnEn) post-processing 
    on a single initialization time and a single forecast lead time.
    The AnEn considers in-place analogs and supplemental locations.
    
    Author: 
    Yingkai Sha <yingkaisha@gmail.com> 
    
    Reference: 
    Hamill, T.M., Scheuerer, M. and Bates, G.T., 2015. Analog probabilistic precipitation forecasts using 
    GEFS reforecasts and climatology-calibrated precipitation analyses. Monthly Weather Review, 143(8), pp.3300-3309.

    Parameters
    ----------
    day_ini : int
        Date of the year. e.g., the 1st day of the year is day_ini=0
    year_analog : array
        A array of (historical) years for seaching analog days. e.g. np.array([2000, 2001,])
    fcst_apcp : array
        Forecasted total precipitation (APCP). shape=(grids,)
    fcst_pwat : array
        Forecasted precipitable water (PWAT). shape=(grids,)
    APCP : tuple
        A tuple of historical APCP forecasts. The indexing of tuple shall be corresponded to the `year_analog` indexing.
        Each tuple element is an array; its shape=(grids_for_anen_search,)
        e.g. 
        APCP=(apcp_2000_array, apcp_2001_array,) 
        # for year_analog=np.array([2000, 2001])
        # apcp_2000_array.shape = (grids_for_anen_search,)
        # !!! Note: grids_for_anen_search != grids; it contains forecast grids and supplemental locations.
    PWAT : tuple
        The same as of APCP, but for PWAT forecasts.
    ERA5 : tuple
        The same as of APCP, but for historical analysis or observations
    EN : int
        Number of AnEn members
    SL_xy : tuple
        A tuple of Supplemental Locations (SL).
        Each tuple element is an array that contains SL indices for all grids and one month
        e.g.
        SL_xy=(SL_Jan, SL_Feb, SL_Mar,)
        # SL_Jan.shape = (grids, number_of_SLs)
    flag_leap_year : bool
        If day_ini is in a leap year; then flag_leap_year=True 
        False otherwise

    Returns
    -------
    AnEn : array
        AnEn members. shape=(grids, number_of_members)
        
    """
    
    # params that can be adjusted
    #EN = 75 # number of AnEn members
    N_SL = 20 # number of suplemental locations
    N_grids = fcst_apcp.shape[0] # number of grid points
    window_day = 30 # time window (by days) of the analog search 
    shape_ravel = (2*window_day+1, N_SL)
    
    # output
    AnEn = np.empty((N_grids, EN))
    
    # allocate single day, grid, and year
    day_n = np.empty((EN,), np.int_)
    ind_n = np.empty((EN,), np.int_)
    year_n = np.empty((EN,), np.int_)
    record_n = np.ones((EN,))
    day_per_sl = 0
    #
    # datetime related variables
    day_365 = np.arange(365)
    day_366 = np.arange(366)
    
    if flag_leap_year:
        mon_ind_all = [ 
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
        1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  
        3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  
        6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6, 
        7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  
        8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  
        9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
        11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11]
    else:
        mon_ind_all = [ 
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
        1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  
        3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  
        6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6, 
        7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  
        8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  
        9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
        11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11]
    
    # selecting initialization month, and corresponded SL for all grids
    ind_month = mon_ind_all[day_ini]
    SL_mon = SL_xy[ind_month]
        
    # loop over grid points
    for n in range(N_grids):
        # fcst values
        apcp_new = fcst_apcp[n,]
        pwat_new = fcst_pwat[n,]
            
        # initialize analog selction records
        record_n[:] = 9999
            
        # loop over reforecast
        for year_ind, year_ana in enumerate(year_analog):
                
            # the 91-day window of analog search
            if year_ana%4 == 0:
                flag_analog_days = search_nearby_days(day_ini, window=window_day, leap_year=True)
                # True/Flase flag to indices
                day_base = day_366[flag_analog_days==1]
                
            else:
                flag_analog_days = search_nearby_days(day_ini, window=window_day, leap_year=False)
                day_base = day_365[flag_analog_days==1]
                
            # sumplemental locations in "month" and "grid n"
            inds_to_inds = SL_mon[n, :] # the input are "indexed" vectors, here create an index for the indexed, so it's named as inds_to_inds 
            
            # loop over a single SL (including itself)
            for s in range(shape_ravel[1]):
                
                # sl inds to actual inds
                ind_real = inds_to_inds[s]
                    
                # in place analog search
                if s == 0:
                    # loop over the time window
                    for d in range(shape_ravel[0]):
                            
                        day_real = int(day_base[d])
                        apcp_old = APCP[year_ind][day_real, ind_real]
                        pwat_old = PWAT[year_ind][day_real, ind_real]

                        # analog criteria
                        record_temp = 0.76*np.abs(apcp_old - apcp_new) + 0.24*np.abs(pwat_old - pwat_new)
                        
                        # if in place analog hit the new record
                        if record_temp < record_n[-1]:
                            # searchosrt positions
                            ind_analog = np.searchsorted(record_n, record_temp)
                            # shift one from the position to free space
                            day_n[ind_analog:] = shift_one(day_n[ind_analog:])
                            ind_n[ind_analog:] = shift_one(ind_n[ind_analog:])
                            year_n[ind_analog:] = shift_one(year_n[ind_analog:])
                            record_n[ind_analog:] = shift_one(record_n[ind_analog:])
                            # insert
                            day_n[ind_analog] = day_real
                            ind_n[ind_analog] = ind_real
                            year_n[ind_analog] = year_ind
                            record_n[ind_analog] = record_temp
                    
                # SL analog search (one analog per year)
                else:
                    record_per_sl = 9999
                    # loop over the time window
                    for d in range(shape_ravel[0]):
                            
                        day_real = int(day_base[d])
                        apcp_old = APCP[year_ind][day_real, ind_real]
                        pwat_old = PWAT[year_ind][day_real, ind_real]
                        
                        # analog criteria of 0.7*APCP + 0.3*PWAT
                        record_temp = 0.76*np.abs(apcp_old - apcp_new) + 0.24*np.abs(pwat_old - pwat_new)
                        # update the best analog of this sl
                        if record_temp < record_per_sl:
                            record_per_sl = record_temp
                            day_per_sl = day_real
                        
                    # the best SL analog is allowed to participate
                    if record_per_sl < record_n[-1]:
                        # searchosrt positions
                        ind_analog = np.searchsorted(record_n, record_per_sl)

                        # shift one from the position to free space
                        day_n[ind_analog:] = shift_one(day_n[ind_analog:])
                        ind_n[ind_analog:] = shift_one(ind_n[ind_analog:])
                        year_n[ind_analog:] = shift_one(year_n[ind_analog:])
                        record_n[ind_analog:] = shift_one(record_n[ind_analog:])

                        # insert
                        day_n[ind_analog] = day_per_sl
                        ind_n[ind_analog] = ind_real
                        year_n[ind_analog] = year_ind
                        record_n[ind_analog] = record_per_sl
                            
        # back to the grid point loop
        # assigning ERA5 based on the (multi-year) reforecast search
        for en in range(EN):
            AnEn[n, en] = ERA5[year_n[en]][day_n[en], ind_n[en]]
    return AnEn

@nb.njit()
def CDF_estimate(X):
    """
    The program computes CDFs from ensemble members or equally likely historical sequences.
    CDFs will be computed based on a fixed set of probability bins
    
    Author: 
    Yingkai Sha <yingkaisha@gmail.com> 

    Parameters
    ----------
    X : array
        gridded data. shape=(number_of_members, lead_time, grids)

    Returns
    -------
    CDF : array
        A set of CDFs. 
        
    """
    q_bins = np.array([0.25, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99])  
    _, N_fcst, N_grids = X.shape
    CDF = np.empty((9, N_fcst, N_grids))
    
    for lead in range(N_fcst):
        for n in range(N_grids):
            CDF[:, lead, n] = np.quantile(X[:, lead, n], q_bins)
    return CDF

@nb.njit()
def total_divergence(CDF1, CDF2):
    """
    The program computes the total divergence of two sets of CDFs
    
    Author: 
    Yingkai Sha <yingkaisha@gmail.com> 
    
    Reference: 
    Scheuerer, M. and Hamill, T.M., 2018. Generating calibrated ensembles of physically 
    realistic, high-resolution precipitation forecast fields based on GEFS model output. 
    Journal of Hydrometeorology, 19(10), pp.1651-1670.

    Parameters
    ----------
    CDF1 : array
        A set of CDFs. shape=(probability_bins, lead_time, grids)
    CDF2 : array
        The other set of CDFs. shape=(probability_bins, lead_time, grids)

    Returns
    -------
    TD : float
        Total divergence. 
        
    """
    
    _, N_fcst, N_grids = CDF1.shape
    TD = 0
    for lead in range(N_fcst):
        for n in range(N_grids):
            TD += np.sum(np.abs(CDF1[:, lead, n] - CDF2[:, lead, n]))
    return TD

@nb.njit()
def MDSS_main(ERA5_mdss, AnEn_out, factor=5, max_trial=5000):
    """
    The program takes historical weather sequences and post-processed univariate members as inputs;
    it returns indices of the historical sequences. The indices are used to subset the historical 
    sequences, producing dependence templates.
    
    Author: 
    Yingkai Sha <yingkaisha@gmail.com> 
    
    Reference: 
    Scheuerer, M. and Hamill, T.M., 2018. Generating calibrated ensembles of physically 
    realistic, high-resolution precipitation forecast fields based on GEFS model output. 
    Journal of Hydrometeorology, 19(10), pp.1651-1670.

    Parameters
    ----------
    ERA5_mdss : array
        Historical weather sequences. shape=(number_of_sequence, lead_time, grids)
    AnEn_out : array
        Univeriate ensemble members. shape=(number_of_members, lead_time, grids)
    factor : float
        A heuristic factor that controls the speed of MDSS convergence. 
        Larger means lower convergence with better quality.
    max_trial : int
        Maximum number of heuristic trials
        If MDSS iterations exceed the max_trial, a random selection will be performed.

    Returns
    -------
    flag_clean : array
        Indicies of historical sequences.
    count_trial : int
        Number of MDSS iterations 
        
    """
    
    # Collect shape info
    ## K: number of ensemble members
    ## N: number of sequences
    K, N_lead_, N_grid_ = AnEn_out.shape
    N, N_leads, N_grids = ERA5_mdss.shape
    
    # assertion: lead_time and grids must match.
    assert N_lead_ == N_leads
    assert N_grid_ == N_grids
    
    # Compute CDFs
    CDF_fcst = CDF_estimate(AnEn_out)
    CDF_ERA5 = CDF_estimate(ERA5_mdss)

    # initial total divergence
    record = total_divergence(CDF_ERA5, CDF_fcst)

    # allocate all available indices
    ind_pick = np.arange(N, dtype=np.int_)
    
    # allocate the output indices
    flag_clean = np.ones((N,), dtype=np.bool_)
    
    # allocation indices for each iteration
    flag_trial = np.copy(flag_clean) # temporal trial pick
    
    # Inner params
    N_0 = N                      # Current number of sequences 
    flag_single = False          # Bool flag: True means a single heuristic drop; False means a batch drop
    bad_list = np.array([9999,]) # Indices discarded by the MDSS iterations
    count_trial = 0              # Iteration counter

    # MDSS main iteration loop
    while N_0 > K and count_trial < max_trial:
        
        # all available sequence based on previous iterations
        ind_candidate = ind_pick[flag_clean]

        # Estimate the number to dropout
        size_ = int((N_0-K)/(factor))

        # If size > 5, perform batch heuristic dropout
        if size_ > 5:
            # generate dropout indices
            ind_ = np.random.choice(np.arange(N_0), size_, replace=False)

            # update droput indices to the trial record
            for i in ind_:
                flag_trial[ind_candidate[i]] = False

        # If size_ <= 5, perform single heuristic dropout
        else:
            flag_single = True

            # generate the index to dropout
            j = np.random.choice(ind_candidate, size=1)
            
            # If the single index has been picked before, pick another one
            while np.any(bad_list == j):
                j = np.random.choice(ind_candidate, size=1)
            
            # update the index to the trial record
            flag_trial[j] = False
        
        # Compute the CDF of selected historical sequence after dropout
        ERA5_sub = ERA5_mdss[flag_trial, ...]
        CDF_ERA5 = CDF_estimate(ERA5_sub)
        
        # Compute the current total divergence
        record_temp = total_divergence(CDF_ERA5, CDF_fcst)
        
        # If divergence is reduced, update to the final output (clean)
        if record_temp < record:
            record = record_temp
            
            if flag_single:
                flag_clean[j] = False
            else:
                flag_clean = np.copy(flag_trial)

        # if divergence not reduced, revert to the previous clean record
        else:
            if flag_single:                
                flag_trial[j] = True
                # if it is a single dropout, add this index to the bad_list
                bad_list = np.append(bad_list, j)
            else:
                flag_trial = np.copy(flag_clean)
                
        # Update counters at the end of the loop
        N_0 = np.sum(flag_clean)
        count_trial += 1
    
    # Return required number of indices regardless of how many has been dropout
    # Return count_trial to check if the MDSS is converged.
    return flag_clean, count_trial

@nb.njit()
def schaake_shuffle(fcst, traj):
    """
    The program shuffles univariate ensemble members based on historical weather sequences;
    it produces spatiotemorally consistent ensemble members as output
    
    Author: 
    Yingkai Sha <yingkaisha@gmail.com> 
    
    Reference: 
    Clark, M., Gangopadhyay, S., Hay, L., Rajagopalan, B. and Wilby, R., 2004. 
    The Schaake shuffle: A method for reconstructing space–time variability in 
    forecasted precipitation and temperature fields. 
    Journal of Hydrometeorology, 5(1), pp.243-262.

    Parameters
    ----------
    fcst : array
        Univeriate ensemble members. shape=(number_of_members, lead_time, grids)
    AnEn_out : array
        Historical weather sequences. shape=(number_of_sequence, lead_time, grids)
    
    Returns
    -------
    output : array
        spatiotemorally consistent ensemble members 
        
    """
    num_traj, N_lead, N_grids = traj.shape
    
    output = np.empty((num_traj, N_lead, N_grids))
    
    for l in range(N_lead):
        for n in range(N_grids):
            
            temp_traj = traj[:, l, n]
            temp_fcst = fcst[:, l, n]
            
            reverse_b_func = np.searchsorted(np.sort(temp_traj), temp_traj)
            
            output[:, l, n] = np.sort(temp_fcst)[reverse_b_func]
    return output


@nb.njit()
def dressing_norm_flat(fcst, folds, k1, k2):
    '''
    AnEn dressing
    '''
    N_lead, N_grids, EN = fcst.shape
    out = np.empty((N_lead, N_grids, EN*folds))
    
    for en in range(EN):
        for lead in range(N_lead):
            for ix in range(N_grids):
                for f in range(folds):
                    base_ = fcst[lead, ix, en]
                    std_ = k1 + k2*base_
                    dress_ = np.random.normal(loc=0.0, scale=std_)
                    out[lead, ix, en*folds+f] = fcst[lead, ix, en] + dress_          
    return out

@nb.njit()
def dressing_norm_2d(fcst, land_mask, folds, k1, k2):
    '''
    AnEn dressing
    '''
    EN, N_lead, Nx, Ny = fcst.shape
    out = np.empty((EN*folds, N_lead, Nx, Ny))
    
    for en in range(EN):
        for lead in range(N_lead):
            for ix in range(Nx):
                for iy in range(Ny):
                    if land_mask[ix, iy]:
                        for f in range(folds):
                            base_ = fcst[en, lead, ix, iy]
                            std_ = k1 + k2*base_
                            dress_ = np.random.normal(loc=0.0, scale=std_)
                            out[en*folds+f, lead, ix, iy] = fcst[en, lead, ix, iy] + dress_          
    return out



def cnn_precip_fix(data):
    '''
    Precipitation values predicted by the CNNs are
    not non-negative. Fixing this issue with logical_and.
    '''
    
    flag_0 = data<0.5
    flag_1 = np.logical_and(data<0.6, data>=0.5)
    flag_2 = np.logical_and(data<0.7, data>=0.6)
    flag_3 = np.logical_and(data<0.8, data>=0.7)
    flag_4 = np.logical_and(data<0.9, data>=0.8)

    data[flag_0] = 0.0
    data[flag_1] = 0.1
    data[flag_2] = 0.3
    data[flag_3] = 0.5
    data[flag_4] = 0.7
    data[data>40] = 40
    return data