# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 14:27:46 2020

@author: IsoldeGlissenaar
"""

import numpy as np
import xarray as xr
from datetime import date
from func_reproj import reproject
from func_scatteredInterpolant import scatteredInterpolant_CS
import scipy.io
import h5py
import sit_functions

def snowdepth_warren(data, petty=False):
    # Function snowdepth_warren(data, petty)
    # Input: data, petty
    # Output: is_snow_depth_w99
    #
    # get Warren '99 snow depth and interpolates to freeboard observations
    #
    # given coordinates and time of satellite altimetry observations, retrieves
    # Warren '99 climatology snow depth and interpolates to location of 
    # altimetry observation. If petty=True, also returns redistributed snow
    # depth according to the piecewise function in Petty et al. (2020).
    # Used to find snow depth at satellite altimetry observation locations.
    
    w99 = np.load('C:/Users/zq19140/OneDrive - University of Bristol/Documents/Scripts/Data/Warren_snowdepth/mean_warren99_snowdepth_20102018_'+str("{:02d}".format(data.month[0]))+'.npy')
    lat_w99 = np.load('C:/Users/zq19140/OneDrive - University of Bristol/Documents/Scripts/Data/Warren_snowdepth/warren_lat.npy')
    lon_w99 = np.load('C:/Users/zq19140/OneDrive - University of Bristol/Documents/Scripts/Data/Warren_snowdepth/warren_lon.npy')
    snow_depth_w99 = np.transpose(np.array([lat_w99.flatten(), lon_w99.flatten(), w99.flatten()]))
    snow_depth_w99[:,0], snow_depth_w99[:,1] = reproject(snow_depth_w99[:,0], snow_depth_w99[:,1])
    is_snow_depth_w99 = scatteredInterpolant_CS(snow_depth_w99, data)
    
    if petty:
        coord = np.transpose(np.array([data.lat.values, data.lon.values]))
        U = np.unique(data.day, axis=0)
        is_snow_depth_w99_petty = np.zeros(len(data))
        #For each day in the freeboard values
        for i in range(0,len(U)):
            #Select freeboard data for that day
            idx = (data.day==U[i])
            fb_day1 = data.freeboard.values[idx]
        
            #Find fb's per snow depth gridcell
            no_cell_snowdepth = np.arange(0,len(snow_depth_w99),1)
            snowdepthcell_fb = scipy.interpolate.griddata(snow_depth_w99[:,0:2], no_cell_snowdepth, coord[idx,:], 'nearest')
            snowdepth_cells = np.unique(snowdepthcell_fb)
            
            #Loop through large scale snow depth grid cells
            fb_sd_grid = np.zeros(len(fb_day1))
            for n in range(len(snowdepth_cells)):
                loc_snowdepths = (snowdepthcell_fb==snowdepth_cells[n])
                sd_cell = snow_depth_w99[snowdepth_cells,2]
                hsl = sd_cell[n]
                fb_sd_grid[loc_snowdepths] = sit_functions.petty(fb_day1[loc_snowdepths], 
                                                                           hsl)

            is_snow_depth_w99_petty[idx] = fb_sd_grid
        return is_snow_depth_w99, is_snow_depth_w99_petty
    
    return is_snow_depth_w99




def snowdepth_pmw(data, petty=False):
    # Function snowdepth_pmw(data, petty)
    # Input: data, petty
    # Output: is_snow_depth_pmw
    #
    # get PMW snow depth and interpolates to freeboard observations
    #
    # given coordinates and time of satellite altimetry observations, retrieves
    # passive microwave snow depth and interpolates to location of 
    # altimetry observation. If petty=True, also returns redistributed snow
    # depth according to the piecewise function in Petty et al. (2020).
    # Used to find snow depth at satellite altimetry observation locations.
    
    sd_folder='C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/DMSP SSMIS 25 km Snow Depth Eastern Canadian Arctic/'
    U = np.unique(data.day, axis=0)
    
    is_snow_depth_pmw = np.zeros(len(data))
    snow_depth_petty = np.zeros(len(data))
    #For each day in the freeboard values
    for i in range(0,len(U)):
        #Open DMSP SSMIS 25 km snow depth
        sd_id= str("{:.0f}".format(data.year[0]))+str("{:02d}".format(data.month[0]))+str("{:02d}".format(int(U[i])))
        sd_filename = sd_folder + 'snowdepth_fromTB_' + sd_id + '.mat'
        try:
            snow_depth_Tb = scipy.io.loadmat(sd_filename)['snow_depth_Tb']
        except FileNotFoundError:
            print(f'PMW snow depth on {data.year[0]}/{data.month[0]}/{int(U[i])} not available')
            sd_id= str("{:.0f}".format(data.year[0]))+str("{:02d}".format(data.month[0]))+str("{:02d}".format(int(U[i]+1)))
            sd_filename = sd_folder + 'snowdepth_fromTB_' + sd_id + '.mat'
            snow_depth_Tb = scipy.io.loadmat(sd_filename)['snow_depth_Tb']
            
        snow_depth_Tb[:,0], snow_depth_Tb[:,1] = reproject(snow_depth_Tb[:,0], snow_depth_Tb[:,1])
       
        #Filter for GRV
        ids = snow_depth_Tb[:,3]>-0.045
        snow_depth_Tb = snow_depth_Tb[ids,:]
        idx = (data.day==U[i])
        
        if petty:
            coord = np.transpose(np.array([data.lat.values, data.lon.values]))
            fb_day1 = data.freeboard.values[idx]
            #Find fb's per snow depth gridcell
            no_cell_snowdepth = np.arange(0,len(snow_depth_Tb),1)
            snowdepthcell_fb = scipy.interpolate.griddata(snow_depth_Tb[:,0:2], no_cell_snowdepth, coord[idx,:], 'nearest')
            snowdepth_cells = np.unique(snowdepthcell_fb)
            
            #Loop through large scale snow depth grid cells
            fb_sd_grid = np.zeros(len(fb_day1))
            for n in range(len(snowdepth_cells)):
                loc_snowdepths = (snowdepthcell_fb==snowdepth_cells[n])
                sd_cell = snow_depth_Tb[snowdepth_cells,2]
                hsl = sd_cell[n]
                fb_sd_grid[loc_snowdepths] = sit_functions.petty(fb_day1[loc_snowdepths], 
                                                                            hsl,)
            snow_depth_petty[idx] = fb_sd_grid
        
        is_snow_depth_pmw[idx] = scatteredInterpolant_CS(snow_depth_Tb, data[idx])
    if petty:
        return is_snow_depth_pmw, snow_depth_petty
    else:
        return is_snow_depth_pmw

    

def snowdepth_snowmodel(data, petty=False):
    # Function snowdepth_snowmodel(data, petty)
    # Input: data, petty
    # Output: is_snow_depth_snowmodel
    #
    # get SnowModel-LG snow depth and interpolates to freeboard observations
    #
    # given coordinates and time of satellite altimetry observations, retrieves
    # SnowModel-LG (Liston et al., 2020) snow depth and interpolates to location 
    # of altimetry observation. If petty=True, also returns redistributed snow
    # depth according to the piecewise function in Petty et al. (2020).
    # Used to find snow depth at satellite altimetry observation locations.
    
    if data.year[0]<=2018:
        sd = xr.open_dataset('C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/SnowModel/1980-2018/snowdepth_merra.nc')
    if data.year[0]>2018:
        sd = xr.open_dataset('C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/SnowModel/SM_merra2_ease_01Aug2018-31Jul2020.nc')
    is_snow_depth_snowmodel = np.zeros(len(data))
    month2 = np.array(date(int(data.year.values[0]), int(data.month.values[0]), 1), dtype='datetime64[M]')
    t = np.where(sd.time==month2)[0][0]
    snowdepth = np.transpose(np.array([np.reshape(sd.lat.values, (361*361)), 
                             np.reshape(sd.lon.values, (361*361)), 
                             np.reshape(sd.snowdepth[:,:,t].values, (361*361))]))
    # snowdepth=snowdepth[snowdepth[:,2]>0,:]
    snowdepth[:,2][snowdepth[:,2]==0]=np.nan  
    snowdepth[:,0], snowdepth[:,1] = reproject(snowdepth[:,0], snowdepth[:,1])
    is_snow_depth_snowmodel[:] = scatteredInterpolant_CS(snowdepth, data)
    coord = np.transpose(np.array([data.lat.values, data.lon.values]))
    
    if petty:
        #Find fb's per snow depth gridcell
        no_cell_snowdepth = np.arange(0,len(snowdepth),1)
        snowdepthcell_fb = scipy.interpolate.griddata(snowdepth[:,0:2], no_cell_snowdepth, coord, 'nearest')
        snowdepth_cells = np.unique(snowdepthcell_fb)
        
        #Redistribute snow according to Petty et al. (in review)
        is_snow_depth_petty = np.zeros(len(data))
        for n in range(len(snowdepth_cells)):
            loc_snowdepths = (snowdepthcell_fb==snowdepth_cells[n])
            sd_cell = snowdepth[snowdepth_cells,2]
            hsl = sd_cell[n]
            is_snow_depth_petty[loc_snowdepths] = sit_functions.petty(data.freeboard.values[loc_snowdepths], 
                                                                           hsl)
        
    if petty:
        return is_snow_depth_snowmodel, is_snow_depth_petty
    return is_snow_depth_snowmodel



def snowdepth_nesosim(data, petty=False):
    # Function snowdepth_nesosim(data, petty)
    # Input: data, petty
    # Output: is_snow_depth_nesosim
    #
    # get NESOSIM snow depth and interpolates to freeboard observations
    #
    # given coordinates and time of satellite altimetry observations, retrieves
    # NESOSIM (Petty et al., 2020) snow depth and interpolates to location of 
    # altimetry observation. If petty=True, also returns redistributed snow
    # depth according to the piecewise function in Petty et al. (2020).
    # Used to find snow depth at satellite altimetry observation locations.
    
    nesosim = xr.open_dataset('C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/NESOSIMv11prelim/NESOSIMv11_0109'+str(data.year[0]-1)+'-3004'+str(data.year[0])+'.nc')
    month2 = int(data.month.values[0])
    nesosim_months = np.array([int(str(nesosim.day.values[d])[4:6]) for d in range(0,len(nesosim.day.values))])
    nesosim_days = np.array([int(str(nesosim.day.values[d])[6:8]) for d in range(0,len(nesosim.day.values))])
                 
    days = np.unique(data.day, axis=0)
    is_snow_depth_nesosim = np.zeros(len(data))
    snow_depth_petty = np.zeros(len(data))
    #For each day in the freeboard values
    for i in range(0,len(days)):
        loc = np.where((nesosim_months==month2)&(nesosim_days==days[i]))
        snowdepth = np.transpose(np.array([np.reshape(nesosim.latitude.values, (90*90)), 
                                           np.reshape(nesosim.longitude.values, (90*90)), 
                                           np.reshape(nesosim.snow_depth[loc[0][0],:,:].values, (90*90))]))
        snowdepth[:,0], snowdepth[:,1] = reproject(snowdepth[:,0], snowdepth[:,1])
        
        idx = (data.day==days[i])
        
        if petty:
            coord = np.transpose(np.array([data.lat.values, data.lon.values]))
            fb_day1 = data.freeboard.values[idx]
            #Find fb's per snow depth gridcell
            no_cell_snowdepth = np.arange(0,len(snowdepth),1)
            snowdepthcell_fb = scipy.interpolate.griddata(snowdepth[:,0:2], no_cell_snowdepth, coord[idx,:], 'nearest')
            snowdepth_cells = np.unique(snowdepthcell_fb)
            
            #Loop through large scale snow depth grid cells
            fb_sd_grid = np.zeros(len(fb_day1))
            for n in range(len(snowdepth_cells)):
                loc_snowdepths = (snowdepthcell_fb==snowdepth_cells[n])
                sd_cell = snowdepth[snowdepth_cells,2]
                hsl = sd_cell[n]
                fb_sd_grid[loc_snowdepths] = sit_functions.petty(fb_day1[loc_snowdepths], 
                                                                            hsl)
            snow_depth_petty[idx] = fb_sd_grid
        
        is_snow_depth_nesosim[idx] = scatteredInterpolant_CS(snowdepth, data[idx])
    if petty:
        return is_snow_depth_nesosim, snow_depth_petty
    else:
        return is_snow_depth_nesosim



def snowdepth_amsr(data, petty=False):
    # Function snowdepth_amsr(data, petty)
    # Input: data, petty
    # Output: is_snow_depth_amsr
    #
    # get AMSR snow depth and interpolates to freeboard observations
    #
    # given coordinates and time of satellite altimetry observations, retrieves
    # AMSR snow depth and interpolates to location of altimetry observation. 
    # If petty=True, also returns redistributed snow depth according to the 
    # piecewise function in Petty et al. (2020).
    # Used to find snow depth at satellite altimetry observation locations.
    
    sd_folder = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/AMSR2/'
    U = np.unique(data.day,axis=0)
    
    is_snow_depth_amsr = np.zeros(len(data))
    snow_depth_petty = np.zeros(len(data))
    #For each day in the freeboard values
    for i in range(0,len(U)):
        #Open AMSR2 12.5 km snow depth
        sd_id=str("{:.0f}".format(data.year[0]))+'_'+str("{:02d}".format(data.month[0]))+'/AMSR_U2_L3_SeaIce12km_B04_'+str("{:.0f}".format(data.year[0]))+str("{:02d}".format(data.month[0]))+str("{:02d}".format(int(U[i])))
        sd_filename = sd_folder + sd_id + '.he5'
        try:
            snow_depth_amsr = h5py.File(sd_filename,'r')
        except FileNotFoundError:
            print(f'AMSR snow depth on {data.year[0]}/{data.month[0]}/{int(U[i])} not available')
            sd_id=str("{:.0f}".format(data.year[0]))+'_'+str("{:02d}".format(data.month[0]))+'/AMSR_U2_L3_SeaIce12km_B04_'+str("{:.0f}".format(data.year[0]))+str("{:02d}".format(data.month[0]))+str("{:02d}".format(int(U[i]+1)))
            sd_filename = sd_folder + sd_id + '.he5'
            snow_depth_amsr = h5py.File(sd_filename,'r')
            
        sd = snow_depth_amsr['/HDFEOS/GRIDS/NpPolarGrid12km/Data Fields/SI_12km_NH_SNOWDEPTH_5DAY'][:,:].flatten()
        lat = snow_depth_amsr['/HDFEOS/GRIDS/NpPolarGrid12km/lat'][:,:].flatten()
        lon = snow_depth_amsr['/HDFEOS/GRIDS/NpPolarGrid12km/lon'][:,:].flatten()
            
        lat, lon = reproject(lat, lon)
        coords_amsr = np.array([lat,lon]).transpose()
        sd = sd.astype(float)
        sd[sd>100] = np.nan
        sd = sd/100 #to meters
        idx = (data.day==U[i])

        if petty:
            coord = np.transpose(np.array([data.lat.values, data.lon.values]))
            fb_day1 = data.freeboard.values[idx]
            #Find fb's per snow depth gridcell
            no_cell_snowdepth = np.arange(0,len(sd),1)
            snowdepthcell_fb = scipy.interpolate.griddata(coords_amsr, no_cell_snowdepth, coord[idx,:], 'nearest')
            snowdepth_cells = np.unique(snowdepthcell_fb)
            
            #Loop through large scale snow depth grid cells
            fb_sd_grid = np.zeros(len(fb_day1))
            for n in range(len(snowdepth_cells)):
                loc_snowdepths = (snowdepthcell_fb==snowdepth_cells[n])
                sd_cell = sd[snowdepth_cells]
                hsl = sd_cell[n]
                fb_sd_grid[loc_snowdepths] = sit_functions.petty(fb_day1[loc_snowdepths], 
                                                                            hsl)
            snow_depth_petty[idx] = fb_sd_grid
        
        amsr_array = np.array([lat, lon, sd]).transpose()
        is_snow_depth_amsr[idx] = scatteredInterpolant_CS(amsr_array, data[idx])
    if petty:
        return is_snow_depth_amsr, snow_depth_petty
    else:
        return is_snow_depth_amsr


    
def effective_snowdepth(data, snowdepth):  
    # Function effective_snowdepth(data, snowdepth)
    # Input: data, snowdepth
    # Output: eff_snowdepth
    #
    # get effective redistribution of snow depth according to KC08
    #
    # given total freeboard and snow depth, redistributes the snow depth
    # according to the sigmoidal function in Kwok & Cunningham (2008) 
    # to get the 'effective snow depth'. 
    
    def sq(x):
        return 1.462244*x**2
    
    def lin(x):
        return x-0.17
    
    def sq2(x):
        return -1.88095*(x-1.3)**2+0.998095
    
    fb_sd = data.freeboard/snowdepth
    eff_loc = (fb_sd<1.3)
    
    sd_eff_sq = (fb_sd[eff_loc]<0.3)
    sd_eff_lin = (fb_sd[eff_loc]>=0.3)&(fb_sd[eff_loc]<1.0)
    sd_eff_sq2 = (fb_sd[eff_loc]>=1.0)
    
    
    eff_sd_ratio = np.zeros(len(fb_sd[eff_loc]))
    eff_sd_ratio[sd_eff_sq] = sq(fb_sd[eff_loc][sd_eff_sq])
    eff_sd_ratio[sd_eff_lin] = lin(fb_sd[eff_loc][sd_eff_lin])
    eff_sd_ratio[sd_eff_sq2] = sq2(fb_sd[eff_loc][sd_eff_sq2])
    
    eff_sd = eff_sd_ratio*snowdepth[eff_loc]
    
    eff_snowdepth = np.copy(snowdepth)
    eff_snowdepth[eff_loc] = eff_sd
    return eff_snowdepth


def snowdensity_snowmodel(data):
    # Function snowdensity_snowmodel(data)
    # Input: data
    # Output: snowdensity_snowmodel
    #
    # get SnowModel-LG snow density and interpolates to freeboard observations
    #
    # given coordinates and time of satellite altimetry observations, retrieves
    # SnowModel-LG (Liston et al., 2020) snow density and interpolates to location 
    # of altimetry observation. 
    # Used to find snow density at satellite altimetry observation locations.
    
    if data.year[0]<=2018:
        sd = xr.open_dataset('C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/SnowModel/1980-2018/snowdepth_merra.nc')
    if data.year[0]>2018:
        sd = xr.open_dataset('C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/SnowModel/SM_merra2_ease_01Aug2018-31Jul2020.nc')
    snowdensity_snowmodel = np.zeros(len(data))
    month2 = np.array(date(int(data.year.values[0]), int(data.month.values[0]), 1), dtype='datetime64[M]')
    t = np.where(sd.time==month2)[0][0]
    snowdensity = np.transpose(np.array([np.reshape(sd.lat.values, (361*361)), 
                             np.reshape(sd.lon.values, (361*361)), 
                             np.reshape(sd.snowdensity[:,:,t].values, (361*361))]))
    snowdensity[:,2][snowdensity[:,2]==0]=330 #np.nan  
    snowdensity[:,0], snowdensity[:,1] = reproject(snowdensity[:,0], snowdensity[:,1])
    snowdensity_snowmodel[:] = scatteredInterpolant_CS(snowdensity, data)
    return snowdensity_snowmodel


def snowdensity_warren(data):
    # Function snowdensity_warren(data)
    # Input: data
    # Output: is_snow_density
    #
    # get Warren '99 snow density and interpolates to freeboard observations
    #
    # given coordinates and time of satellite altimetry observations, retrieves
    # Warren '99 snow density. 
    # Used to find snow density at satellite altimetry observation locations.
    
    snow_density = scipy.io.loadmat('C:/Users/zq19140/OneDrive - University of Bristol/Documents/Jack_onedrive/Existing Sea Ice Thickness Datasets/Landy et al 2017 ICESat GLAS13/snow_density')['snow_density']

    sep01 = date.toordinal(date(2003, 9, 1))
    datnum = np.array([[date.fromordinal(int(sep01 + (i-1))).year, 
                        date.fromordinal(int(sep01 + (i-1))).month, 
                        date.fromordinal(int(sep01 + (i-1))).day] for i in snow_density[:,0]])

    datenumber = np.zeros(len(data))
    for i in range(len(data)):
        datenumber[i] = date.toordinal(date(data.year.values[0], data.month.values[0], int(data.day.values[i])))
    
    Lia = np.isin(data.month.values[0],datnum[:,1])*np.isin(data.day.values,datnum[:,2])
    Locb = np.zeros(len(data), dtype=int)
    for i in range (len(data)):
        if Lia[i]==True:
            Locb[i] = int(np.where((datnum[:,1]==data.month.values[0]) & (datnum[:,2]==data.day.values[i]))[0])
        else:
            Locb[i] = -999
    is_snow_density = np.zeros(len(data))
    is_snow_density[Lia] = snow_density[Locb[Locb>0],1]
    return is_snow_density