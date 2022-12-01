# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 13:46:07 2021

@author: IsoldeGlissenaar

Uses Envisat ESA CCI freeboard observations. Can be retrieved from:
https://catalogue.ceda.ac.uk/uuid/54e2ee0803764b4e84c906da3f16d81b
"""

import sys
sys.path.append('../functions')
from func_reproj import reproject
from func_snowdepth import snowdepth_warren, snowdepth_pmw, snowdepth_snowmodel, snowdensity_warren, snowdensity_snowmodel
import sit_functions

import xarray as xr
import numpy as np
import pandas as pd
import scipy.io
from scipy import spatial
import time
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

start = time.time()

#Select years + month
years = np.arange(2003,2013,1)
month = 3

for year in years:
    direc = docs+'SatelliteData/Envisat/swath/'+str(year)+'/'+str("{:02d}".format(month))+'/'
    envi = xr.open_mfdataset(direc+'*')
    
    envi = envi.where((envi.lon>-110)&(envi.lon<-20))
    
    # Get track number
    lat_roll = np.roll(envi.lat.values,-1)
    lon_roll = np.roll(envi.lon.values,-1)
    steps_lat = lat_roll-envi.lat.values
    a = (np.abs(steps_lat)>0.2)
    steps_lon = lon_roll-envi.lon.values
    b = (np.abs(steps_lon)>0.2)
    steps = (a|b)
    track = np.zeros(len(envi.lat.values))
    p=1
    for i in range(len(envi.lat.values)):
        if steps[i]==True:
            p=p+1
        track[i] = p
      
    
    # Create DataFrame
    envi_time = envi.time.values[~np.isnan(envi.freeboard)]
    envi_radar_freeboard = envi.radar_freeboard.values[~np.isnan(envi.lat)]
    envi_lat = envi.lat.values[~np.isnan(envi.lat)]
    envi_lon = envi.lon.values[~np.isnan(envi.lat)]
    envi_radar_freeboard_unc = envi.radar_freeboard_uncertainty.values[~np.isnan(envi.lat)]
    envi_track = track[~np.isnan(envi.lat)]
    
    data = {'lat':envi_lat,
            'lon':envi_lon,
            'time':envi_time,
            'radar_freeboard':envi_radar_freeboard,
            'radar_freeboard_uncertainty':envi_radar_freeboard_unc,
            'track':envi_track}
     
    envi = pd.DataFrame(data)
    
    coord = np.transpose(np.array([envi.lat.values.flatten(), envi.lon.values.flatten()]))
    coord[:,0], coord[:,1] = reproject(coord[:,0], coord[:,1])
    
    # # Land mask, remove land contaminated samples
    HB_coast = scipy.io.loadmat(docs+'Jack_onedrive/Existing Sea Ice Thickness Datasets/Landy et al 2017 ICESat GLAS13/HB_coast')['HB_coast']
    HB_coast[:,0], HB_coast[:,1] = reproject(HB_coast[:,0], HB_coast[:,1])
    
    MdlKDT = spatial.KDTree(HB_coast)
    land_idx = MdlKDT.query(coord, k=1, distance_upper_bound=10000)
    flat_land_idx = (land_idx[0]<10000)
    land_pts = np.zeros(len(envi))
    land_pts[flat_land_idx] = 1
    envi = envi[land_pts<1]
    coord = coord[land_pts<1,:]
    
    # intialise data of lists.
    data = {'freeboard':envi.radar_freeboard.values,
            'freeboard_uncertainty': envi.radar_freeboard_uncertainty.values,
            'track':envi.track.values,
            'lat':coord[:,0],
            'lon':coord[:,1],
            'year':year,
            'month':month,
            'day':pd.to_datetime(envi.time.values).day}
     
    # Create DataFrame
    envi = pd.DataFrame(data)
    #%%
    # Get snow density (snowmodel or warren)
    envi['snowdensity'] = snowdensity_snowmodel(envi)
    
    envi_w99 = envi.copy()
    envi_pmw = envi.copy()
    envi_sm = envi.copy()
    
    #%%
    # Snowdepth 
    envi_w99['snowdepth'] = snowdepth_warren(envi_w99, petty=False)
    envi_pmw['snowdepth'] = snowdepth_pmw(envi_pmw, petty=False)
    envi_sm['snowdepth'] = snowdepth_snowmodel(envi_sm, petty=False)
    
    #%%
    #Calculate sea ice thickness
    ow_density=1023.9; si_density_v1=915.1
    
    envi_w99['sit'],__,__ = sit_functions.sit_radar(envi_w99)
    envi_pmw['sit'],__,__ = sit_functions.sit_radar(envi_pmw)
    envi_sm['sit'],__,__ = sit_functions.sit_radar(envi_sm)
    
    envi_w99 = envi_w99[envi_w99.sit>0]
    envi_pmw = envi_pmw[envi_pmw.sit>0]
    envi_sm = envi_sm[envi_sm.sit>0]
    
    print('finished sea ice thickness')
    #%%
    #Uncertainty
    sit_uncertainty_w99 = sit_functions.uncertainty(envi_w99, icesat2=True)
    sit_uncertainty_pmw = sit_functions.uncertainty(envi_pmw, icesat2=True)
    sit_uncertainty_snowmodel = sit_functions.uncertainty(envi_sm, icesat2=True)
    
    #%%
    ## Create X km grid of ice thickness over study area
    reslat_plt, reslon_plt, reslat, reslon = sit_functions.create_grid(resolution=25)
    
    ## Grid Filters
    # Inverse linear mean filter (weighted by sea level uncertainty + distance)
    # per month
    search_radius = 50 #km
    sea_ice_grid_w99, sit_uncer_w99 = sit_functions.gridding(reslat, reslon, envi_w99, sit_uncertainty_w99, search_radius)
    sea_ice_grid_pmw, sit_uncer_pmw = sit_functions.gridding(reslat, reslon, envi_pmw, sit_uncertainty_pmw, search_radius)
    sea_ice_grid_snowmodel, sit_uncer_sm = sit_functions.gridding(reslat, reslon, envi_sm, sit_uncertainty_snowmodel, search_radius)
        
    print('finished grid')
    #%%
    #Save data
    # intialise data of lists.    
    data = {'lat':reslat_plt,
            'lon':reslon_plt,
            'sea_ice_thickness_w99':sea_ice_grid_w99,
            'sea_ice_thickness_pmw':sea_ice_grid_pmw,
            'sea_ice_thickness_snowmodel':sea_ice_grid_snowmodel,
            'sit_uncer_w99':sit_uncer_w99,
            'sit_uncer_pmw':sit_uncer_pmw,
            'sit_uncer_sm':sit_uncer_sm}
    
    # Create DataFrame + Save
    envi_sit = pd.DataFrame(data)
    envi_sit.to_csv('../../data/processed/envi/sit_all_'+str(year)+'_'+str("{:02d}".format(month))+'_smdens.csv',
                  index=False)
    
    print('finished '+str(year))

end = time.time()
print('Script took ',end - start,' s')
