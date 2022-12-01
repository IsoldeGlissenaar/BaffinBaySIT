# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 15:19:19 2021

@author: IsoldeGlissenaar

Uses CryoSat-2 ESA CCI freeboard observations. Can be retrieved from:
https://catalogue.ceda.ac.uk/uuid/5b6033bfb7f241e89132a83fdc3d5364
"""

import sys
sys.path.append('../functions')
from func_reproj import reproject, backproject
from func_snowdepth import snowdepth_warren, snowdepth_pmw, snowdepth_snowmodel, snowdensity_warren, snowdensity_snowmodel
import sit_functions

import numpy as np
import pandas as pd
import scipy.io
from scipy import spatial
import time
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

start = time.time()

#Select month and years
month = 3
years = np.arange(2011,2018,1)

for year in years:
    cs2 = pd.read_csv(docs+'SatelliteData/CryoSat-2_CCI/dataframe/cs2_cci_'+str(year)+str("{:02d}".format(month))+'.csv')
    coord = np.transpose(np.array([cs2.lat.values, cs2.lon.values]))
    coord[:,0], coord[:,1] = reproject(coord[:,0], coord[:,1])
    
    # # Land mask, remove land contaminated samples
    HB_coast = scipy.io.loadmat(docs+'Jack_onedrive/Existing Sea Ice Thickness Datasets/Landy et al 2017 ICESat GLAS13/HB_coast')['HB_coast']
    HB_coast[:,0], HB_coast[:,1] = reproject(HB_coast[:,0], HB_coast[:,1])
    
    MdlKDT = spatial.KDTree(HB_coast)
    land_idx = MdlKDT.query(coord, k=1, distance_upper_bound=10000)
    flat_land_idx = (land_idx[0]<10000)
    land_pts = np.zeros(len(cs2))
    land_pts[flat_land_idx] = 1
    cs2 = cs2[land_pts<1]
    coord = coord[land_pts<1,:]
    
    # intialise data of lists.
    data = {'freeboard':cs2.freeboard.values,
            'freeboard_uncertainty': cs2.freeboard_uncertainty.values,
            'lat':coord[:,0],
            'lon':coord[:,1],
            'year':year,
            'month':month,
            'day':cs2.day.values}
     
    # Create DataFrame
    cs2 = pd.DataFrame(data)
    #%%
    #Find track no.
    coord_wgs = np.zeros(coord.shape)
    coord_wgs[:,0], coord_wgs[:,1] = backproject(coord[:,0],coord[:,1])
    
    # Find track
    lat_roll = np.roll(coord_wgs[:,0],-1)
    lon_roll = np.roll(coord_wgs[:,1],-1)
    steps_lat = lat_roll-coord_wgs[:,0]
    a = (np.abs(steps_lat)>0.5)
    steps_lon = lon_roll-coord_wgs[:,1]
    b = (np.abs(steps_lon)>0.5)
    steps = (a|b)
    track = np.zeros(len(coord_wgs[:,0]))
    p=1
    for i in range(len(coord_wgs[:,0])):
        if steps[i]==True:
            p=p+1
        track[i] = p
    
    cs2['track'] = track
    
    #%%
    #Get snow density (snowmodel or warren)
    cs2['snowdensity'] = snowdensity_snowmodel(cs2)
    
    cs2_w99 = cs2.copy()
    cs2_pmw = cs2.copy()
    cs2_sm = cs2.copy()

    #%%
    # Snowdepth 
    cs2_w99['snowdepth'] = snowdepth_warren(cs2_w99, petty=False)
    cs2_pmw['snowdepth'] = snowdepth_pmw(cs2_pmw, petty=False)
    cs2_sm['snowdepth'] = snowdepth_snowmodel(cs2_sm, petty=False)
    
    #%%
    #Calculate sea ice thickness
    ow_density=1023.9; si_density_v1=915.1
    
    cs2_w99['sit'],__,__ = sit_functions.sit_radar(cs2_w99)
    cs2_pmw['sit'],__,__ = sit_functions.sit_radar(cs2_pmw)
    cs2_sm['sit'],__,__ = sit_functions.sit_radar(cs2_sm)
    
    cs2_w99.sit[cs2_w99.sit.values<=0] = np.nan
    cs2_pmw.sit[cs2_pmw.sit.values<=0] = np.nan
    cs2_sm.sit[cs2_sm.sit.values<=0] = np.nan
    # cs2_w99 = cs2_w99.where(cs2_w99.sit>0)
    # cs2_pmw = cs2_pmw.where(cs2_pmw.sit>0)
    # cs2_sm = cs2_sm.where(cs2_sm.sit>0)
    
    print('finished sea ice thickness')
    #%%
    #Uncertainty
    sit_uncertainty_w99 = sit_functions.uncertainty(cs2_w99, icesat2=True)
    sit_uncertainty_pmw = sit_functions.uncertainty(cs2_pmw, icesat2=True)
    sit_uncertainty_snowmodel = sit_functions.uncertainty(cs2_sm, icesat2=True)
    
    #%%
    ## Create X km grid of ice thickness over study area
    reslat_plt, reslon_plt, reslat, reslon = sit_functions.create_grid(resolution=12)
    
    ## Grid Filters
    # Inverse linear mean filter (weighted by sea level uncertainty + distance)
    # per month
    search_radius = 24 #km
    sea_ice_grid_w99, sit_uncer_w99 = sit_functions.gridding(reslat, reslon, cs2_w99, sit_uncertainty_w99, search_radius)
    sea_ice_grid_pmw, sit_uncer_pmw = sit_functions.gridding(reslat, reslon, cs2_pmw, sit_uncertainty_pmw, search_radius)
    sea_ice_grid_snowmodel, sit_uncer_sm = sit_functions.gridding(reslat, reslon, cs2_sm, sit_uncertainty_snowmodel, search_radius)
    
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
    cs2_sit = pd.DataFrame(data)
    cs2_sit.to_csv('../../data/processed/cs2_cci/sit_all_'+str(year)+'_'+str("{:02d}".format(month))+'_smdens.csv',
                   index=False)
    
    print('finished '+str(year))

end = time.time()
print('Script took ',end - start,' s')
