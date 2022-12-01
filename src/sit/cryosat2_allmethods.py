# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 14:09:39 2021

@author: IsoldeGlissenaar

Uses CryoSat-2 alongtrack LARM freeboard according to Landy et al. (2020)
"""


import sys
sys.path.append('../functions')
from func_reproj import reproject
from func_snowdepth import snowdepth_warren, snowdepth_pmw, snowdepth_snowmodel, snowdensity_warren, snowdensity_snowmodel, snowdepth_nesosim
import sit_functions

import numpy as np
import pandas as pd
import scipy.io
from scipy import spatial
import time
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

start = time.time()

#Select years and month
years = np.arange(2011,2021,1)
month = 3

for year in years:
    cs2 = pd.read_csv(docs+'SatelliteData/CryoSat-2/dataframe/ubristol_trajectory_rfb_'+str(year)+'_'+str("{:02d}".format(month))+'.txt.csv')
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
            'sea_level_uncertainty': cs2.sea_level_uncertainty.values,
            'track':cs2.track.values,
            'sar_mode':cs2.sar_mode.values,
            'lat':coord[:,0],
            'lon':coord[:,1],
            'year':year,
            'month':month,
            'day':cs2.day.values}
     
    # Create DataFrame
    cs2 = pd.DataFrame(data)
    #%%
    #Get snow density (snowmodel or warren)
    cs2['snowdensity'] = snowdensity_snowmodel(cs2)
    
    cs2_w99 = cs2.copy()
    cs2_pmw = cs2.copy()
    cs2_sm = cs2.copy()
    cs2_nesosim = cs2.copy()
    
    #%%
    # Snowdepth 
    cs2_w99['snowdepth'] = snowdepth_warren(cs2_w99, petty=False)
    cs2_pmw['snowdepth'] = snowdepth_pmw(cs2_pmw, petty=False)
    cs2_sm['snowdepth'] = snowdepth_snowmodel(cs2_sm, petty=False)
    cs2_nesosim['snowdepth'] = snowdepth_nesosim(cs2_nesosim, petty=False)
    
    #%%
    #Calculate sea ice thickness
    ow_density=1023.9; si_density_v1=915.1
    
    cs2_w99['sit'],__,__ = sit_functions.sit_radar(cs2_w99)
    cs2_pmw['sit'],__,__ = sit_functions.sit_radar(cs2_pmw)
    cs2_sm['sit'],__,__ = sit_functions.sit_radar(cs2_sm)
    cs2_nesosim['sit'],__,__ = sit_functions.sit_radar(cs2_nesosim)
    
    cs2_w99 = cs2_w99[cs2_w99.sit>0]
    cs2_pmw = cs2_pmw[cs2_pmw.sit>0]
    cs2_sm = cs2_sm[cs2_sm.sit>0]
    cs2_nesosim = cs2_nesosim[cs2_nesosim.sit>0]
    
    print('finished sea ice thickness')
    #%%
    #Uncertainty
    sit_uncertainty_w99 = sit_functions.uncertainty(cs2_w99, speckle=True)
    sit_uncertainty_pmw = sit_functions.uncertainty(cs2_pmw, speckle=True)
    sit_uncertainty_snowmodel = sit_functions.uncertainty(cs2_sm, speckle=True)
    sit_uncertainty_nesosim = sit_functions.uncertainty(cs2_nesosim, speckle=True)
    
    #%%
    ## Create X km grid of ice thickness over study area
    reslat_plt, reslon_plt, reslat, reslon = sit_functions.create_grid(resolution=12)
    
    ## Grid Filters
    # Inverse linear mean filter (weighted by sea level uncertainty + distance)
    # per month
    search_radius = 24 #km
    sea_ice_grid_w99, sit_uncer_w99 = sit_functions.gridding(reslat, reslon, cs2_w99, sit_uncertainty_w99, search_radius, speckle=True)
    sea_ice_grid_pmw, sit_uncer_pmw = sit_functions.gridding(reslat, reslon, cs2_pmw, sit_uncertainty_pmw, search_radius, speckle=True)
    sea_ice_grid_snowmodel, sit_uncer_sm = sit_functions.gridding(reslat, reslon, cs2_sm, sit_uncertainty_snowmodel, search_radius, speckle=True)
    sea_ice_grid_nesosim, sit_uncer_nesosim = sit_functions.gridding(reslat, reslon, cs2_nesosim, sit_uncertainty_nesosim, search_radius, speckle=True)
        
    print('finished grid')
    
    #%%
    #Save data
    # intialise data of lists.
    data = {'lat':reslat_plt,
            'lon':reslon_plt,
            'sea_ice_thickness_w99':sea_ice_grid_w99,
            'sea_ice_thickness_pmw':sea_ice_grid_pmw,
            'sea_ice_thickness_snowmodel':sea_ice_grid_snowmodel,   
            'sea_ice_thickness_nesosim':sea_ice_grid_nesosim,
            'sit_uncer_w99':sit_uncer_w99,
            'sit_uncer_pmw':sit_uncer_pmw,
            'sit_uncer_sm':sit_uncer_sm,
            'sit_uncer_nesosim':sit_uncer_nesosim}
    
    # Create DataFrame + save
    cs2_sit = pd.DataFrame(data)
    cs2_sit.to_csv('../../data/processed/cs2/sit_all_'+str(year)+'_'+str("{:02d}".format(month))+'_smdens.csv',
                   index=False)
    
    print('finished '+str(year))

end = time.time()
print('Script took ',end - start,' s')
