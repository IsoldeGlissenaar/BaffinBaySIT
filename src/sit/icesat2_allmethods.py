# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 08:23:55 2021

@author: IsoldeGlissenaar

Uses ICESat-2 L3A ATL10 alongtrack freeboard product. Can be retrieved from
NSIDC: https://nsidc.org/data/icesat-2/data-sets
"""

import sys
sys.path.append('../functions')
from func_reproj import reproject
import sit_functions
from func_snowdepth import snowdepth_warren, snowdepth_pmw, snowdepth_snowmodel, effective_snowdepth, snowdensity_warren, snowdensity_snowmodel, snowdepth_nesosim

import numpy as np
import pandas as pd
import scipy.interpolate
import scipy.io
from datetime import datetime, timedelta
import time
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

start = time.time()

#Select month and years
month = 3
years = np.arange(2019,2021,1)

for year in years:
    is2 = pd.read_csv(docs+'/SatelliteData/ICESat-2/'+str(year)+'_'+str("{:02d}".format(month))+'/dataframe/ATL10_'+str(year)+str("{:02d}".format(month))+'_noland.csv')
    
    day = np.zeros(len(is2))
    for i in range(len(is2)):
        date2 = datetime(2018, 1, 1) + timedelta(seconds=is2.time.values[i])
        day[i] = date2.day
    
    coord = np.transpose(np.array([is2.lat.values, is2.lon.values]))
    coord[:,0], coord[:,1] = reproject(coord[:,0], coord[:,1])
    
    is2.freeboard_uncertainty.values[is2.freeboard_uncertainty.values>1e38]=np.nan
    
    # intialise data of lists.
    data = {'freeboard':is2.freeboard.values,
            'freeboard_uncertainty': is2.freeboard_uncertainty.values,
            'lat':coord[:,0],
            'lon':coord[:,1],
            'year':year,
            'month':month,
            'day':day}
     
    # Create DataFrame
    icesat2 = pd.DataFrame(data)
    #%%
    # Find track no.
    lat_roll = np.roll(is2.lat.values,-1)
    lon_roll = np.roll(is2.lon.values,-1)
    steps_lat = lat_roll-is2.lat.values
    a = (np.abs(steps_lat)>0.5)
    steps_lon = lon_roll-is2.lon.values
    b = (np.abs(steps_lon)>0.5)
    steps = (a|b)
    track = np.zeros(len(is2.lat.values))
    p=1
    for i in range(len(is2.lon.values)):
        if steps[i]==True:
            p=p+1
        track[i] = p
    
    icesat2['track'] = track
    
    #%%
    #Get snow density (snowmodel or warren)
    icesat2['snowdensity'] = snowdensity_warren(icesat2)
    
    icesat2_w99 = icesat2.copy()
    icesat2_w99_eff = icesat2.copy()
    icesat2_w99_petty = icesat2.copy()
    icesat2_pmw = icesat2.copy()
    icesat2_pmw_eff = icesat2.copy()
    icesat2_pmw_petty = icesat2.copy()
    icesat2_sm = icesat2.copy()
    icesat2_sm_eff = icesat2.copy()
    icesat2_sm_petty = icesat2.copy()
    
    icesat2_nesosim_petty = icesat2.copy()    
    
    #%%
    # Snowdepth 
    icesat2_w99['snowdepth'], icesat2_w99_petty['snowdepth'] = snowdepth_warren(icesat2_w99, petty=True)
    icesat2_pmw['snowdepth'], icesat2_pmw_petty['snowdepth'] = snowdepth_pmw(icesat2_pmw, petty=True) 
    icesat2_sm['snowdepth'], icesat2_sm_petty['snowdepth'] = snowdepth_snowmodel(icesat2_sm, petty=True)
    __, icesat2_nesosim_petty['snowdepth'] = snowdepth_nesosim(icesat2_nesosim_petty, petty=True)
         
    icesat2_w99_eff['snowdepth'] = effective_snowdepth(icesat2_w99_eff, icesat2_w99['snowdepth'])
    icesat2_pmw_eff['snowdepth'] = effective_snowdepth(icesat2_pmw_eff, icesat2_pmw['snowdepth'])
    icesat2_sm_eff['snowdepth'] = effective_snowdepth(icesat2_sm_eff, icesat2_sm['snowdepth'])
    
    
    #%%
    # Set snowdepth = fb where SIT<0
    ow_density=1023.9; si_density_v1=915.1
    minn = ow_density/(ow_density-icesat2.snowdensity.values)*icesat2.freeboard.values
    
    icesat2_w99.snowdepth.values[icesat2_w99.snowdepth.values>minn] = icesat2_w99.freeboard.values[icesat2_w99.snowdepth.values>minn]
    icesat2_pmw.snowdepth.values[icesat2_pmw.snowdepth.values>minn] = icesat2_pmw.freeboard.values[icesat2_pmw.snowdepth.values>minn]
    icesat2_sm.snowdepth.values[icesat2_sm.snowdepth.values>minn] = icesat2_sm.freeboard.values[icesat2_sm.snowdepth.values>minn]
    
    #%%
    #Calculate sea ice thickness
    
    icesat2_w99['sit'],__,__ = sit_functions.sit_laser(icesat2_w99)
    icesat2_w99_eff['sit'],__,__ = sit_functions.sit_laser(icesat2_w99_eff)
    icesat2_w99_petty['sit'],__,__ = sit_functions.sit_laser(icesat2_w99_petty)
    
    icesat2_pmw['sit'],__,__ = sit_functions.sit_laser(icesat2_pmw)
    icesat2_pmw_eff['sit'],__,__ = sit_functions.sit_laser(icesat2_pmw_eff)
    icesat2_pmw_petty['sit'],__,__ = sit_functions.sit_laser(icesat2_pmw_petty)
    
    icesat2_sm['sit'],__,__ = sit_functions.sit_laser(icesat2_sm)
    icesat2_sm_eff['sit'],__,__ = sit_functions.sit_laser(icesat2_sm_eff)
    icesat2_sm_petty['sit'],__,__ = sit_functions.sit_laser(icesat2_sm_petty)
    
    icesat2_nesosim_petty['sit'],__,__ = sit_functions.sit_laser(icesat2_nesosim_petty)

        
    print('finished sea ice thickness')
    #%%
    #Uncertainty
    sit_uncertainty_w99 = sit_functions.uncertainty(icesat2_w99, icesat2=True)
    sit_uncertainty_w99_eff = sit_functions.uncertainty(icesat2_w99_eff, icesat2=True)
    sit_uncertainty_w99_petty = sit_functions.uncertainty(icesat2_w99_petty, icesat2=True)
    
    sit_uncertainty_pmw = sit_functions.uncertainty(icesat2_pmw, icesat2=True)
    sit_uncertainty_pmw_eff = sit_functions.uncertainty(icesat2_pmw_eff, icesat2=True)
    sit_uncertainty_pmw_petty = sit_functions.uncertainty(icesat2_pmw_petty, icesat2=True)
    
    sit_uncertainty_sm = sit_functions.uncertainty(icesat2_sm, icesat2=True)
    sit_uncertainty_sm_eff = sit_functions.uncertainty(icesat2_sm_eff, icesat2=True)
    sit_uncertainty_sm_petty = sit_functions.uncertainty(icesat2_sm_petty, icesat2=True)
    
    sit_uncertainty_nesosim_petty = sit_functions.uncertainty(icesat2_nesosim_petty, icesat2=True)
    
    #%%
    ## Create X km grid of ice thickness over study area
    reslat_plt, reslon_plt, reslat, reslon = sit_functions.create_grid(resolution=12)
    
    ## Grid Filters
    # Inverse linear mean filter (weighted by sea level uncertainty + distance)
    # per month
    search_radius = 24 #km
    sea_ice_grid_w99, sit_uncer_w99 = sit_functions.gridding(reslat, reslon, icesat2_w99, sit_uncertainty_w99, search_radius)
    sea_ice_grid_w99_eff, sit_uncer_w99_eff = sit_functions.gridding(reslat, reslon, icesat2_w99_eff, sit_uncertainty_w99_eff, search_radius)
    sea_ice_grid_w99_petty, sit_uncer_w99_petty = sit_functions.gridding(reslat, reslon, icesat2_w99_petty, sit_uncertainty_w99_petty, search_radius)
    
    sea_ice_grid_pmw, sit_uncer_pmw = sit_functions.gridding(reslat, reslon, icesat2_pmw, sit_uncertainty_pmw, search_radius)
    sea_ice_grid_pmw_eff, sit_uncer_pmw_eff = sit_functions.gridding(reslat, reslon, icesat2_pmw_eff, sit_uncertainty_pmw_eff, search_radius)
    sea_ice_grid_pmw_petty, sit_uncer_pmw_petty = sit_functions.gridding(reslat, reslon, icesat2_pmw_petty, sit_uncertainty_pmw_petty, search_radius)
    
    sea_ice_grid_snowmodel, sit_uncer_snowmodel = sit_functions.gridding(reslat, reslon, icesat2_sm, sit_uncertainty_sm, search_radius)
    sea_ice_grid_snowmodel_eff, sit_uncer_snowmodel_eff = sit_functions.gridding(reslat, reslon, icesat2_sm_eff, sit_uncertainty_sm_eff, search_radius)
    sea_ice_grid_snowmodel_petty, sit_uncer_snowmodel_petty = sit_functions.gridding(reslat, reslon, icesat2_sm_petty, sit_uncertainty_sm_petty, search_radius)
    
    sea_ice_grid_nesosim_petty, sit_uncer_nesosim_petty = sit_functions.gridding(reslat, reslon, icesat2_nesosim_petty, sit_uncertainty_nesosim_petty, search_radius)
            
    print('finished grid')
    
    #%%
    #Save data
    # intialise data of lists.
    data = {'lat':reslat_plt,
            'lon':reslon_plt,
            'sea_ice_thickness_w99':sea_ice_grid_w99,
            'sea_ice_thickness_w99_eff':sea_ice_grid_w99_eff,
            'sea_ice_thickness_w99_petty':sea_ice_grid_w99_petty,
            'sea_ice_thickness_pmw':sea_ice_grid_pmw,
            'sea_ice_thickness_pmw_eff':sea_ice_grid_pmw_eff,
            'sea_ice_thickness_pmw_petty':sea_ice_grid_pmw_petty,
            'sea_ice_thickness_snowmodel':sea_ice_grid_snowmodel,
            'sea_ice_thickness_snowmodel_eff':sea_ice_grid_snowmodel_eff,
            'sea_ice_thickness_snowmodel_petty':sea_ice_grid_snowmodel_petty,
            'sea_ice_thickness_nesosim_petty':sea_ice_grid_nesosim_petty,            
            
            'sit_uncer_w99':sit_uncer_w99,
            'sit_uncer_w99_eff':sit_uncer_w99_eff,
            'sit_uncer_w99_petty':sit_uncer_w99_petty,
            'sit_uncer_pmw':sit_uncer_pmw,
            'sit_uncer_pmw_eff':sit_uncer_pmw_eff,
            'sit_uncer_pmw_petty':sit_uncer_pmw_petty,
            'sit_uncer_snowmodel':sit_uncer_snowmodel,
            'sit_uncer_snowmodel_eff':sit_uncer_snowmodel_eff,
            'sit_uncer_snowmodel_petty':sit_uncer_snowmodel_petty,
            'sit_uncer_nesosim_petty':sit_uncer_nesosim_petty
            }
     
    # Create DataFrame + save
    is2_sit = pd.DataFrame(data)
    is2_sit.to_csv('../../data/processed/is2/sit_all_'+str(year)+'_'+str("{:02d}".format(month))+'_warrendens.csv',
                   index=False)
    
    print('finished '+str(year))

end = time.time()
print('Script took ',end - start,' s')
