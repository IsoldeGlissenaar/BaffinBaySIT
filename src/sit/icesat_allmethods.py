# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 11:34:11 2021

@author: IsoldeGlissenaar

Uses ICESat along track freeboard determined in icesat_freeboard.py
"""

import sys
sys.path.append('../functions')
from func_reproj import reproject, backproject
import sit_functions
from func_snowdepth import snowdepth_warren, snowdepth_pmw, snowdepth_snowmodel, effective_snowdepth, snowdensity_warren, snowdensity_snowmodel

import numpy as np
import pandas as pd
import scipy.interpolate
from scipy import spatial
import scipy.io
import time
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

start = time.time()

#Select years and month
years = np.arange(2003,2010,1)
month = 3

for year in years:
    
    #Import data from ice_thickness_from_ICESat script
    data = np.load('../../data/interim/ICESat_freeboard/data_cell2.npy')
    sea_level = np.load('../../data/interim/ICESat_freeboard/sea_level_v2_cell2.npy')
    sea_level_uncertainty = np.load('../../data/interim/ICESat_freeboard/sea_level_uncertainty_cell2.npy')
    freeboard = np.load('../../data/interim/ICESat_freeboard/freeboard_cell2.npy')
    dist_closest_lead = np.load('../../data/interim/ICESat_freeboard/dist_closest_lead_cell2.npy') 

    # Include only March of given year
    ids = (data[:,3]==month)*(data[:,2]==year) 
    sea_level = sea_level[ids]
    sea_level_uncertainty = sea_level_uncertainty[ids]
    dist_closest_lead = dist_closest_lead[ids]
    data = data[ids,:]
    freeboard = freeboard[ids]
    
    # # Land mask, remove land contaminated samples
    coord = data[:,7:9]
    HB_coast = scipy.io.loadmat(docs+'/Jack_onedrive/Existing Sea Ice Thickness Datasets/Landy et al 2017 ICESat GLAS13/HB_coast')['HB_coast']
    HB_coast[:,0], HB_coast[:,1] = reproject(HB_coast[:,0], HB_coast[:,1])
    
    MdlKDT = spatial.KDTree(HB_coast)
    land_idx = MdlKDT.query(coord, k=1, distance_upper_bound=10000)
    flat_land_idx = (land_idx[0]<10000)
    land_pts = np.zeros(len(freeboard))
    land_pts[flat_land_idx] = 1
    
    #Remove freeboard=0 values
    nonzero_fb = (freeboard!=0)&(land_pts<1)&(dist_closest_lead<=200000)
    sea_level = sea_level[nonzero_fb]
    sea_level_uncertainty = sea_level_uncertainty[nonzero_fb]
    data = data[nonzero_fb,:]
    freeboard = freeboard[nonzero_fb]
    coord = coord[nonzero_fb,:]
    dist_closest_lead = dist_closest_lead[nonzero_fb]
    
    day = data[:,4]
    
    # intialise data of lists.
    data2 = {'freeboard':freeboard,
             'sea_level': sea_level,
             'sea_level_uncertainty': sea_level_uncertainty,
             'geoid_height': data[:,13],
             'icesat_height': data[:,9],
             'lat':data[:,7],
             'lon':data[:,8],
             'year':year,
             'month':month,
             'day':day}
     
    # Create DataFrame
    icesat = pd.DataFrame(data2)
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
    
    icesat['track'] = track
    
    #%%
    #Get snow density (snowmodel or warren)
    icesat['snowdensity'] = snowdensity_warren(icesat)
    
    #fb uncertainty
    icesat['freeboard_uncertainty'] = sit_functions.freeboard_uncertainty(icesat, dist_closest_lead)
    
    print('finished snow density')
    #%%
    icesat_w99 = icesat.copy()
    icesat_w99_eff = icesat.copy()
    icesat_w99_petty = icesat.copy()
    
    icesat_pmw = icesat.copy()
    icesat_pmw_eff = icesat.copy()
    icesat_pmw_petty = icesat.copy()
    
    icesat_sm = icesat.copy()
    icesat_sm_eff = icesat.copy()
    icesat_sm_petty = icesat.copy()
    
    #%%
    # Snowdepth 
    icesat_w99['snowdepth'], icesat_w99_petty['snowdepth'] = snowdepth_warren(icesat, petty=True)
    icesat_pmw['snowdepth'], icesat_pmw_petty['snowdepth'] = snowdepth_pmw(icesat, petty=True) 
    icesat_sm['snowdepth'], icesat_sm_petty['snowdepth'] = snowdepth_snowmodel(icesat, petty=True)
    
    icesat_w99_eff['snowdepth'] = effective_snowdepth(icesat_w99, icesat_w99.snowdepth)
    icesat_pmw_eff['snowdepth'] = effective_snowdepth(icesat_pmw, icesat_pmw.snowdepth)
    icesat_sm_eff['snowdepth'] = effective_snowdepth(icesat_sm, icesat_sm.snowdepth)
    
    print('finished snow depth')
    #%%
    #Set snowdepth=fb when SIT<0
    ow_density=1023.9; si_density_v1=915.1
    minn = ow_density/(ow_density-icesat.snowdensity.values)*icesat.freeboard.values
    
    icesat_w99.snowdepth.values[icesat_w99.snowdepth.values>minn] = icesat.freeboard.values[icesat_w99.snowdepth.values>minn]
    icesat_w99_eff.snowdepth.values[icesat_w99_eff.snowdepth.values>minn] = icesat.freeboard.values[icesat_w99_eff.snowdepth.values>minn]
    icesat_w99_petty.snowdepth.values[icesat_w99_petty.snowdepth.values>minn] = icesat.freeboard.values[icesat_w99_petty.snowdepth.values>minn]
    
    icesat_pmw.snowdepth.values[icesat_pmw.snowdepth.values>minn] = icesat.freeboard.values[icesat_pmw.snowdepth.values>minn]
    icesat_pmw_eff.snowdepth.values[icesat_pmw_eff.snowdepth.values>minn] = icesat.freeboard.values[icesat_pmw_eff.snowdepth.values>minn]
    icesat_pmw_petty.snowdepth.values[icesat_pmw_petty.snowdepth.values>minn] = icesat.freeboard.values[icesat_pmw_petty.snowdepth.values>minn]
    
    icesat_sm.snowdepth.values[icesat_sm.snowdepth.values>minn] = icesat.freeboard.values[icesat_sm.snowdepth.values>minn]
    icesat_sm_eff.snowdepth.values[icesat_sm_eff.snowdepth.values>minn] = icesat.freeboard.values[icesat_sm_eff.snowdepth.values>minn]
    icesat_sm_petty.snowdepth.values[icesat_sm_petty.snowdepth.values>minn] = icesat.freeboard.values[icesat_sm_petty.snowdepth.values>minn]
    
    #%%
    #Calculate sea ice thickness
    
    icesat_w99['sit'],__,__ = sit_functions.sit_laser(icesat_w99)
    icesat_w99_eff['sit'],__,__ = sit_functions.sit_laser(icesat_w99_eff)
    icesat_w99_petty['sit'],__,__ = sit_functions.sit_laser(icesat_w99_petty)
    
    icesat_pmw['sit'],__,__ = sit_functions.sit_laser(icesat_pmw)
    icesat_pmw_eff['sit'],__,__ = sit_functions.sit_laser(icesat_pmw_eff)
    icesat_pmw_petty['sit'],__,__ = sit_functions.sit_laser(icesat_pmw_petty)
    
    icesat_sm['sit'],__,__ = sit_functions.sit_laser(icesat_sm)
    icesat_sm_eff['sit'],__,__ = sit_functions.sit_laser(icesat_sm_eff)
    icesat_sm_petty['sit'],__,__ = sit_functions.sit_laser(icesat_sm_petty)
    
    print('finished sea ice thickness')
    #%%
    ## Uncertainty ##
    sit_uncertainty_w99 = sit_functions.uncertainty(icesat_w99, icesat2=True)
    sit_uncertainty_w99_eff = sit_functions.uncertainty(icesat_w99_eff, icesat2=True)
    sit_uncertainty_w99_petty = sit_functions.uncertainty(icesat_w99_petty, icesat2=True)
    
    sit_uncertainty_pmw = sit_functions.uncertainty(icesat_pmw, icesat2=True)
    sit_uncertainty_pmw_eff = sit_functions.uncertainty(icesat_pmw_eff, icesat2=True)
    sit_uncertainty_pmw_petty = sit_functions.uncertainty(icesat_pmw_petty, icesat2=True)
    
    sit_uncertainty_sm = sit_functions.uncertainty(icesat_sm, icesat2=True)
    sit_uncertainty_sm_eff = sit_functions.uncertainty(icesat_sm_eff, icesat2=True)
    sit_uncertainty_sm_petty = sit_functions.uncertainty(icesat_sm_petty, icesat2=True)
    
    #%%
    ## Create X km grid of ice thickness over study area
    reslat_plt, reslon_plt, reslat, reslon = sit_functions.create_grid(resolution=25)
    
    ## Grid Filters
    # Inverse linear mean filter (weighted by sea level uncertainty + distance)
    # per month
    search_radius = 50 #km
    sea_ice_grid_w99, sit_uncer_w99 = sit_functions.gridding(reslat, reslon, icesat_w99, sit_uncertainty_w99, search_radius)
    sea_ice_grid_w99_eff, sit_uncer_w99_eff = sit_functions.gridding(reslat, reslon, icesat_w99_eff, sit_uncertainty_w99_eff, search_radius)
    sea_ice_grid_w99_petty, sit_uncer_w99_petty = sit_functions.gridding(reslat, reslon, icesat_w99_petty, sit_uncertainty_w99_petty, search_radius)
    
    sea_ice_grid_pmw, sit_uncer_pmw = sit_functions.gridding(reslat, reslon, icesat_pmw, sit_uncertainty_pmw, search_radius)
    sea_ice_grid_pmw_eff, sit_uncer_pmw_eff = sit_functions.gridding(reslat, reslon, icesat_pmw_eff, sit_uncertainty_pmw_eff, search_radius)
    sea_ice_grid_pmw_petty, sit_uncer_pmw_petty = sit_functions.gridding(reslat, reslon, icesat_pmw_petty, sit_uncertainty_pmw_petty, search_radius)
    
    sea_ice_grid_snowmodel, sit_uncer_sm = sit_functions.gridding(reslat, reslon, icesat_sm, sit_uncertainty_sm, search_radius)
    sea_ice_grid_snowmodel_eff, sit_uncer_sm_eff = sit_functions.gridding(reslat, reslon, icesat_sm_eff, sit_uncertainty_sm_eff, search_radius)
    sea_ice_grid_snowmodel_petty, sit_uncer_sm_petty = sit_functions.gridding(reslat, reslon, icesat_sm_petty, sit_uncertainty_sm_petty, search_radius)
    
    print('finished grid')
        
    #%%
    # #Save data
    # # intialise data of lists.
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
    
            'sit_uncer_w99':sit_uncer_w99,
            'sit_uncer_w99_eff':sit_uncer_w99_eff,
            'sit_uncer_w99_petty':sit_uncer_w99_petty,
            'sit_uncer_pmw':sit_uncer_pmw,
            'sit_uncer_pmw_eff':sit_uncer_pmw_eff,
            'sit_uncer_pmw_petty':sit_uncer_pmw_petty,
            'sit_uncer_snowmodel':sit_uncer_sm,
            'sit_uncer_snowmodel_eff':sit_uncer_sm_eff,
            'sit_uncer_snowmodel_petty':sit_uncer_sm_petty
            }
     
    # # Create DataFrame + save
    is2_sit = pd.DataFrame(data)
    is2_sit.to_csv('../../data/processed/is/sit_all_'+str(year)+'_'+str("{:02d}".format(month))+'_warrendens.csv',
                    index=False)
    
    print('finished '+str(year))

end = time.time()
print('Script took ',end - start,' s')
