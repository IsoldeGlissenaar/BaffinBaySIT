# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 11:41:50 2021

@author: IsoldeGlissenaar

Converts and combines SIT gridded data from CSV to NetCDF 
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir
import xarray as xr
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

#%%
'''ICESat Warren snow dens'''

direc = '../../data/processed/is/'
is_warren = np.zeros((3045,7)); is_warren[:,:] = np.nan
is_warren_petty = np.zeros((3045,7)); is_warren_petty[:,:] = np.nan
is_warren_eff = np.zeros((3045,7)); is_warren_eff[:,:] = np.nan
is_sit_pmw = np.zeros((3045,7)); is_sit_pmw[:,:] = np.nan
is_sit_pmw_effsd = np.zeros((3045,7)); is_sit_pmw_effsd[:,:] = np.nan
is_sit = np.zeros((3045,7)); is_sit[:,:] = np.nan
is_sit_stroeve = np.zeros((3045,7)); is_sit_stroeve[:,:] = np.nan
is_sit_stroeve_effsd = np.zeros((3045,7)); is_sit_stroeve_effsd[:,:] = np.nan
is_sit_stroeve_petty = np.zeros((3045,7)); is_sit_stroeve_petty[:,:] = np.nan

is_warren_u = np.zeros((3045,7)); is_warren_u[:,:] = np.nan
is_warren_petty_u = np.zeros((3045,7)); is_warren_petty_u[:,:] = np.nan
is_warren_eff_u = np.zeros((3045,7)); is_warren_eff_u[:,:] = np.nan
is_sit_pmw_u = np.zeros((3045,7)); is_sit_pmw_u[:,:] = np.nan
is_sit_pmw_effsd_u = np.zeros((3045,7)); is_sit_pmw_effsd_u[:,:] = np.nan
is_sit_u = np.zeros((3045,7)); is_sit_u[:,:] = np.nan
is_sit_stroeve_u = np.zeros((3045,7)); is_sit_stroeve_u[:,:] = np.nan
is_sit_stroeve_effsd_u = np.zeros((3045,7)); is_sit_stroeve_effsd_u[:,:] = np.nan
is_sit_stroeve_petty_u = np.zeros((3045,7)); is_sit_stroeve_petty_u[:,:] = np.nan

is_freeboard = np.zeros((3045,7))

i=0
for filename in listdir(direc):
    if filename.endswith('03_warrendens_inclfb.csv'):
        is_temp = pd.read_csv(direc+filename)
        is_warren[:,i] = is_temp.sea_ice_thickness_w99.values
        is_warren_petty[:,i] = is_temp.sea_ice_thickness_w99_petty.values
        is_warren_eff[:,i] = is_temp.sea_ice_thickness_w99_eff.values
        is_sit_pmw[:,i] = is_temp.sea_ice_thickness_pmw.values
        is_sit_pmw_effsd[:,i] = is_temp.sea_ice_thickness_pmw_eff.values
        is_sit[:,i] = is_temp.sea_ice_thickness_pmw_petty.values
        is_sit_stroeve[:,i] = is_temp.sea_ice_thickness_snowmodel.values
        is_sit_stroeve_effsd[:,i] = is_temp.sea_ice_thickness_snowmodel_eff.values
        is_sit_stroeve_petty[:,i] = is_temp.sea_ice_thickness_snowmodel_petty.values
        
        is_warren_u[:,i] = is_temp.sit_uncer_w99.values
        is_warren_petty_u[:,i] = is_temp.sit_uncer_w99_petty.values
        is_warren_eff_u[:,i] = is_temp.sit_uncer_w99_eff.values
        is_sit_pmw_u[:,i] = is_temp.sit_uncer_pmw.values
        is_sit_pmw_effsd_u[:,i] = is_temp.sit_uncer_pmw_eff.values
        is_sit_u[:,i] = is_temp.sit_uncer_pmw_petty.values
        is_sit_stroeve_u[:,i] = is_temp.sit_uncer_snowmodel.values
        is_sit_stroeve_effsd_u[:,i] = is_temp.sit_uncer_snowmodel_eff.values
        is_sit_stroeve_petty_u[:,i] = is_temp.sit_uncer_snowmodel_petty.values
        
        is_freeboard[:,i] = is_temp.freeboard.values
        i=i+1
        
lat25 = is_temp.lat.values
lon25 = is_temp.lon.values      

is_years = np.arange(2003,2010,1).astype(int)
n = np.arange(0,3045,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], is_warren),
        sit_w99_sig     =(["n","year"], is_warren_eff),
        sit_w99_piece   =(["n","year"], is_warren_petty),
        sit_pmw         =(["n","year"], is_sit_pmw),
        sit_pmw_sig     =(["n","year"], is_sit_pmw_effsd),
        sit_pmw_piece   =(["n","year"], is_sit),
        sit_sm          =(["n","year"], is_sit_stroeve),
        sit_sm_sig      =(["n","year"], is_sit_stroeve_effsd),
        sit_sm_piece    =(["n","year"], is_sit_stroeve_petty),
        sit_w99_uncertainty       =(["n","year"], is_warren_u),
        sit_w99_sig_uncertainty   =(["n","year"], is_warren_eff_u),
        sit_w99_piece_uncertainty =(["n","year"], is_warren_petty_u),
        sit_pmw_uncertainty       =(["n","year"], is_sit_pmw_u),
        sit_pmw_sig_uncertainty   =(["n","year"], is_sit_pmw_effsd_u),
        sit_pmw_piece_uncertainty =(["n","year"], is_sit_u),
        sit_sm_uncertainty        =(["n","year"], is_sit_stroeve_u),
        sit_sm_sig_uncertainty    =(["n","year"], is_sit_stroeve_effsd_u),
        sit_sm_piece_uncertainty  =(["n","year"], is_sit_stroeve_petty_u),
        freeboard       =(["n","year"], is_freeboard)
        ),
    coords=dict(
        lon=(["n"],lon25),
        lat=(["n"],lat25),
        year = is_years,
        ),
    attrs=dict(description="Gridded ICESat sea ice thickness [m] (Warren snowdensity)"))
    
data.to_netcdf(direc+'SIT_ICESat_March_2003_2009_Warrensnowdens_inclfb.nc')

#%%
"""ICESat - SnowModel snow density""" 

direc = '../../data/processed/is/'
is_warren = np.zeros((3045,7)); is_warren[:,:] = np.nan
is_warren_petty = np.zeros((3045,7)); is_warren_petty[:,:] = np.nan
is_warren_eff = np.zeros((3045,7)); is_warren_eff[:,:] = np.nan
is_sit_pmw = np.zeros((3045,7)); is_sit_pmw[:,:] = np.nan
is_sit_pmw_effsd = np.zeros((3045,7)); is_sit_pmw_effsd[:,:] = np.nan
is_sit = np.zeros((3045,7)); is_sit[:,:] = np.nan
is_sit_stroeve = np.zeros((3045,7)); is_sit_stroeve[:,:] = np.nan
is_sit_stroeve_effsd = np.zeros((3045,7)); is_sit_stroeve_effsd[:,:] = np.nan
is_sit_stroeve_petty = np.zeros((3045,7)); is_sit_stroeve_petty[:,:] = np.nan

is_warren_u = np.zeros((3045,7)); is_warren_u[:,:] = np.nan
is_warren_petty_u = np.zeros((3045,7)); is_warren_petty_u[:,:] = np.nan
is_warren_eff_u = np.zeros((3045,7)); is_warren_eff_u[:,:] = np.nan
is_sit_pmw_u = np.zeros((3045,7)); is_sit_pmw_u[:,:] = np.nan
is_sit_pmw_effsd_u = np.zeros((3045,7)); is_sit_pmw_effsd_u[:,:] = np.nan
is_sit_u = np.zeros((3045,7)); is_sit_u[:,:] = np.nan
is_sit_stroeve_u = np.zeros((3045,7)); is_sit_stroeve_u[:,:] = np.nan
is_sit_stroeve_effsd_u = np.zeros((3045,7)); is_sit_stroeve_effsd_u[:,:] = np.nan
is_sit_stroeve_petty_u = np.zeros((3045,7)); is_sit_stroeve_petty_u[:,:] = np.nan

is_freeboard = np.zeros((3045,7))

i=0
for filename in listdir(direc):
    if filename.endswith('03_smdens_inclfb.csv'):
        is_temp = pd.read_csv(direc+filename)
        is_warren[:,i] = is_temp.sea_ice_thickness_w99.values
        is_warren_petty[:,i] = is_temp.sea_ice_thickness_w99_petty.values
        is_warren_eff[:,i] = is_temp.sea_ice_thickness_w99_eff.values
        is_sit_pmw[:,i] = is_temp.sea_ice_thickness_pmw.values
        is_sit_pmw_effsd[:,i] = is_temp.sea_ice_thickness_pmw_eff.values
        is_sit[:,i] = is_temp.sea_ice_thickness_pmw_petty.values
        is_sit_stroeve[:,i] = is_temp.sea_ice_thickness_snowmodel.values
        is_sit_stroeve_effsd[:,i] = is_temp.sea_ice_thickness_snowmodel_eff.values
        is_sit_stroeve_petty[:,i] = is_temp.sea_ice_thickness_snowmodel_petty.values
        
        is_warren_u[:,i] = is_temp.sit_uncer_w99.values
        is_warren_petty_u[:,i] = is_temp.sit_uncer_w99_petty.values
        is_warren_eff_u[:,i] = is_temp.sit_uncer_w99_eff.values
        is_sit_pmw_u[:,i] = is_temp.sit_uncer_pmw.values
        is_sit_pmw_effsd_u[:,i] = is_temp.sit_uncer_pmw_eff.values
        is_sit_u[:,i] = is_temp.sit_uncer_pmw_petty.values
        is_sit_stroeve_u[:,i] = is_temp.sit_uncer_snowmodel.values
        is_sit_stroeve_effsd_u[:,i] = is_temp.sit_uncer_snowmodel_eff.values
        is_sit_stroeve_petty_u[:,i] = is_temp.sit_uncer_snowmodel_petty.values
        
        is_freeboard[:,i] = is_temp.freeboard.values
        i=i+1
        
lat25 = is_temp.lat.values
lon25 = is_temp.lon.values      

is_years = np.arange(2003,2010,1).astype(int)
n = np.arange(0,3045,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], is_warren),
        sit_w99_sig     =(["n","year"], is_warren_eff),
        sit_w99_piece   =(["n","year"], is_warren_petty),
        sit_pmw         =(["n","year"], is_sit_pmw),
        sit_pmw_sig     =(["n","year"], is_sit_pmw_effsd),
        sit_pmw_piece   =(["n","year"], is_sit),
        sit_sm          =(["n","year"], is_sit_stroeve),
        sit_sm_sig      =(["n","year"], is_sit_stroeve_effsd),
        sit_sm_piece    =(["n","year"], is_sit_stroeve_petty),
        sit_w99_uncertainty       =(["n","year"], is_warren_u),
        sit_w99_sig_uncertainty   =(["n","year"], is_warren_eff_u),
        sit_w99_piece_uncertainty =(["n","year"], is_warren_petty_u),
        sit_pmw_uncertainty       =(["n","year"], is_sit_pmw_u),
        sit_pmw_sig_uncertainty   =(["n","year"], is_sit_pmw_effsd_u),
        sit_pmw_piece_uncertainty =(["n","year"], is_sit_u),
        sit_sm_uncertainty        =(["n","year"], is_sit_stroeve_u),
        sit_sm_sig_uncertainty    =(["n","year"], is_sit_stroeve_effsd_u),
        sit_sm_piece_uncertainty  =(["n","year"], is_sit_stroeve_petty_u),
        freeboard       =(["n","year"], is_freeboard)
        ),
    coords=dict(
        lon=(["n"],lon25),
        lat=(["n"],lat25),
        year = is_years,
        ),
    attrs=dict(description="Gridded ICESat sea ice thickness [m] (SnowModel snowdensity)"))
    
data.to_netcdf(direc+'SIT_ICESat_March_2003_2009_SnowModelsnowdens_inclfb.nc')

#%%
"""CryoSat-2 Warren snowdensity"""

direc = '../../data/processed/cs2/'
cs2_sit = np.zeros((13357,10)); cs2_sit[:,:] = np.nan
cs2_warren = np.zeros((13357,10)); cs2_warren[:,:] = np.nan
cs2_stroeve = np.zeros((13357,10)); cs2_stroeve[:,:] = np.nan
cs2_nesosim = np.zeros((13357,10)); cs2_nesosim[:,:] = np.nan

cs2_sit_u = np.zeros((13357,10)); cs2_sit_u[:,:] = np.nan
cs2_warren_u = np.zeros((13357,10)); cs2_warren_u[:,:] = np.nan
cs2_stroeve_u = np.zeros((13357,10)); cs2_stroeve_u[:,:] = np.nan
cs2_nesosim_u = np.zeros((13357,10)); cs2_nesosim_u[:,:] = np.nan
cs2_years = np.zeros((10)).astype(int)
i=0
for filename in listdir(direc):
    if filename.endswith('03_warrendens.csv'):
        year = filename[8:12] #[14:18]
        f = pd.read_csv(direc+filename)
        cs2_sit[:,i] = f.sea_ice_thickness_pmw.values
        cs2_warren[:,i] = f.sea_ice_thickness_w99.values
        cs2_stroeve[:,i] = f.sea_ice_thickness_snowmodel.values
        cs2_nesosim[:,i] = f.sea_ice_thickness_nesosim.values
        
        cs2_sit_u[:,i] = f.sit_uncer_pmw.values
        cs2_warren_u[:,i] = f.sit_uncer_w99.values
        cs2_stroeve_u[:,i] = f.sit_uncer_sm.values
        cs2_nesosim_u[:,i] = f.sit_uncer_nesosim.values
        cs2_years[i] = year
        i=i+1
lat12 = f.lat.values
lon12 = f.lon.values

n = np.arange(0,13357,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], cs2_warren),
        sit_pmw         =(["n","year"], cs2_sit),
        sit_sm          =(["n","year"], cs2_stroeve),
        sit_w99_uncertainty       =(["n","year"], cs2_warren_u),
        sit_pmw_uncertainty       =(["n","year"], cs2_sit_u),
        sit_sm_uncertainty        =(["n","year"], cs2_stroeve_u)
        ),
    coords=dict(
        lon=(["n"],lon12),
        lat=(["n"],lat12),
        year = cs2_years,
        ),
    attrs=dict(description="Gridded CryoSat-2 LARM sea ice thickness [m] (Warren snowdensity"))
    
data.to_netcdf(direc+'SIT_CryoSat2_March_2011_2020_Warrensnowdens.nc')

#%%
"""CryoSat-2 SnowModel snowdensity"""

direc = '../../data/processed/cs2/'
cs2_sit = np.zeros((13357,10)); cs2_sit[:,:] = np.nan
cs2_warren = np.zeros((13357,10)); cs2_warren[:,:] = np.nan
cs2_stroeve = np.zeros((13357,10)); cs2_stroeve[:,:] = np.nan

cs2_sit_u = np.zeros((13357,10)); cs2_sit_u[:,:] = np.nan
cs2_warren_u = np.zeros((13357,10)); cs2_warren_u[:,:] = np.nan
cs2_stroeve_u = np.zeros((13357,10)); cs2_stroeve_u[:,:] = np.nan
cs2_years = np.zeros((10)).astype(int)
i=0
for filename in listdir(direc):
    if filename.endswith('03_smdens.csv'):
        year = filename[8:12] #[14:18]
        f = pd.read_csv(direc+filename)
        cs2_sit[:,i] = f.sea_ice_thickness_pmw.values
        cs2_warren[:,i] = f.sea_ice_thickness_w99.values
        cs2_stroeve[:,i] = f.sea_ice_thickness_snowmodel.values
        
        cs2_sit_u[:,i] = f.sit_uncer_pmw.values
        cs2_warren_u[:,i] = f.sit_uncer_w99.values
        cs2_stroeve_u[:,i] = f.sit_uncer_sm.values
        cs2_years[i] = year
        i=i+1
lat12 = f.lat.values
lon12 = f.lon.values

n = np.arange(0,13357,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], cs2_warren),
        sit_pmw         =(["n","year"], cs2_sit),
        sit_sm          =(["n","year"], cs2_stroeve),
        sit_w99_uncertainty       =(["n","year"], cs2_warren_u),
        sit_pmw_uncertainty       =(["n","year"], cs2_sit_u),
        sit_sm_uncertainty        =(["n","year"], cs2_stroeve_u)  
        ),
    coords=dict(
        lon=(["n"],lon12),
        lat=(["n"],lat12),
        year = cs2_years,
        ),
    attrs=dict(description="Gridded CryoSat-2 LARM sea ice thickness [m] (SnowModel snowdensity"))
    
data.to_netcdf(direc+'SIT_CryoSat2_March_2011_2020_SnowModelsnowdens.nc')
#%%
"""CryoSat-2 CCI Warren snowdensity"""

direc = '../../data/processed/cs2_cci/'
cs2_sit = np.zeros((13357,10)); cs2_sit[:,:] = np.nan
cs2_warren = np.zeros((13357,10)); cs2_warren[:,:] = np.nan
cs2_stroeve = np.zeros((13357,10)); cs2_stroeve[:,:] = np.nan

cs2_sit_u = np.zeros((13357,10)); cs2_sit_u[:,:] = np.nan
cs2_warren_u = np.zeros((13357,10)); cs2_warren_u[:,:] = np.nan
cs2_stroeve_u = np.zeros((13357,10)); cs2_stroeve_u[:,:] = np.nan
cs2_years = np.zeros((10)).astype(int)
i=0
for filename in listdir(direc):
    if filename.endswith('03_warrendens.csv'):
        year = filename[8:12] #[14:18]
        f = pd.read_csv(direc+filename)
        cs2_sit[:,i] = f.sea_ice_thickness_pmw.values
        cs2_warren[:,i] = f.sea_ice_thickness_w99.values
        cs2_stroeve[:,i] = f.sea_ice_thickness_snowmodel.values
        
        cs2_sit_u[:,i] = f.sit_uncer_pmw.values
        cs2_warren_u[:,i] = f.sit_uncer_w99.values
        cs2_stroeve_u[:,i] = f.sit_uncer_sm.values
        cs2_years[i] = year
        i=i+1
lat12 = f.lat.values
lon12 = f.lon.values

n = np.arange(0,13357,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], cs2_warren),
        sit_pmw         =(["n","year"], cs2_sit),
        sit_sm          =(["n","year"], cs2_stroeve),
        sit_w99_uncertainty       =(["n","year"], cs2_warren_u),
        sit_pmw_uncertainty       =(["n","year"], cs2_sit_u),
        sit_sm_uncertainty        =(["n","year"], cs2_stroeve_u)  
        ),
    coords=dict(
        lon=(["n"],lon12),
        lat=(["n"],lat12),
        year = cs2_years,
        ),
    attrs=dict(description="Gridded CryoSat-2 CCI sea ice thickness [m] (Warren snowdensity"))
    
data.to_netcdf(direc+'SIT_CryoSat2_March_2011_2017_Warrensnowdens.nc')

#%%
"""CryoSat-2 CCI SnowModel snowdensity"""

direc = '../../data/processed/cs2_cci/'
cs2_sit = np.zeros((13357,10)); cs2_sit[:,:] = np.nan
cs2_warren = np.zeros((13357,10)); cs2_warren[:,:] = np.nan
cs2_stroeve = np.zeros((13357,10)); cs2_stroeve[:,:] = np.nan

cs2_sit_u = np.zeros((13357,10)); cs2_sit_u[:,:] = np.nan
cs2_warren_u = np.zeros((13357,10)); cs2_warren_u[:,:] = np.nan
cs2_stroeve_u = np.zeros((13357,10)); cs2_stroeve_u[:,:] = np.nan
cs2_years = np.zeros((10)).astype(int)
i=0
for filename in listdir(direc):
    if filename.endswith('03_smdens.csv'):
        year = filename[8:12] #[14:18]
        f = pd.read_csv(direc+filename)
        cs2_sit[:,i] = f.sea_ice_thickness_pmw.values
        cs2_warren[:,i] = f.sea_ice_thickness_w99.values
        cs2_stroeve[:,i] = f.sea_ice_thickness_snowmodel.values
        
        cs2_sit_u[:,i] = f.sit_uncer_pmw.values
        cs2_warren_u[:,i] = f.sit_uncer_w99.values
        cs2_stroeve_u[:,i] = f.sit_uncer_sm.values
        cs2_years[i] = year
        i=i+1
lat12 = f.lat.values
lon12 = f.lon.values

n = np.arange(0,13357,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], cs2_warren),
        sit_pmw         =(["n","year"], cs2_sit),
        sit_sm          =(["n","year"], cs2_stroeve),
        sit_w99_uncertainty       =(["n","year"], cs2_warren_u),
        sit_pmw_uncertainty       =(["n","year"], cs2_sit_u),
        sit_sm_uncertainty        =(["n","year"], cs2_stroeve_u)  
        ),
    coords=dict(
        lon=(["n"],lon12),
        lat=(["n"],lat12),
        year = cs2_years,
        ),
    attrs=dict(description="Gridded CryoSat-2 CCI sea ice thickness [m] (SnowModel snowdensity"))
    
data.to_netcdf(direc+'SIT_CryoSat2_March_2011_2017_SnowModelsnowdens.nc')

#%%
"""Envisat Warren snowdensity"""

direc = '../../data/processed/envi/'
envi_sit = np.zeros((3045,10)); envi_sit[:,:] = np.nan
envi_warren = np.zeros((3045,10)); envi_warren[:,:] = np.nan
envi_stroeve = np.zeros((3045,10)); envi_stroeve[:,:] = np.nan

envi_sit_u = np.zeros((3045,10)); envi_sit_u[:,:] = np.nan
envi_warren_u = np.zeros((3045,10)); envi_warren_u[:,:] = np.nan
envi_stroeve_u = np.zeros((3045,10)); envi_stroeve_u[:,:] = np.nan
envi_years = np.zeros((10)).astype(int)
i=0
for filename in listdir(direc):
    if filename.endswith('03_warrendens.csv'):
        year = filename[8:12] #[14:18]
        f = pd.read_csv(direc+filename)
        envi_sit[:,i] = f.sea_ice_thickness_pmw.values
        envi_warren[:,i] = f.sea_ice_thickness_w99.values
        envi_stroeve[:,i] = f.sea_ice_thickness_snowmodel.values
        
        envi_sit_u[:,i] = f.sit_uncer_pmw.values
        envi_warren_u[:,i] = f.sit_uncer_w99.values
        envi_stroeve_u[:,i] = f.sit_uncer_sm.values
        envi_years[i] = year
        i=i+1

n = np.arange(0,3045,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], envi_warren),
        sit_pmw         =(["n","year"], envi_sit),
        sit_sm          =(["n","year"], envi_stroeve),
        sit_w99_uncertainty       =(["n","year"], envi_warren_u),
        sit_pmw_uncertainty       =(["n","year"], envi_sit_u),
        sit_sm_uncertainty        =(["n","year"], envi_stroeve_u)  
        ),
    coords=dict(
        lon=(["n"],lon25),
        lat=(["n"],lat25),
        year = envi_years,
        ),
    attrs=dict(description="Gridded Envisat sea ice thickness [m] (Warren snowdensity"))
    
data.to_netcdf(direc+'SIT_Envisat_March_2003_2012_Warrensnowdens.nc')

#%%
"""Envisat SnowModel snowdensity"""

direc = '../../data/processed/envi/'
envi_sit = np.zeros((3045,10)); envi_sit[:,:] = np.nan
envi_warren = np.zeros((3045,10)); envi_warren[:,:] = np.nan
envi_stroeve = np.zeros((3045,10)); envi_stroeve[:,:] = np.nan

envi_sit_u = np.zeros((3045,10)); envi_sit_u[:,:] = np.nan
envi_warren_u = np.zeros((3045,10)); envi_warren_u[:,:] = np.nan
envi_stroeve_u = np.zeros((3045,10)); envi_stroeve_u[:,:] = np.nan
envi_years = np.zeros((10)).astype(int)
i=0
for filename in listdir(direc):
    if filename.endswith('03_smdens.csv'):
        year = filename[8:12] #[14:18]
        f = pd.read_csv(direc+filename)
        envi_sit[:,i] = f.sea_ice_thickness_pmw.values
        envi_warren[:,i] = f.sea_ice_thickness_w99.values
        envi_stroeve[:,i] = f.sea_ice_thickness_snowmodel.values
        
        envi_sit_u[:,i] = f.sit_uncer_pmw.values
        envi_warren_u[:,i] = f.sit_uncer_w99.values
        envi_stroeve_u[:,i] = f.sit_uncer_sm.values
        envi_years[i] = year
        i=i+1

n = np.arange(0,3045,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], envi_warren),
        sit_pmw         =(["n","year"], envi_sit),
        sit_sm          =(["n","year"], envi_stroeve),
        sit_w99_uncertainty       =(["n","year"], envi_warren_u),
        sit_pmw_uncertainty       =(["n","year"], envi_sit_u),
        sit_sm_uncertainty        =(["n","year"], envi_stroeve_u)  
        ),
    coords=dict(
        lon=(["n"],lon25),
        lat=(["n"],lat25),
        year = envi_years,
        ),
    attrs=dict(description="Gridded Envisat sea ice thickness [m] (SnowModel snowdensity"))
    
data.to_netcdf(direc+'SIT_Envisat_March_2003_2012_SnowModelsnowdens.nc')

#%%
'''ICESat-2 Warren snowdensity'''
direc = '../../data/processed/is2/'
is2_warren = np.zeros((13357,2)); is2_warren[:,:] = np.nan
is2_warren_petty = np.zeros((13357,2)); is2_warren_petty[:,:] = np.nan
is2_warren_eff = np.zeros((13357,2)); is2_warren_eff[:,:] = np.nan
is2_pmw = np.zeros((13357,2)); is2_pmw[:,:] = np.nan
is2_pmw_eff = np.zeros((13357,2)); is2_pmw_eff[:,:] = np.nan
is2_pmw_petty = np.zeros((13357,2)); is2_pmw_petty[:,:] = np.nan
is2_sm = np.zeros((13357,2)); is2_sm[:,:] = np.nan
is2_sm_eff = np.zeros((13357,2)); is2_sm_eff[:,:] = np.nan
is2_sm_petty = np.zeros((13357,2)); is2_sm_petty[:,:] = np.nan
is2_nesosim_petty = np.zeros((13357,2)); is2_nesosim_petty[:,:] = np.nan

is2_warren_u = np.zeros((13357,2)); is2_warren_u[:,:] = np.nan
is2_warren_petty_u = np.zeros((13357,2)); is2_warren_petty_u[:,:] = np.nan
is2_warren_eff_u = np.zeros((13357,2)); is2_warren_eff_u[:,:] = np.nan
is2_pmw_u = np.zeros((13357,2)); is2_pmw_u[:,:] = np.nan
is2_pmw_eff_u = np.zeros((13357,2)); is2_pmw_eff_u[:,:] = np.nan
is2_pmw_petty_u = np.zeros((13357,2)); is2_pmw_petty_u[:,:] = np.nan
is2_sm_u = np.zeros((13357,2)); is2_sm_u[:,:] = np.nan
is2_sm_eff_u = np.zeros((13357,2)); is2_sm_eff_u[:,:] = np.nan
is2_sm_petty_u = np.zeros((13357,2)); is2_sm_petty_u[:,:] = np.nan
is2_nesosim_petty_u = np.zeros((13357,2)); is2_nesosim_petty_u[:,:] = np.nan
is2_years = np.zeros((2)).astype(int)
i=0
for filename in listdir(direc):
    if filename.endswith('03_warrendens.csv'):
        is_temp = pd.read_csv(direc+filename)
        is2_warren[:,i] = is_temp.sea_ice_thickness_w99.values
        is2_warren_petty[:,i] = is_temp.sea_ice_thickness_w99_petty.values
        is2_warren_eff[:,i] = is_temp.sea_ice_thickness_w99_eff.values
        is2_pmw[:,i] = is_temp.sea_ice_thickness_pmw.values
        is2_pmw_eff[:,i] = is_temp.sea_ice_thickness_pmw_eff.values
        is2_pmw_petty[:,i] = is_temp.sea_ice_thickness_pmw_petty.values  
        is2_sm[:,i] = is_temp.sea_ice_thickness_snowmodel.values 
        is2_sm_eff[:,i] = is_temp.sea_ice_thickness_snowmodel_eff.values 
        is2_sm_petty[:,i] = is_temp.sea_ice_thickness_snowmodel_petty.values 
        is2_nesosim_petty[:,i] = is_temp.sea_ice_thickness_nesosim_petty.values 
        
        is2_warren_u[:,i] = is_temp.sit_uncer_w99.values
        is2_warren_petty_u[:,i] = is_temp.sit_uncer_w99_petty.values
        is2_warren_eff_u[:,i] = is_temp.sit_uncer_w99_eff.values
        is2_pmw_u[:,i] = is_temp.sit_uncer_pmw.values
        is2_pmw_eff_u[:,i] = is_temp.sit_uncer_pmw_eff.values
        is2_pmw_petty_u[:,i] = is_temp.sit_uncer_pmw_petty.values  
        is2_sm_u[:,i] = is_temp.sit_uncer_snowmodel.values 
        is2_sm_eff_u[:,i] = is_temp.sit_uncer_snowmodel_eff.values 
        is2_sm_petty_u[:,i] = is_temp.sit_uncer_snowmodel_petty.values 
        is2_nesosim_petty_u[:,i] = is_temp.sit_uncer_nesosim_petty.values 
        is2_years[i] = filename[8:12]
        i=i+1

n = np.arange(0,13357,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], is2_warren),
        sit_w99_sig     =(["n","year"], is2_warren_eff),
        sit_w99_piece   =(["n","year"], is2_warren_petty),
        sit_pmw         =(["n","year"], is2_pmw),
        sit_pmw_sig     =(["n","year"], is2_pmw_eff),
        sit_pmw_piece   =(["n","year"], is2_pmw_petty),
        sit_sm          =(["n","year"], is2_sm),
        sit_sm_sig      =(["n","year"], is2_sm_eff),
        sit_sm_piece    =(["n","year"], is2_sm_petty),
        sit_w99_uncertainty       =(["n","year"], is2_warren_u),
        sit_w99_sig_uncertainty   =(["n","year"], is2_warren_eff_u),
        sit_w99_piece_uncertainty =(["n","year"], is2_warren_petty_u),
        sit_pmw_uncertainty       =(["n","year"], is2_pmw_u),
        sit_pmw_sig_uncertainty   =(["n","year"], is2_pmw_eff_u),
        sit_pmw_piece_uncertainty =(["n","year"], is2_pmw_petty_u),
        sit_sm_uncertainty        =(["n","year"], is2_sm_u),
        sit_sm_sig_uncertainty    =(["n","year"], is2_sm_eff_u),
        sit_sm_piece_uncertainty  =(["n","year"], is2_sm_petty_u)   
        ),
    coords=dict(
        lon=(["n"],lon12),
        lat=(["n"],lat12),
        year = is2_years,
        ),
    attrs=dict(description="Gridded ICESat-2 sea ice thickness [m] (Warren snowdensity)"))
    
data.to_netcdf(direc+'SIT_ICESat2_March_2019_2020_Warrensnowdens.nc')

#%%
'''ICESat-2 SnowModel snowdensity'''
direc = '../../data/processed/is2/'
is2_warren = np.zeros((13357,2)); is2_warren[:,:] = np.nan
is2_warren_petty = np.zeros((13357,2)); is2_warren_petty[:,:] = np.nan
is2_warren_eff = np.zeros((13357,2)); is2_warren_eff[:,:] = np.nan
is2_pmw = np.zeros((13357,2)); is2_pmw[:,:] = np.nan
is2_pmw_eff = np.zeros((13357,2)); is2_pmw_eff[:,:] = np.nan
is2_pmw_petty = np.zeros((13357,2)); is2_pmw_petty[:,:] = np.nan
is2_sm = np.zeros((13357,2)); is2_sm[:,:] = np.nan
is2_sm_eff = np.zeros((13357,2)); is2_sm_eff[:,:] = np.nan
is2_sm_petty = np.zeros((13357,2)); is2_sm_petty[:,:] = np.nan

is2_warren_u = np.zeros((13357,2)); is2_warren_u[:,:] = np.nan
is2_warren_petty_u = np.zeros((13357,2)); is2_warren_petty_u[:,:] = np.nan
is2_warren_eff_u = np.zeros((13357,2)); is2_warren_eff_u[:,:] = np.nan
is2_pmw_u = np.zeros((13357,2)); is2_pmw_u[:,:] = np.nan
is2_pmw_eff_u = np.zeros((13357,2)); is2_pmw_eff_u[:,:] = np.nan
is2_pmw_petty_u = np.zeros((13357,2)); is2_pmw_petty_u[:,:] = np.nan
is2_sm_u = np.zeros((13357,2)); is2_sm_u[:,:] = np.nan
is2_sm_eff_u = np.zeros((13357,2)); is2_sm_eff_u[:,:] = np.nan
is2_sm_petty_u = np.zeros((13357,2)); is2_sm_petty_u[:,:] = np.nan
is2_years = np.zeros((2)).astype(int)
i=0
for filename in listdir(direc):
    if filename.endswith('03_smdens.csv'):
        is_temp = pd.read_csv(direc+filename)
        is2_warren[:,i] = is_temp.sea_ice_thickness_w99.values
        is2_warren_petty[:,i] = is_temp.sea_ice_thickness_w99_petty.values
        is2_warren_eff[:,i] = is_temp.sea_ice_thickness_w99_eff.values
        is2_pmw[:,i] = is_temp.sea_ice_thickness_pmw.values
        is2_pmw_eff[:,i] = is_temp.sea_ice_thickness_pmw_eff.values
        is2_pmw_petty[:,i] = is_temp.sea_ice_thickness_pmw_petty.values  
        is2_sm[:,i] = is_temp.sea_ice_thickness_snowmodel.values 
        is2_sm_eff[:,i] = is_temp.sea_ice_thickness_snowmodel_eff.values 
        is2_sm_petty[:,i] = is_temp.sea_ice_thickness_snowmodel_petty.values 
        
        is2_warren_u[:,i] = is_temp.sit_uncer_w99.values
        is2_warren_petty_u[:,i] = is_temp.sit_uncer_w99_petty.values
        is2_warren_eff_u[:,i] = is_temp.sit_uncer_w99_eff.values
        is2_pmw_u[:,i] = is_temp.sit_uncer_pmw.values
        is2_pmw_eff_u[:,i] = is_temp.sit_uncer_pmw_eff.values
        is2_pmw_petty_u[:,i] = is_temp.sit_uncer_pmw_petty.values  
        is2_sm_u[:,i] = is_temp.sit_uncer_snowmodel.values 
        is2_sm_eff_u[:,i] = is_temp.sit_uncer_snowmodel_eff.values 
        is2_sm_petty_u[:,i] = is_temp.sit_uncer_snowmodel_petty.values 
        is2_years[i] = filename[8:12]
        i=i+1

n = np.arange(0,13357,1).astype(int)

data = xr.Dataset(
    data_vars = dict(
        sit_w99         =(["n","year"], is2_warren),
        sit_w99_sig     =(["n","year"], is2_warren_eff),
        sit_w99_piece   =(["n","year"], is2_warren_petty),
        sit_pmw         =(["n","year"], is2_pmw),
        sit_pmw_sig     =(["n","year"], is2_pmw_eff),
        sit_pmw_piece   =(["n","year"], is2_pmw_petty),
        sit_sm          =(["n","year"], is2_sm),
        sit_sm_sig      =(["n","year"], is2_sm_eff),
        sit_sm_piece    =(["n","year"], is2_sm_petty),
        sit_w99_uncertainty       =(["n","year"], is2_warren_u),
        sit_w99_sig_uncertainty   =(["n","year"], is2_warren_eff_u),
        sit_w99_piece_uncertainty =(["n","year"], is2_warren_petty_u),
        sit_pmw_uncertainty       =(["n","year"], is2_pmw_u),
        sit_pmw_sig_uncertainty   =(["n","year"], is2_pmw_eff_u),
        sit_pmw_piece_uncertainty =(["n","year"], is2_pmw_petty_u),
        sit_sm_uncertainty        =(["n","year"], is2_sm_u),
        sit_sm_sig_uncertainty    =(["n","year"], is2_sm_eff_u),
        sit_sm_piece_uncertainty  =(["n","year"], is2_sm_petty_u)    
        ),
    coords=dict(
        lon=(["n"],lon12),
        lat=(["n"],lat12),
        year = is2_years,
        ),
    attrs=dict(description="Gridded ICESat-2 sea ice thickness [m] (SnowModel snowdensity)"))


    
data.to_netcdf(direc+'SIT_ICESat2_March_2019_2020_SnowModelsnowdens.nc')
