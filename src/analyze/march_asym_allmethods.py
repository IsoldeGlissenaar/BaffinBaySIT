# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 13:44:51 2021

@author: IsoldeGlissenaar

Time series of asymmetry in SIT for all satellites and processing methods. 
Plots SIT asymmetry time series envelopes for all processing methods and satellites.
Prints SIT asymmetry trends over ICESat & CryoSat-2 LARM for all combinations
of processing methods.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir
from shapely.geometry import Point # Point class
from shapely.geometry import shape # shape() is a function to convert geo objects through the interface
import shapefile
import xarray as xr
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

#%%
#Import SIT products
icesat_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is/SIT_ICESat_March_2003_2009_SnowModelsnowdens.nc')
icesat_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is/SIT_ICESat_March_2003_2009_Warrensnowdens.nc')

envi_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/envi/SIT_Envisat_March_2003_2012_SnowModelsnowdens.nc')
envi_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/envi/SIT_Envisat_March_2003_2012_Warrensnowdens.nc')

cs2_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2/SIT_CryoSat2_March_2011_2020_SnowModelsnowdens.nc')
cs2_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2/SIT_CryoSat2_March_2011_2020_Warrensnowdens.nc')

cs2_cci_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2_cci/SIT_CryoSat2_March_2011_2017_SnowModelsnowdens.nc')
cs2_cci_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2_cci/SIT_CryoSat2_March_2011_2017_Warrensnowdens.nc')

icesat2_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is2/SIT_ICESat2_March_2019_2020_SnowModelsnowdens.nc')
icesat2_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is2/SIT_ICESat2_March_2019_2020_Warrensnowdens.nc')
 

#%%
''' Import CIS charts data'''

#Load CIS chart 
cis_sit = np.load('../../data/interim/CIS/march_mean_sit.npy')
latcis = np.load('../../data/interim/CIS/march_mean_lat.npy')
loncis = np.load('../../data/interim/CIS/march_mean_lon.npy')
cis_years = np.load('../../data/interim/CIS/march_mean_yr.npy')


#%%
#Find east-west location gridcells
shp = shapefile.Reader(docs+'/QGIS/4quad.shp') #open the shapefile
all_shapes = shp.shapes() # get all the polygons
all_records = shp.records()     

#Get east/west location gridcells on 25 km grid (ICESat & Envisat)
len_f = icesat_sm.dims['n']
eastwest25 = np.zeros(len_f)     
for i in range(len_f):
    pt = (icesat_sm.lon.values[i], icesat_sm.lat.values[i])  
    for k in range (len(all_shapes)):             
        boundary = all_shapes[k]
        if Point(pt).within(shape(boundary)): # make a point and see if it's in the polygon
            eastwest25[i] = all_records[k][0]

#Get east/west location gridcells on 12 km grid (ICESat-2 & CryoSat-2)
len_f = cs2_sm.dims['n']
eastwest12 = np.zeros(len_f)     
for i in range(len_f):
    pt = (cs2_sm.lon.values[i],cs2_sm.lat.values[i])  
    for k in range (len(all_shapes)):             
        boundary = all_shapes[k]
        if Point(pt).within(shape(boundary)): # make a point and see if it's in the polygon
            eastwest12[i] = all_records[k][0]
            
#Get east/west location gridcells on CIS chart grid
len_f = len(cis_sit[0,:])
eastwest_cis = np.zeros(len_f)     
for i in range(len_f):
    pt = (loncis[i],latcis[i])  
    for k in range (len(all_shapes)):             
        boundary = all_shapes[k]
        if Point(pt).within(shape(boundary)): # make a point and see if it's in the polygon
            eastwest_cis[i] = all_records[k][0]
            
            
#%%

def asymmetry(sit, eastwest):
    west = sit[(eastwest==3)|(eastwest==4),:]
    mean_west = np.nanmean(west,axis=0)
    east = sit[(eastwest==1)|(eastwest==2),:]
    mean_east = np.nanmean(east,axis=0)
    return (mean_west-mean_east)

cis_asym = asymmetry(np.transpose(cis_sit), eastwest_cis)
     
#%%
"""ICESat"""

eastwest25_2d = np.transpose(np.tile(eastwest25, (7,1)))
is_warren_west = icesat_warren.where((eastwest25_2d==3)|(eastwest25_2d==4)).mean('n')
is_warren_east = icesat_warren.where((eastwest25_2d==1)|(eastwest25_2d==2)).mean('n')
is_sm_west = icesat_sm.where((eastwest25_2d==3)|(eastwest25_2d==4)).mean('n')
is_sm_east = icesat_sm.where((eastwest25_2d==1)|(eastwest25_2d==2)).mean('n')

is_asym = np.array([is_warren_west.sit_w99.values-is_warren_east.sit_w99.values,
                    is_warren_west.sit_w99_sig.values-is_warren_east.sit_w99_sig.values,
                    is_warren_west.sit_w99_piece.values-is_warren_east.sit_w99_piece.values,
                    is_warren_west.sit_pmw.values-is_warren_east.sit_pmw.values,
                    is_warren_west.sit_pmw_sig.values-is_warren_east.sit_pmw_sig.values,
                    is_warren_west.sit_pmw_piece.values-is_warren_east.sit_pmw_piece.values,
                    is_warren_west.sit_sm.values-is_warren_east.sit_sm.values,
                    is_warren_west.sit_sm_sig.values-is_warren_east.sit_sm_sig.values,
                    is_warren_west.sit_sm_piece.values-is_warren_east.sit_sm_piece.values,
                    
                    is_sm_west.sit_w99.values-is_sm_east.sit_w99.values,
                    is_sm_west.sit_w99_sig.values-is_sm_east.sit_w99_sig.values,
                    is_sm_west.sit_w99_piece.values-is_sm_east.sit_w99_piece.values,
                    is_sm_west.sit_pmw.values-is_sm_east.sit_pmw.values,
                    is_sm_west.sit_pmw_sig.values-is_sm_east.sit_pmw_sig.values,
                    is_sm_west.sit_pmw_piece.values-is_sm_east.sit_pmw_piece.values,
                    is_sm_west.sit_sm.values-is_sm_east.sit_sm.values,
                    is_sm_west.sit_sm_sig.values-is_sm_east.sit_sm_sig.values,
                    is_sm_west.sit_sm_piece.values-is_sm_east.sit_sm_piece.values])

is_asym_u = np.array([is_warren_west.sit_w99_uncertainty.values+is_warren_east.sit_w99_uncertainty.values,
                      is_warren_west.sit_w99_sig_uncertainty.values+is_warren_east.sit_w99_sig_uncertainty.values,
                      is_warren_west.sit_w99_piece_uncertainty.values+is_warren_east.sit_w99_piece_uncertainty.values,
                      is_warren_west.sit_pmw_uncertainty.values+is_warren_east.sit_pmw_uncertainty.values,
                      is_warren_west.sit_pmw_sig_uncertainty.values+is_warren_east.sit_pmw_sig_uncertainty.values,
                      is_warren_west.sit_pmw_piece_uncertainty.values+is_warren_east.sit_pmw_piece_uncertainty.values,
                      is_warren_west.sit_sm_uncertainty.values+is_warren_east.sit_sm_uncertainty.values,
                      is_warren_west.sit_sm_sig_uncertainty.values+is_warren_east.sit_sm_sig_uncertainty.values,
                      is_warren_west.sit_sm_piece_uncertainty.values+is_warren_east.sit_sm_piece_uncertainty.values,
                      
                      is_sm_west.sit_w99_uncertainty.values+is_sm_east.sit_w99_uncertainty.values,
                      is_sm_west.sit_w99_sig_uncertainty.values+is_sm_east.sit_w99_sig_uncertainty.values,
                      is_sm_west.sit_w99_piece_uncertainty.values+is_sm_east.sit_w99_piece_uncertainty.values,
                      is_sm_west.sit_pmw_uncertainty.values+is_sm_east.sit_pmw_uncertainty.values,
                      is_sm_west.sit_pmw_sig_uncertainty.values+is_sm_east.sit_pmw_sig_uncertainty.values,
                      is_sm_west.sit_pmw_piece_uncertainty.values+is_sm_east.sit_pmw_piece_uncertainty.values,
                      is_sm_west.sit_sm_uncertainty.values+is_sm_east.sit_sm_uncertainty.values,
                      is_sm_west.sit_sm_sig_uncertainty.values+is_sm_east.sit_sm_sig_uncertainty.values,
                      is_sm_west.sit_sm_piece_uncertainty.values+is_sm_east.sit_sm_piece_uncertainty.values])

min_methods = np.nanmin(is_asym, axis=0)
max_methods = np.nanmax(is_asym, axis=0)

is_methods_max = is_asym+is_asym_u
max_methods_u = np.nanmax(is_methods_max, axis=0)

is_methods_min = is_asym-is_asym_u
min_methods_u = np.nanmin(is_methods_min, axis=0)

#%%
"""Envisat"""

eastwest25_2d = np.transpose(np.tile(eastwest25, (10,1)))
envi_warren_west = envi_warren.where((eastwest25_2d==3)|(eastwest25_2d==4)).mean('n')
envi_warren_east = envi_warren.where((eastwest25_2d==1)|(eastwest25_2d==2)).mean('n')
envi_sm_west = envi_sm.where((eastwest25_2d==3)|(eastwest25_2d==4)).mean('n')
envi_sm_east = envi_sm.where((eastwest25_2d==1)|(eastwest25_2d==2)).mean('n')

envi_asym = np.array([envi_warren_west.sit_w99.values-envi_warren_east.sit_w99.values,
                      envi_warren_west.sit_pmw.values-envi_warren_east.sit_pmw.values,
                      envi_warren_west.sit_sm.values-envi_warren_east.sit_sm.values,
                      envi_sm_west.sit_w99.values-envi_sm_east.sit_w99.values,
                      envi_sm_west.sit_pmw.values-envi_sm_east.sit_pmw.values,
                      envi_sm_west.sit_sm.values-envi_sm_east.sit_sm.values])

envi_asym_u = np.array([envi_warren_west.sit_w99_uncertainty.values+envi_warren_east.sit_w99_uncertainty.values,
                        envi_warren_west.sit_pmw_uncertainty.values+envi_warren_east.sit_pmw_uncertainty.values,
                        envi_warren_west.sit_sm_uncertainty.values+envi_warren_east.sit_sm_uncertainty.values,
                      
                        envi_sm_west.sit_w99_uncertainty.values+envi_sm_east.sit_w99_uncertainty.values,
                        envi_sm_west.sit_pmw_uncertainty.values+envi_sm_east.sit_pmw_uncertainty.values,
                        envi_sm_west.sit_sm_uncertainty.values+envi_sm_east.sit_sm_uncertainty.values])

min_methods_envi = np.nanmin(envi_asym, axis=0)
max_methods_envi = np.nanmax(envi_asym, axis=0)

envi_methods_max = envi_asym+envi_asym_u
max_methods_envi_u = np.nanmax(envi_methods_max, axis=0)

envi_methods_min = envi_asym-envi_asym_u
min_methods_envi_u = np.nanmin(envi_methods_min, axis=0)

#%%
"""CryoSat-2"""

eastwest12_2d = np.transpose(np.tile(eastwest12, (10,1)))
cs2_warren_west = cs2_warren.where((eastwest12_2d==3)|(eastwest12_2d==4)).mean('n')
cs2_warren_east = cs2_warren.where((eastwest12_2d==1)|(eastwest12_2d==2)).mean('n')
cs2_sm_west = cs2_sm.where((eastwest12_2d==3)|(eastwest12_2d==4)).mean('n')
cs2_sm_east = cs2_sm.where((eastwest12_2d==1)|(eastwest12_2d==2)).mean('n')

cs2_asym = np.array([cs2_warren_west.sit_w99.values-cs2_warren_east.sit_w99.values,
                     cs2_warren_west.sit_pmw.values-cs2_warren_east.sit_pmw.values,
                     cs2_warren_west.sit_sm.values-cs2_warren_east.sit_sm.values,
                     cs2_sm_west.sit_w99.values-cs2_sm_east.sit_w99.values,
                     cs2_sm_west.sit_pmw.values-cs2_sm_east.sit_pmw.values,
                     cs2_sm_west.sit_sm.values-cs2_sm_east.sit_sm.values])

cs2_asym_u = np.array([cs2_warren_west.sit_w99_uncertainty.values+cs2_warren_east.sit_w99_uncertainty.values,
                       cs2_warren_west.sit_pmw_uncertainty.values+cs2_warren_east.sit_pmw_uncertainty.values,
                       cs2_warren_west.sit_sm_uncertainty.values+cs2_warren_east.sit_sm_uncertainty.values,
                      
                       cs2_sm_west.sit_w99_uncertainty.values+cs2_sm_east.sit_w99_uncertainty.values,
                       cs2_sm_west.sit_pmw_uncertainty.values+cs2_sm_east.sit_pmw_uncertainty.values,
                       cs2_sm_west.sit_sm_uncertainty.values+cs2_sm_east.sit_sm_uncertainty.values])

min_methods_cs2 = np.nanmin(cs2_asym, axis=0)
max_methods_cs2 = np.nanmax(cs2_asym, axis=0)

cs2_methods_max = cs2_asym+cs2_asym_u
max_methods_cs2_u = np.nanmax(cs2_methods_max, axis=0)

cs2_methods_min = cs2_asym-cs2_asym_u
min_methods_cs2_u = np.nanmin(cs2_methods_min, axis=0)

#%%
"""CryoSat-2 CCI"""

eastwest12_2d = np.transpose(np.tile(eastwest12, (10,1)))
cs2_cci_warren_west = cs2_cci_warren.where((eastwest12_2d==3)|(eastwest12_2d==4)).mean('n')
cs2_cci_warren_east = cs2_cci_warren.where((eastwest12_2d==1)|(eastwest12_2d==2)).mean('n')
cs2_cci_sm_west = cs2_cci_sm.where((eastwest12_2d==3)|(eastwest12_2d==4)).mean('n')
cs2_cci_sm_east = cs2_cci_sm.where((eastwest12_2d==1)|(eastwest12_2d==2)).mean('n')

cs2_cci_asym = np.array([cs2_cci_warren_west.sit_w99.values-cs2_cci_warren_east.sit_w99.values,
                         cs2_cci_warren_west.sit_pmw.values-cs2_cci_warren_east.sit_pmw.values,
                         cs2_cci_warren_west.sit_sm.values-cs2_cci_warren_east.sit_sm.values,
                         cs2_cci_sm_west.sit_w99.values-cs2_cci_sm_east.sit_w99.values,
                         cs2_cci_sm_west.sit_pmw.values-cs2_cci_sm_east.sit_pmw.values,
                         cs2_cci_sm_west.sit_sm.values-cs2_cci_sm_east.sit_sm.values])

cs2_cci_asym_u = np.array([cs2_cci_warren_west.sit_w99_uncertainty.values+cs2_cci_warren_east.sit_w99_uncertainty.values,
                           cs2_cci_warren_west.sit_pmw_uncertainty.values+cs2_cci_warren_east.sit_pmw_uncertainty.values,
                           cs2_cci_warren_west.sit_sm_uncertainty.values+cs2_cci_warren_east.sit_sm_uncertainty.values,
                      
                           cs2_cci_sm_west.sit_w99_uncertainty.values+cs2_cci_sm_east.sit_w99_uncertainty.values,
                           cs2_cci_sm_west.sit_pmw_uncertainty.values+cs2_cci_sm_east.sit_pmw_uncertainty.values,
                           cs2_cci_sm_west.sit_sm_uncertainty.values+cs2_cci_sm_east.sit_sm_uncertainty.values])

min_methods_cs2_cci = np.nanmin(cs2_cci_asym, axis=0)
max_methods_cs2_cci = np.nanmax(cs2_cci_asym, axis=0)

cs2_cci_methods_max = cs2_cci_asym+cs2_cci_asym_u
max_methods_cs2_cci_u = np.nanmax(cs2_cci_methods_max, axis=0)

cs2_cci_methods_min = cs2_cci_asym-cs2_cci_asym_u
min_methods_cs2_cci_u = np.nanmin(cs2_cci_methods_min, axis=0)

#%%
"""ICESat-2"""

eastwest12_2d = np.transpose(np.tile(eastwest12, (2,1)))
is2_warren_west = icesat2_warren.where((eastwest12_2d==3)|(eastwest12_2d==4)).mean('n')
is2_warren_east = icesat2_warren.where((eastwest12_2d==1)|(eastwest12_2d==2)).mean('n')
is2_sm_west = icesat2_sm.where((eastwest12_2d==3)|(eastwest12_2d==4)).mean('n')
is2_sm_east = icesat2_sm.where((eastwest12_2d==1)|(eastwest12_2d==2)).mean('n')

is2_asym = np.array([is2_warren_west.sit_w99.values-is2_warren_east.sit_w99.values,
                     is2_warren_west.sit_w99_sig.values-is2_warren_east.sit_w99_sig.values,
                     is2_warren_west.sit_w99_piece.values-is2_warren_east.sit_w99_piece.values,
                     is2_warren_west.sit_pmw.values-is2_warren_east.sit_pmw.values,
                     is2_warren_west.sit_pmw_sig.values-is2_warren_east.sit_pmw_sig.values,
                     is2_warren_west.sit_pmw_piece.values-is2_warren_east.sit_pmw_piece.values,
                     is2_warren_west.sit_sm.values-is2_warren_east.sit_sm.values,
                     is2_warren_west.sit_sm_sig.values-is2_warren_east.sit_sm_sig.values,
                     is2_warren_west.sit_sm_piece.values-is2_warren_east.sit_sm_piece.values,
                     
                     is2_sm_west.sit_w99.values-is2_sm_east.sit_w99.values,
                     is2_sm_west.sit_w99_sig.values-is2_sm_east.sit_w99_sig.values,
                     is2_sm_west.sit_w99_piece.values-is2_sm_east.sit_w99_piece.values,
                     is2_sm_west.sit_pmw.values-is2_sm_east.sit_pmw.values,
                     is2_sm_west.sit_pmw_sig.values-is2_sm_east.sit_pmw_sig.values,
                     is2_sm_west.sit_pmw_piece.values-is2_sm_east.sit_pmw_piece.values,
                     is2_sm_west.sit_sm.values-is2_sm_east.sit_sm.values,
                     is2_sm_west.sit_sm_sig.values-is2_sm_east.sit_sm_sig.values,
                     is2_sm_west.sit_sm_piece.values-is2_sm_east.sit_sm_piece.values])

is2_asym_u = np.array([is2_warren_west.sit_w99_uncertainty.values+is2_warren_east.sit_w99_uncertainty.values,
                       is2_warren_west.sit_w99_sig_uncertainty.values+is2_warren_east.sit_w99_sig_uncertainty.values,
                       is2_warren_west.sit_w99_piece_uncertainty.values+is2_warren_east.sit_w99_piece_uncertainty.values,
                       is2_warren_west.sit_pmw_uncertainty.values+is2_warren_east.sit_pmw_uncertainty.values,
                       is2_warren_west.sit_pmw_sig_uncertainty.values+is2_warren_east.sit_pmw_sig_uncertainty.values,
                       is2_warren_west.sit_pmw_piece_uncertainty.values+is2_warren_east.sit_pmw_piece_uncertainty.values,
                       is2_warren_west.sit_sm_uncertainty.values+is2_warren_east.sit_sm_uncertainty.values,
                       is2_warren_west.sit_sm_sig_uncertainty.values+is2_warren_east.sit_sm_sig_uncertainty.values,
                       is2_warren_west.sit_sm_piece_uncertainty.values+is2_warren_east.sit_sm_piece_uncertainty.values,
                      
                       is2_sm_west.sit_w99_uncertainty.values+is2_sm_east.sit_w99_uncertainty.values,
                       is2_sm_west.sit_w99_sig_uncertainty.values+is2_sm_east.sit_w99_sig_uncertainty.values,
                       is2_sm_west.sit_w99_piece_uncertainty.values+is2_sm_east.sit_w99_piece_uncertainty.values,
                       is2_sm_west.sit_pmw_uncertainty.values+is2_sm_east.sit_pmw_uncertainty.values,
                       is2_sm_west.sit_pmw_sig_uncertainty.values+is2_sm_east.sit_pmw_sig_uncertainty.values,
                       is2_sm_west.sit_pmw_piece_uncertainty.values+is2_sm_east.sit_pmw_piece_uncertainty.values,
                       is2_sm_west.sit_sm_uncertainty.values+is2_sm_east.sit_sm_uncertainty.values,
                       is2_sm_west.sit_sm_sig_uncertainty.values+is2_sm_east.sit_sm_sig_uncertainty.values,
                       is2_sm_west.sit_sm_piece_uncertainty.values+is2_sm_east.sit_sm_piece_uncertainty.values])

min_methods_is2 = np.nanmin(is2_asym, axis=0)
max_methods_is2 = np.nanmax(is2_asym, axis=0)

is2_methods_max = is2_asym+is2_asym_u
max_methods_is2_u = np.nanmax(is2_methods_max, axis=0)

is2_methods_min = is2_asym-is2_asym_u
min_methods_is2_u = np.nanmin(is2_methods_min, axis=0)

#%%
#Plot
fig, ax = plt.subplots(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)

plt.plot(icesat_warren.year, min_methods, label='ICESat')
plt.fill_between(icesat_warren.year, min_methods, max_methods, alpha=0.5)
plt.plot(icesat_warren.year, max_methods, c='C0')

plt.plot(envi_warren.year, min_methods_envi, c='C1', label='Envisat')
plt.fill_between(envi_warren.year, min_methods_envi, max_methods_envi, alpha=0.5)
plt.plot(envi_warren.year, max_methods_envi,c='C1')

plt.plot(cs2_warren.year, min_methods_cs2, c='#43853B', label='CryoSat-2 LARM')
plt.fill_between(cs2_warren.year, min_methods_cs2, max_methods_cs2, alpha=0.5)
plt.plot(cs2_warren.year, max_methods_cs2,c='#43853B')
plt.plot(cs2_warren.year, min_methods_cs2_cci, c='#60D652', label='CryoSat-2 CCI')
plt.fill_between(cs2_warren.year, min_methods_cs2_cci, max_methods_cs2_cci, color='#60D652', alpha=0.5)
plt.plot(cs2_warren.year, max_methods_cs2_cci,c='#60D652')

plt.plot(icesat2_warren.year, min_methods_is2,  c='C3', label='ICESat-2')
plt.fill_between(icesat2_warren.year, min_methods_is2, max_methods_is2, color='C3', alpha=0.5)
plt.plot(icesat2_warren.year, max_methods_is2, c='C3')

plt.plot(cis_years, cis_asym, c='C4', linestyle='--', label='CIS')

plt.axhline(0, c='k', linestyle='--', linewidth=1)
plt.grid(':', linewidth=0.5)
lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylabel('west-east [m]')
plt.xlim(left=2000)
plt.xticks(np.arange(2000,2022,2), rotation=45)
ax.set_xticks(np.arange(2001,2021,2), minor=True)
# plt.savefig('../../figures/fig10.png',bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)

#%%

from scipy import stats
## Trends diff processing techniques

''' W99 density '''

print('--------')
print('Warren snowdensity')
## Warren
years = np.arange(2003,2021,1)
sit_warren = np.full(len(years), np.nan)
sit_warren[icesat_warren.year.values-years[0]] = is_warren_west.sit_w99.values-is_warren_east.sit_w99.values
sit_warren[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_w99.values-cs2_warren_east.sit_w99.values
# sit_warren[is2_years-years[0]] = is2_warren_asym
x=years[:]
y=sit_warren[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_warren, intercept, r_value, p_value_warren, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99: ',np.round(trend_warren*100,2),' cm/yr (p=',np.round(p_value_warren,2),')')            

## Warren effsd
years = np.arange(2003,2021,1)
sit_w99eff = np.full(len(years), np.nan)
sit_w99eff[icesat_warren.year.values-years[0]] = is_warren_west.sit_w99_sig.values-is_warren_east.sit_w99_sig.values
sit_w99eff[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_w99.values-cs2_warren_east.sit_w99.values
# sit_w99eff[is2_years-years[0]] = is2_warren_e_asym
x=years[:]
y=sit_w99eff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_w99eff, intercept, r_value, p_value_w99eff, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99 Effsd: ',np.round(trend_w99eff*100,2),' cm/yr (p=',np.round(p_value_w99eff,2),')')        

## Warren Petty
years = np.arange(2003,2021,1)
sit_w99petty = np.full(len(years), np.nan)
sit_w99petty[icesat_warren.year.values-years[0]] = is_warren_west.sit_w99_piece.values-is_warren_east.sit_w99_piece.values
sit_w99petty[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_w99.values-cs2_warren_east.sit_w99.values
# sit_w99petty[is2_years-years[0]] = is2_warren_pet_asym
x=years[:]
y=sit_w99petty[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_w99petty, intercept, r_value, p_value_w99petty, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99 Petty: ',np.round(trend_w99petty*100,2),' cm/yr (p=',np.round(p_value_w99petty,2),')')        


## PMW
years = np.arange(2003,2021,1)
sit_pmw = np.full(len(years), np.nan)
sit_pmw[icesat_warren.year.values-years[0]] = is_warren_west.sit_pmw.values-is_warren_east.sit_pmw.values
sit_pmw[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_pmw.values-cs2_warren_east.sit_pmw.values
# sit_pmw[is2_years-years[0]] = is2_pmw_asym
x=years[:]
y=sit_pmw[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmw, intercept_pmw, r_value, p_value_pmw, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW: ',np.round(trend_pmw*100,2),' cm/yr (p=',np.round(p_value_pmw,2),')')            

## PMW effsd
years = np.arange(2003,2021,1)
sit_pmweff = np.full(len(years), np.nan)
sit_pmweff[icesat_warren.year.values-years[0]] = is_warren_west.sit_pmw_sig.values-is_warren_east.sit_pmw_sig.values
sit_pmweff[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_pmw.values-cs2_warren_east.sit_pmw.values
# sit_pmweff[is2_years-years[0]] = is2_pmw_e_asym
x=years[:]
y=sit_pmweff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmweff, intercept, r_value, p_value_pmweff, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW Effsd: ',np.round(trend_pmweff*100,2),' cm/yr (p=',np.round(p_value_pmweff,2),')')        

# PMW Petty
years = np.arange(2003,2021,1)
sit_pmwpetty = np.full(len(years), np.nan)
sit_pmwpetty[icesat_warren.year.values-years[0]] = is_warren_west.sit_pmw_piece.values-is_warren_east.sit_pmw_piece.values
sit_pmwpetty[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_pmw.values-cs2_warren_east.sit_pmw.values
# sit_pmwpetty[is2_years-years[0]] = is2_pmw_pet_asym
x=years[:]
y=sit_pmwpetty[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmwpetty, intercept_pmwpetty, r_value, p_value_pmwpetty, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW Petty: ',np.round(trend_pmwpetty*100,2),' cm/yr (p=',np.round(p_value_pmwpetty,2),')')        


## Stroeve
years = np.arange(2003,2021,1)
sit_stroe = np.full(len(years), np.nan)
sit_stroe[icesat_warren.year.values-years[0]] = is_warren_west.sit_sm.values-is_warren_east.sit_sm.values
sit_stroe[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_sm.values-cs2_warren_east.sit_sm.values
x=years[:]
y=sit_stroe[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroe, intercept_stroe, r_value, p_value_stroe, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve: ',np.round(trend_stroe*100,2),' cm/yr (p=',np.round(p_value_stroe,2),')')            

## Stroeve effsd
years = np.arange(2003,2021,1)
sit_stroeeff = np.full(len(years), np.nan)
sit_stroeeff[icesat_warren.year.values-years[0]] = is_warren_west.sit_sm_sig.values-is_warren_east.sit_sm_sig.values
sit_stroeeff[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_sm.values-cs2_warren_east.sit_sm.values
x=years[:]
y=sit_stroeeff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroeeff, intercept, r_value, p_value_stroeeff, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve Effsd: ',np.round(trend_stroeeff*100,2),' cm/yr (p=',np.round(p_value_stroeeff,2),')')            

## Stroeve Petty
years = np.arange(2003,2021,1)
sit_stroepett = np.full(len(years), np.nan)
sit_stroepett[icesat_warren.year.values-years[0]] = is_warren_west.sit_sm_piece.values-is_warren_east.sit_sm_piece.values
sit_stroepett[cs2_warren.year.values-years[0]] = cs2_warren_west.sit_sm.values-cs2_warren_east.sit_sm.values
x=years[:]
y=sit_stroepett[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroepett, intercept, r_value, p_value_stroepett, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve Petty: ',np.round(trend_stroepett*100,2),' cm/yr (p=',np.round(p_value_stroepett,2),')')  

#%%
''' SnowModel density '''

print('--------')
print('SnowModel snowdensity')
## Warren
years = np.arange(2003,2021,1)
sit_warren = np.full(len(years), np.nan)
sit_warren[icesat_sm.year.values-years[0]] = is_sm_west.sit_w99.values-is_sm_east.sit_w99.values
sit_warren[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_w99.values-cs2_sm_east.sit_w99.values
# sit_sm[is2_years-years[0]] = is2_sm_asym
x=years[:]
y=sit_warren[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_warren, intercept, r_value, p_value_warren, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99: ',np.round(trend_warren*100,2),' cm/yr (p=',np.round(p_value_warren,2),')')            

## Warren effsd
years = np.arange(2003,2021,1)
sit_w99eff = np.full(len(years), np.nan)
sit_w99eff[icesat_sm.year.values-years[0]] = is_sm_west.sit_w99_sig.values-is_sm_east.sit_w99_sig.values
sit_w99eff[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_w99.values-cs2_sm_east.sit_w99.values
# sit_w99eff[is2_years-years[0]] = is2_sm_e_asym
x=years[:]
y=sit_w99eff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_w99eff, intercept, r_value, p_value_w99eff, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99 Effsd: ',np.round(trend_w99eff*100,2),' cm/yr (p=',np.round(p_value_w99eff,2),')')        

## Warren Petty
years = np.arange(2003,2021,1)
sit_w99petty = np.full(len(years), np.nan)
sit_w99petty[icesat_sm.year.values-years[0]] = is_sm_west.sit_w99_piece.values-is_sm_east.sit_w99_piece.values
sit_w99petty[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_w99.values-cs2_sm_east.sit_w99.values
# sit_w99petty[is2_years-years[0]] = is2_sm_pet_asym
x=years[:]
y=sit_w99petty[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_w99petty, intercept, r_value, p_value_w99petty, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99 Petty: ',np.round(trend_w99petty*100,2),' cm/yr (p=',np.round(p_value_w99petty,2),')')        


## PMW
years = np.arange(2003,2021,1)
sit_pmw = np.full(len(years), np.nan)
sit_pmw[icesat_sm.year.values-years[0]] = is_sm_west.sit_pmw.values-is_sm_east.sit_pmw.values
sit_pmw[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_pmw.values-cs2_sm_east.sit_pmw.values
# sit_pmw[is2_years-years[0]] = is2_pmw_asym
x=years[:]
y=sit_pmw[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmw, intercept_pmw, r_value, p_value_pmw, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW: ',np.round(trend_pmw*100,2),' cm/yr (p=',np.round(p_value_pmw,2),')')            

## PMW effsd
years = np.arange(2003,2021,1)
sit_pmweff = np.full(len(years), np.nan)
sit_pmweff[icesat_sm.year.values-years[0]] = is_sm_west.sit_pmw_sig.values-is_sm_east.sit_pmw_sig.values
sit_pmweff[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_pmw.values-cs2_sm_east.sit_pmw.values
# sit_pmweff[is2_years-years[0]] = is2_pmw_e_asym
x=years[:]
y=sit_pmweff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmweff, intercept, r_value, p_value_pmweff, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW Effsd: ',np.round(trend_pmweff*100,2),' cm/yr (p=',np.round(p_value_pmweff,2),')')        

# PMW Petty
years = np.arange(2003,2021,1)
sit_pmwpetty = np.full(len(years), np.nan)
sit_pmwpetty[icesat_sm.year.values-years[0]] = is_sm_west.sit_pmw_piece.values-is_sm_east.sit_pmw_piece.values
sit_pmwpetty[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_pmw.values-cs2_sm_east.sit_pmw.values
# sit_pmwpetty[is2_years-years[0]] = is2_pmw_pet_asym
x=years[:]
y=sit_pmwpetty[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmwpetty, intercept_pmwpetty, r_value, p_value_pmwpetty, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW Petty: ',np.round(trend_pmwpetty*100,2),' cm/yr (p=',np.round(p_value_pmwpetty,2),')')        


## Stroeve
years = np.arange(2003,2021,1)
sit_stroe = np.full(len(years), np.nan)
sit_stroe[icesat_sm.year.values-years[0]] = is_sm_west.sit_sm.values-is_sm_east.sit_sm.values
sit_stroe[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_sm.values-cs2_sm_east.sit_sm.values
x=years[:]
y=sit_stroe[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend, intercept, r_value, p_value_stroe, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve: ',np.round(trend_stroe*100,2),' cm/yr (p=',np.round(p_value_stroe,2),')')            

## Stroeve effsd
years = np.arange(2003,2021,1)
sit_stroeeff = np.full(len(years), np.nan)
sit_stroeeff[icesat_sm.year.values-years[0]] = is_sm_west.sit_sm_sig.values-is_sm_east.sit_sm_sig.values
sit_stroeeff[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_sm.values-cs2_sm_east.sit_sm.values
x=years[:]
y=sit_stroeeff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroeeff, intercept, r_value, p_value_stroeeff, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve Effsd: ',np.round(trend_stroeeff*100,2),' cm/yr (p=',np.round(p_value_stroeeff,2),')')            

## Stroeve Petty
years = np.arange(2003,2021,1)
sit_stroepett = np.full(len(years), np.nan)
sit_stroepett[icesat_sm.year.values-years[0]] = is_sm_west.sit_sm_piece.values-is_sm_east.sit_sm_piece.values
sit_stroepett[cs2_sm.year.values-years[0]] = cs2_sm_west.sit_sm.values-cs2_sm_east.sit_sm.values
x=years[:]
y=sit_stroepett[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroepett, intercept, r_value, p_value_stroepett, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve Petty: ',np.round(trend_stroepett*100,2),' cm/yr (p=',np.round(p_value_stroepett,2),')')  


