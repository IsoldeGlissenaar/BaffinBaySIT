# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 11:34:09 2021

@author: IsoldeGlissenaar

Time series SIT for all satellites and processing methods. 
Plots SIT time series envelopes for all processing methods and satellites.
Prints SIT trends over ICESat & CryoSat-2 LARM for all combinations of processing
methods.
"""

import numpy as np
import matplotlib.pyplot as plt
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
""" Mean sea ice thickness methods """

#ICESat-1
eastwest25_2d   = np.transpose(np.tile(eastwest25, (7,1)))
is_warren_mean  = icesat_warren.where(eastwest25_2d>0).mean('n')
is_sm_mean      = icesat_sm.where(eastwest25_2d>0).mean('n')

all_methods = np.array([is_warren_mean.sit_w99.values, is_warren_mean.sit_w99_sig.values, is_warren_mean.sit_w99_piece.values,
                        is_warren_mean.sit_pmw.values, is_warren_mean.sit_pmw_sig.values, is_warren_mean.sit_pmw_piece.values,
                        is_warren_mean.sit_sm.values, is_warren_mean.sit_sm_sig.values, is_warren_mean.sit_sm_piece.values,
                        is_sm_mean.sit_w99.values, is_sm_mean.sit_w99_sig.values, is_sm_mean.sit_w99_piece.values,
                        is_sm_mean.sit_pmw.values, is_sm_mean.sit_pmw_sig.values, is_sm_mean.sit_pmw_piece.values,
                        is_sm_mean.sit_sm.values, is_sm_mean.sit_sm_sig.values, is_sm_mean.sit_sm_piece.values])
min_methods = np.nanmin(all_methods, axis=0)
max_methods = np.nanmax(all_methods, axis=0)

# Uncertainty
all_methods_u = np.array([is_warren_mean.sit_w99_uncertainty.values, is_warren_mean.sit_w99_sig_uncertainty.values, 
                          is_warren_mean.sit_w99_piece_uncertainty.values, is_warren_mean.sit_pmw_uncertainty.values, 
                          is_warren_mean.sit_pmw_sig_uncertainty.values, is_warren_mean.sit_pmw_piece_uncertainty.values,
                          is_warren_mean.sit_sm_uncertainty.values, is_warren_mean.sit_sm_sig_uncertainty.values, 
                          is_warren_mean.sit_sm_piece_uncertainty.values, is_sm_mean.sit_w99_uncertainty.values, 
                          is_sm_mean.sit_w99_sig_uncertainty.values, is_sm_mean.sit_w99_piece_uncertainty.values,
                          is_sm_mean.sit_pmw_uncertainty.values, is_sm_mean.sit_pmw_sig_uncertainty.values, 
                          is_sm_mean.sit_pmw_piece_uncertainty.values, is_sm_mean.sit_sm_uncertainty.values, 
                          is_sm_mean.sit_sm_sig_uncertainty.values, is_sm_mean.sit_sm_piece_uncertainty.values])

all_methods_max = all_methods+all_methods_u
max_methods_u = np.nanmax(all_methods_max, axis=0)
all_methods_min = all_methods-all_methods_u
min_methods_u = np.nanmin(all_methods_min, axis=0)


#%%
#CS2
eastwest12_2d   = np.transpose(np.tile(eastwest12, (10,1)))
cs2_warren_mean  = cs2_warren.where(eastwest12_2d>0).mean('n')
cs2_sm_mean      = cs2_sm.where(eastwest12_2d>0).mean('n')


all_methods_cs2 = np.array([cs2_warren_mean.sit_w99.values,
                            cs2_warren_mean.sit_pmw.values, 
                            cs2_warren_mean.sit_sm.values,
                            cs2_sm_mean.sit_w99.values,
                            cs2_sm_mean.sit_pmw.values,
                            cs2_sm_mean.sit_sm.values])
min_methods_cs2 = np.nanmin(all_methods_cs2, axis=0)
max_methods_cs2 = np.nanmax(all_methods_cs2, axis=0)

# Uncertainty
all_methods_cs2_u = np.array([cs2_warren_mean.sit_w99_uncertainty.values,
                              cs2_warren_mean.sit_pmw_uncertainty.values, 
                              cs2_warren_mean.sit_sm_uncertainty.values,
                              cs2_sm_mean.sit_w99_uncertainty.values,
                              cs2_sm_mean.sit_pmw_uncertainty.values,
                              cs2_sm_mean.sit_sm_uncertainty.values])

all_methods_cs2_max = all_methods_cs2+all_methods_cs2_u
max_methods_cs2_u = np.nanmax(all_methods_cs2_max, axis=0)
all_methods_cs2_min = all_methods_cs2-all_methods_cs2_u
min_methods_cs2_u = np.nanmin(all_methods_cs2_min, axis=0)

# CCI
eastwest12_2d   = np.transpose(np.tile(eastwest12, (10,1)))
cs2_cci_warren_mean  = cs2_cci_warren.where(eastwest12_2d>0).mean('n')
cs2_cci_sm_mean      = cs2_cci_sm.where(eastwest12_2d>0).mean('n')


all_methods_cs2_cci = np.array([cs2_cci_warren_mean.sit_w99.values,
                                cs2_cci_warren_mean.sit_pmw.values, 
                                cs2_cci_warren_mean.sit_sm.values,
                                cs2_cci_sm_mean.sit_w99.values,
                                cs2_cci_sm_mean.sit_pmw.values,
                                cs2_cci_sm_mean.sit_sm.values])
min_methods_cs2_cci = np.nanmin(all_methods_cs2_cci, axis=0)
max_methods_cs2_cci = np.nanmax(all_methods_cs2_cci, axis=0)



#%%
#Envisat
eastwest25_2d   = np.transpose(np.tile(eastwest25, (10,1)))
envi_warren_mean  = envi_warren.where(eastwest25_2d>0).mean('n')
envi_sm_mean      = envi_sm.where(eastwest25_2d>0).mean('n')

all_methods_envi = np.array([envi_warren_mean.sit_w99.values,
                             envi_warren_mean.sit_pmw.values,
                             envi_warren_mean.sit_sm.values,
                             envi_sm_mean.sit_w99.values,
                             envi_sm_mean.sit_pmw.values,
                             envi_sm_mean.sit_sm.values])
min_methods_envi = np.nanmin(all_methods_envi, axis=0)
max_methods_envi = np.nanmax(all_methods_envi, axis=0)

#Uncertainty
all_methods_envi_u = np.array([envi_warren_mean.sit_w99_uncertainty.values,
                               envi_warren_mean.sit_pmw_uncertainty.values, 
                               envi_warren_mean.sit_sm_uncertainty.values,
                               envi_sm_mean.sit_w99_uncertainty.values,
                               envi_sm_mean.sit_pmw_uncertainty.values,
                               envi_sm_mean.sit_sm_uncertainty.values])

all_methods_envi_max = all_methods_envi+all_methods_envi_u
max_methods_envi_u = np.nanmax(all_methods_envi_max, axis=0)
all_methods_envi_min = all_methods_envi-all_methods_envi_u
min_methods_envi_u = np.nanmin(all_methods_envi_min, axis=0)

#%%
#ICESat-2
eastwest12_2d   = np.transpose(np.tile(eastwest12, (2,1)))
is2_warren_mean  = icesat2_warren.where(eastwest12_2d>0).mean('n')
is2_sm_mean      = icesat2_sm.where(eastwest12_2d>0).mean('n')

all_methods_is2 = np.array([is2_warren_mean.sit_w99.values, is2_warren_mean.sit_w99_sig.values, is2_warren_mean.sit_w99_piece.values,
                            is2_warren_mean.sit_pmw.values, is2_warren_mean.sit_pmw_sig.values, is2_warren_mean.sit_pmw_piece.values,
                            is2_warren_mean.sit_sm.values, is2_warren_mean.sit_sm_sig.values, is2_warren_mean.sit_sm_piece.values,
                            is2_sm_mean.sit_w99.values, is2_sm_mean.sit_w99_sig.values, is2_sm_mean.sit_w99_piece.values,
                            is2_sm_mean.sit_pmw.values, is2_sm_mean.sit_pmw_sig.values, is2_sm_mean.sit_pmw_piece.values,
                            is2_sm_mean.sit_sm.values, is2_sm_mean.sit_sm_sig.values, is2_sm_mean.sit_sm_piece.values])
min_methods_is2 = np.nanmin(all_methods_is2, axis=0)
max_methods_is2 = np.nanmax(all_methods_is2, axis=0)

# Uncertainty
all_methods_is2_u = np.array([is2_warren_mean.sit_w99_uncertainty.values, is2_warren_mean.sit_w99_sig_uncertainty.values, 
                              is2_warren_mean.sit_w99_piece_uncertainty.values, is2_warren_mean.sit_pmw_uncertainty.values, 
                              is2_warren_mean.sit_pmw_sig_uncertainty.values, is2_warren_mean.sit_pmw_piece_uncertainty.values,
                              is2_warren_mean.sit_sm_uncertainty.values, is2_warren_mean.sit_sm_sig_uncertainty.values, 
                              is2_warren_mean.sit_sm_piece_uncertainty.values, is2_sm_mean.sit_w99_uncertainty.values, 
                              is2_sm_mean.sit_w99_sig_uncertainty.values, is2_sm_mean.sit_w99_piece_uncertainty.values,
                              is2_sm_mean.sit_pmw_uncertainty.values, is2_sm_mean.sit_pmw_sig_uncertainty.values, 
                              is2_sm_mean.sit_pmw_piece_uncertainty.values, is2_sm_mean.sit_sm_uncertainty.values, 
                              is2_sm_mean.sit_sm_sig_uncertainty.values, is2_sm_mean.sit_sm_piece_uncertainty.values])

all_methods_is2_max = all_methods_is2+all_methods_is2_u
max_methods_is2_u = np.nanmax(all_methods_is2_max, axis=0)
all_methods_is2_min = all_methods_is2-all_methods_is2_u
min_methods_is2_u = np.nanmin(all_methods_is2_min, axis=0)

#%%
#CIS
cis_bb_mean = np.nanmean(cis_sit[:,(eastwest_cis>0)], axis=1)

#%%
#Plot
fig, ax = plt.subplots(dpi=200)

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

plt.plot(cis_years, cis_bb_mean, c='C4', linestyle='--', label='CIS')

plt.grid(':', linewidth=0.5)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylabel('March mean SIT [m]')
plt.xticks(np.arange(2000,2022,2), rotation=45)
plt.ylim(bottom=0)
plt.xlim(left=2000)
ax.set_xticks(np.arange(2001,2021,2), minor=True)
# plt.savefig('../../figures/fig6.png', bbox_inches='tight', dpi=300)

#%%
from scipy import stats
## Trends diff processing techniques

print('-------')
print('Warren snow density')
## Warren
years = np.arange(2003,2021,1)
sit_warren = np.full(len(years), np.nan)
sit_warren[icesat_warren.year.values-years[0]] = is_warren_mean.sit_w99.values
sit_warren[cs2_warren.year.values-years[0]] = cs2_warren_mean.sit_w99.values
# sit_warren[is2_years-years[0]] = is2_warren_bb_mean
x=years[:]
y=sit_warren[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_warren, intercept, r_value, p_value_warren, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99: ',np.round(trend_warren*100,2),' cm/yr (p=',np.round(p_value_warren,2),')')            

## Warren effsd
years = np.arange(2003,2021,1)
sit_w99eff = np.full(len(years), np.nan)
sit_w99eff[icesat_warren.year.values-years[0]] = is_warren_mean.sit_w99_sig.values
sit_w99eff[cs2_warren.year.values-years[0]] =  cs2_warren_mean.sit_w99.values
# sit_w99eff[is2_years-years[0]] = is2_warren_eff_bb_mean
x=years[:]
y=sit_w99eff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_w99eff, intercept, r_value, p_value_w99eff, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99 Effsd: ',np.round(trend_w99eff*100,2),' cm/yr (p=',np.round(p_value_w99eff,2),')')        

## Warren Petty
years = np.arange(2003,2021,1)
sit_w99petty = np.full(len(years), np.nan)
sit_w99petty[icesat_warren.year.values-years[0]] = is_warren_mean.sit_w99_piece.values
sit_w99petty[cs2_warren.year.values-years[0]] = cs2_warren_mean.sit_w99.values
# sit_w99petty[is2_years-years[0]] = is2_warren_petty_bb_mean
x=years[:]
y=sit_w99petty[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_w99petty, intercept, r_value, p_value_w99petty, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99 Petty: ',np.round(trend_w99petty*100,2),' cm/yr (p=',np.round(p_value_w99petty,2),')')        


## PMW
years = np.arange(2003,2021,1)
sit_pmw = np.full(len(years), np.nan)
sit_pmw[icesat_warren.year.values-years[0]] = is_warren_mean.sit_pmw.values
sit_pmw[cs2_warren.year.values-years[0]] = cs2_warren_mean.sit_pmw.values
# sit_pmw[is2_years-years[0]] = is2_pmw_bb_mean
x=years[:]
y=sit_pmw[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmw, intercept, r_value, p_value_pmw, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW: ',np.round(trend_pmw*100,2),' cm/yr (p=',np.round(p_value_pmw,2),')')            

## PMW effsd
years = np.arange(2003,2021,1)
sit_pmweff = np.full(len(years), np.nan)
sit_pmweff[icesat_warren.year.values-years[0]] = is_warren_mean.sit_pmw_sig.values
sit_pmweff[cs2_warren.year.values-years[0]] = cs2_warren_mean.sit_pmw.values
# sit_pmweff[is2_years-years[0]] = is2_pmw_eff_bb_mean
x=years[:]
y=sit_pmweff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmweff, intercept, r_value, p_value_pmweff, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW Effsd: ',np.round(trend_pmweff*100,2),' cm/yr (p=',np.round(p_value_pmweff,2),')')        

## PMW Petty
years = np.arange(2003,2021,1)
sit_pmwpetty = np.full(len(years), np.nan)
sit_pmwpetty[icesat_warren.year.values-years[0]] = is_warren_mean.sit_pmw_piece.values
sit_pmwpetty[cs2_warren.year.values-years[0]] = cs2_warren_mean.sit_pmw.values
# sit_pmwpetty[is2_years-years[0]] = is2_pmw_petty_bb_mean
x=years[:]
y=sit_pmwpetty[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmwpetty, intercept_pmwpetty, r_value, p_value_pmwpetty, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW Petty: ',np.round(trend_pmwpetty*100,2),' cm/yr (p=',np.round(p_value_pmwpetty,2),')')        


## Stroeve
years = np.arange(2003,2021,1)
sit_stroe = np.full(len(years), np.nan)
sit_stroe[icesat_warren.year.values-years[0]] = is_warren_mean.sit_sm.values
sit_stroe[cs2_warren.year.values-years[0]] = cs2_warren_mean.sit_sm.values
x=years[:]
y=sit_stroe[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroe, intercept_stroe, r_value, p_value_stroe, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve: ',np.round(trend_stroe*100,2),' cm/yr (p=',np.round(p_value_stroe,2),')')            

## Stroeve effsd
years = np.arange(2003,2021,1)
sit_stroeeff = np.full(len(years), np.nan)
sit_stroeeff[icesat_warren.year.values-years[0]] = is_warren_mean.sit_sm_sig.values
sit_stroeeff[cs2_warren.year.values-years[0]] = cs2_warren_mean.sit_sm.values
x=years[:]
y=sit_stroeeff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroeeff, intercept, r_value, p_value_stroeeff, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve Effsd: ',np.round(trend_stroeeff*100,2),' cm/yr (p=',np.round(p_value_stroeeff,2),')')            

## Stroeve Petty
years = np.arange(2003,2021,1)
sit_stroepett = np.full(len(years), np.nan)
sit_stroepett[icesat_warren.year.values-years[0]] = is_warren_mean.sit_sm_piece.values
sit_stroepett[cs2_warren.year.values-years[0]] = cs2_warren_mean.sit_sm.values
x=years[:]
y=sit_stroepett[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroepett, intercept, r_value, p_value_stroepett, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve Petty: ',np.round(trend_stroepett*100,2),' cm/yr (p=',np.round(p_value_stroepett,2),')')            

#%%
## Trends diff processing techniques

print('--------')
print('SnowModel snowdensity')
## Warren
years = np.arange(2003,2021,1)
sit_warren = np.full(len(years), np.nan)
sit_warren[icesat_warren.year.values-years[0]] = is_sm_mean.sit_w99.values
sit_warren[cs2_warren.year.values-years[0]] = cs2_sm_mean.sit_w99.values
# sit_warren[is2_years-years[0]] = is2_warren_bb_mean
x=years[:]
y=sit_warren[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_warren, intercept, r_value, p_value_warren, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99: ',np.round(trend_warren*100,2),' cm/yr (p=',np.round(p_value_warren,2),')')            

## Warren effsd
years = np.arange(2003,2021,1)
sit_w99eff = np.full(len(years), np.nan)
sit_w99eff[icesat_warren.year.values-years[0]] = is_sm_mean.sit_w99_sig.values
sit_w99eff[cs2_warren.year.values-years[0]] =  cs2_sm_mean.sit_w99.values
# sit_w99eff[is2_years-years[0]] = is2_warren_eff_bb_mean
x=years[:]
y=sit_w99eff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_w99eff, intercept_w99eff, r_value, p_value_w99eff, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99 Effsd: ',np.round(trend_w99eff*100,2),' cm/yr (p=',np.round(p_value_w99eff,2),')')        

## Warren Petty
years = np.arange(2003,2021,1)
sit_w99petty = np.full(len(years), np.nan)
sit_w99petty[icesat_warren.year.values-years[0]] = is_sm_mean.sit_w99_piece.values
sit_w99petty[cs2_warren.year.values-years[0]] = cs2_sm_mean.sit_w99.values
# sit_w99petty[is2_years-years[0]] = is2_warren_petty_bb_mean
x=years[:]
y=sit_w99petty[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_w99petty, intercept, r_value, p_value_w99petty, std_err = stats.linregress(x[mask], y[mask])
print('Trend W99 Petty: ',np.round(trend_w99petty*100,2),' cm/yr (p=',np.round(p_value_w99petty,2),')')        


## PMW
years = np.arange(2003,2021,1)
sit_pmw = np.full(len(years), np.nan)
sit_pmw[icesat_warren.year.values-years[0]] = is_sm_mean.sit_pmw.values
sit_pmw[cs2_warren.year.values-years[0]] = cs2_sm_mean.sit_pmw.values
# sit_pmw[is2_years-years[0]] = is2_pmw_bb_mean
x=years[:]
y=sit_pmw[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmw, intercept, r_value, p_value_pmw, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW: ',np.round(trend_pmw*100,2),' cm/yr (p=',np.round(p_value_pmw,2),')')            

## PMW effsd
years = np.arange(2003,2021,1)
sit_pmweff = np.full(len(years), np.nan)
sit_pmweff[icesat_warren.year.values-years[0]] = is_sm_mean.sit_pmw_sig.values
sit_pmweff[cs2_warren.year.values-years[0]] = cs2_sm_mean.sit_pmw.values
# sit_pmweff[is2_years-years[0]] = is2_pmw_eff_bb_mean
x=years[:]
y=sit_pmweff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_pmweff, intercept, r_value, p_value_pmweff, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW Effsd: ',np.round(trend_pmweff*100,2),' cm/yr (p=',np.round(p_value_pmweff,2),')')        

## PMW Petty
years = np.arange(2003,2021,1)
sit_pmwpetty = np.full(len(years), np.nan)
sit_pmwpetty[icesat_warren.year.values-years[0]] = is_sm_mean.sit_pmw_piece.values
sit_pmwpetty[cs2_warren.year.values-years[0]] = cs2_sm_mean.sit_pmw.values
# sit_pmwpetty[is2_years-years[0]] = is2_pmw_petty_bb_mean
x=years[:]
y=sit_pmwpetty[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend, intercept, r_value, p_value_pmwpetty, std_err = stats.linregress(x[mask], y[mask])
print('Trend PMW Petty: ',np.round(trend_pmwpetty*100,2),' cm/yr (p=',np.round(p_value_pmwpetty,2),')')        


## Stroeve
years = np.arange(2003,2021,1)
sit_stroe = np.full(len(years), np.nan)
sit_stroe[icesat_warren.year.values-years[0]] = is_sm_mean.sit_sm.values
sit_stroe[cs2_warren.year.values-years[0]] = cs2_sm_mean.sit_sm.values
x=years[:]
y=sit_stroe[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroe, intercept_stroe, r_value, p_value_stroe, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve: ',np.round(trend_stroe*100,2),' cm/yr (p=',np.round(p_value_stroe,2),')')            

## Stroeve effsd
years = np.arange(2003,2021,1)
sit_stroeeff = np.full(len(years), np.nan)
sit_stroeeff[icesat_warren.year.values-years[0]] = is_sm_mean.sit_sm_sig.values
sit_stroeeff[cs2_warren.year.values-years[0]] = cs2_sm_mean.sit_sm.values
x=years[:]
y=sit_stroeeff[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroeeff, intercept, r_value, p_value_stroeeff, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve Effsd: ',np.round(trend_stroeeff*100,2),' cm/yr (p=',np.round(p_value_stroeeff,2),')')            

## Stroeve Petty
years = np.arange(2003,2021,1)
sit_stroepett = np.full(len(years), np.nan)
sit_stroepett[icesat_warren.year.values-years[0]] = is_sm_mean.sit_sm_piece.values
sit_stroepett[cs2_warren.year.values-years[0]] = cs2_sm_mean.sit_sm.values
x=years[:]
y=sit_stroepett[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_stroepett, intercept, r_value, p_value_stroepett, std_err = stats.linregress(x[mask], y[mask])
print('Trend Stroeve Petty: ',np.round(trend_stroepett*100,2),' cm/yr (p=',np.round(p_value_stroepett,2),')')      
#%%
## CIS charts
x=cis_years[:]
y=cis_bb_mean[:]
mask = ~np.isnan(x) & ~np.isnan(y)
trend_cis, intercept, r_value, p_value_cis, std_err = stats.linregress(x[mask], y[mask])
print('-----')
print('Trend CIS: ',np.round(trend_cis*100,2),' cm/yr (p=',np.round(p_value_cis,2),')')            

#%% 
# Mean uncertainty envelope ??????
middle_is = 0.5*(max_methods-min_methods)+min_methods
uncer_envelope_is = 0.5*(max_methods-min_methods)
percent_envelope_is = (uncer_envelope_is)/middle_is*100

# max_methods_cs2[:7] = np.nanmax(np.array([max_methods_cs2[:7],max_methods_cs2_cci[:7]]),axis=0)
# min_methods_cs2[:7] = np.nanmin(np.array([min_methods_cs2[:7],min_methods_cs2_cci[:7]]),axis=0)

middle_envi = 0.5*(max_methods_envi-min_methods_envi)+min_methods_envi
uncer_envelope_envi = 0.5*(max_methods_envi-min_methods_envi)
percent_envelope_envi = (uncer_envelope_envi)/middle_envi*100

middle_cs = 0.5*(max_methods_cs2-min_methods_cs2)+min_methods_cs2
uncer_envelope_cs = 0.5*(max_methods_cs2-min_methods_cs2)
percent_envelope_cs = (uncer_envelope_cs)/middle_cs*100

middle_is2 = 0.5*(max_methods_is2-min_methods_is2)+min_methods_is2
uncer_envelope_is2 = 0.5*(max_methods_is2-min_methods_is2)
percent_envelope_is2 = (uncer_envelope_is2)/middle_is2*100

all_perenv = []
all_perenv = np.append(all_perenv, percent_envelope_envi)
all_perenv = np.append(all_perenv, percent_envelope_is)
all_perenv = np.append(all_perenv, percent_envelope_cs)
all_perenv = np.append(all_perenv, percent_envelope_is2)
print('-----')
print(f'mean uncertainty envelope around the March mean sea ice thickness: {np.round(np.nanmean(all_perenv),1)}%')

#%%
# Mean uncertainty envelope around long-term trends
min_trend = np.nanmin([trend_warren, trend_w99eff, trend_w99petty, 
                      trend_pmw, trend_pmweff, trend_pmwpetty,
                      trend_stroe, trend_stroeeff, trend_stroepett])
max_trend = np.nanmax([trend_warren, trend_w99eff, trend_w99petty, 
                      trend_pmw, trend_pmweff, trend_pmwpetty,
                      trend_stroe, trend_stroeeff, trend_stroepett])

middle_trend = (max_trend+min_trend)/2
uncer_envelope_trend = 0.5*(max_trend-min_trend)
percent_envelope_trend = (uncer_envelope_trend)/middle_trend*100

print('-----')
print(f'mean uncertainty envelope around the long-term trend: {np.round(percent_envelope_trend,1)}%')

