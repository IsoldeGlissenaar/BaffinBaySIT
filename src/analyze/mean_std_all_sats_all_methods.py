# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:37:11 2021

@author: IsoldeGlissenaar

Plot mean and standard deviation in Sea Ice Thickness between different
processing methods for all satellites. Print mean SIT and trends.
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import scipy.stats
import xarray as xr
import cartopy.feature as cfeature
from shapely.geometry import Point # Point class
from shapely.geometry import shape # shape() is a function to convert geo objects through the interface
import shapefile
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')


#%%
is_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is/SIT_ICESat_March_2003_2009_SnowModelsnowdens.nc')
is_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is/SIT_ICESat_March_2003_2009_Warrensnowdens.nc')

envi_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/envi/SIT_Envisat_March_2003_2012_SnowModelsnowdens.nc')
envi_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/envi/SIT_Envisat_March_2003_2012_Warrensnowdens.nc')

cs2_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2/SIT_CryoSat2_March_2011_2020_SnowModelsnowdens.nc')
cs2_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2/SIT_CryoSat2_March_2011_2020_Warrensnowdens.nc')

cs2_cci_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2_cci/SIT_CryoSat2_March_2011_2017_SnowModelsnowdens.nc')
cs2_cci_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2_cci/SIT_CryoSat2_March_2011_2017_Warrensnowdens.nc')

is2_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is2/SIT_ICESat2_March_2019_2020_SnowModelsnowdens.nc')
is2_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is2/SIT_ICESat2_March_2019_2020_Warrensnowdens.nc')
 

#%%
'''ICESat'''

icesat = np.array([is_warren.sit_w99.values, is_warren.sit_w99_sig.values, is_warren.sit_w99_piece.values,
                   is_warren.sit_pmw.values, is_warren.sit_pmw_sig.values, is_warren.sit_pmw_piece.values,
                   is_warren.sit_sm.values, is_warren.sit_sm_sig.values, is_warren.sit_sm_piece.values,
                   is_sm.sit_w99.values, is_sm.sit_w99_sig.values, is_sm.sit_w99_piece.values,
                   is_sm.sit_pmw.values, is_sm.sit_pmw_sig.values, is_sm.sit_pmw_piece.values,
                   is_sm.sit_sm.values, is_sm.sit_sm_sig.values, is_sm.sit_sm_piece.values])

icesat_uncer = np.array([is_warren.sit_w99_uncertainty.values, is_warren.sit_w99_sig_uncertainty.values, is_warren.sit_w99_piece_uncertainty.values,
                         is_warren.sit_pmw_uncertainty.values, is_warren.sit_pmw_sig_uncertainty.values, is_warren.sit_pmw_piece_uncertainty.values,
                         is_warren.sit_sm_uncertainty.values, is_warren.sit_sm_sig_uncertainty.values, is_warren.sit_sm_piece_uncertainty.values,
                         is_sm.sit_w99_uncertainty.values, is_sm.sit_w99_sig_uncertainty.values, is_sm.sit_w99_piece_uncertainty.values,
                         is_sm.sit_pmw_uncertainty.values, is_sm.sit_pmw_sig_uncertainty.values, is_sm.sit_pmw_piece_uncertainty.values,
                         is_sm.sit_sm_uncertainty.values, is_sm.sit_sm_sig_uncertainty.values, is_sm.sit_sm_piece_uncertainty.values])

is_yr_mean = np.nanmean(icesat,axis=2)
is_sats_mean = np.nanmean(is_yr_mean,axis=0)

is_yr_mean_uncer = np.nanmean(icesat_uncer,axis=2)
is_sats_mean_uncer = np.nanmean(is_yr_mean_uncer,axis=0)

is_sats_std = scipy.stats.sem(is_yr_mean,axis=0)

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m', linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon, is_warren.lat, c=is_sats_mean,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='sea ice thickness [m]')
plt.title('March mean ICESat (2003-2009) ',fontsize=8)
# plt.savefig('../../figures/figure5b/icesat_mean.png',dpi=300)

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m', linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon, is_warren.lat, c=is_sats_std,cmap='hot_r',vmin=0,vmax=0.15,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.title('March stdev ICESat ',fontsize=8)
# plt.savefig('../../figures/figure5b/icesat_std.png',dpi=300)

#%%
'''CryoSat-2''' 

cryosat = np.array([cs2_warren.sit_w99.values, cs2_warren.sit_pmw.values, cs2_warren.sit_sm.values, 
                    cs2_sm.sit_w99.values, cs2_sm.sit_pmw.values, cs2_sm.sit_sm.values, ])

cryosat_uncer = np.array([cs2_warren.sit_w99_uncertainty.values, cs2_warren.sit_pmw_uncertainty.values, cs2_warren.sit_sm_uncertainty.values, 
                          cs2_sm.sit_w99_uncertainty.values, cs2_sm.sit_pmw_uncertainty.values, cs2_sm.sit_sm_uncertainty.values])

cs2_yr_mean = np.nanmean(cryosat,axis=2)
cs2_sats_mean = np.nanmean(cs2_yr_mean,axis=0)

cs2_yr_mean_uncer = np.nanmean(cryosat_uncer,axis=2)
cs2_sats_mean_uncer = np.nanmean(cs2_yr_mean_uncer,axis=0)

cs2_sats_std = scipy.stats.sem(cs2_yr_mean,axis=0)

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree())
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(cs2_warren.lon, cs2_warren.lat, c=cs2_sats_mean, cmap='Spectral_r',vmin=0,vmax=2.5,s=1.5,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='sea ice thickness [m]')
plt.title('March mean SIT CS2 (2011-2020) ',fontsize=8)
# plt.savefig('../../figures/figure5b/cryosat_mean.png',dpi=300)

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree())
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(cs2_warren.lon, cs2_warren.lat, c=cs2_sats_std, cmap='hot_r',vmin=0,vmax=0.15,s=1.5,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.title('March std SIT CS2 (2011-2020) ',fontsize=8)
# plt.savefig('../../figures/figure5b/cryosat_std.png',dpi=300)

#%%
'''ENVISAT'''  

envisat = np.array([envi_warren.sit_w99.values, envi_warren.sit_pmw.values, envi_warren.sit_sm.values, 
                    envi_sm.sit_w99.values, envi_sm.sit_pmw.values, envi_sm.sit_sm.values, ])

envisat_uncer = np.array([envi_warren.sit_w99_uncertainty.values, envi_warren.sit_pmw_uncertainty.values, envi_warren.sit_sm_uncertainty.values, 
                          envi_sm.sit_w99_uncertainty.values, envi_sm.sit_pmw_uncertainty.values, envi_sm.sit_sm_uncertainty.values])

envi_yr_mean = np.nanmean(envisat,axis=2)
envi_sats_mean = np.nanmean(envi_yr_mean,axis=0)

envi_yr_mean_uncer = np.nanmean(envisat_uncer,axis=2)
envi_sats_mean_uncer = np.nanmean(envi_yr_mean_uncer,axis=0)

envi_sats_std = scipy.stats.sem(envi_yr_mean,axis=0)

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree())
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(envi_warren.lon, envi_warren.lat, c=envi_sats_mean, cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='sea ice thickness [m]')
plt.title('March mean SIT Envisat (2003-2012) ',fontsize=8)
# plt.savefig('../../figures/figure5b/envisat_mean.png',dpi=300)

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
# ax.set_extent([-90,-40,58,85],crs=ccrs.PlateCarree()) 
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree())
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(envi_warren.lon, envi_warren.lat, c=envi_sats_std, cmap='hot_r',vmin=0,vmax=0.15,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
# ax.add_patch(poly)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.title('March std SIT Envisat (2003-2012) ',fontsize=8)
# plt.savefig('../../figures/figure5b/envisat_std.png',dpi=300)

#%%
''' ICESat-2'''

icesat2 = np.array([is2_warren.sit_w99.values, is2_warren.sit_w99_sig.values, is2_warren.sit_w99_piece.values,
                    is2_warren.sit_pmw.values, is2_warren.sit_pmw_sig.values, is2_warren.sit_pmw_piece.values,
                    is2_warren.sit_sm.values, is2_warren.sit_sm_sig.values, is2_warren.sit_sm_piece.values,
                    is2_sm.sit_w99.values, is2_sm.sit_w99_sig.values, is2_sm.sit_w99_piece.values,
                    is2_sm.sit_pmw.values, is2_sm.sit_pmw_sig.values, is2_sm.sit_pmw_piece.values,
                    is2_sm.sit_sm.values, is2_sm.sit_sm_sig.values, is2_sm.sit_sm_piece.values])

icesat2_uncer = np.array([is2_warren.sit_w99_uncertainty.values, is2_warren.sit_w99_sig_uncertainty.values, is2_warren.sit_w99_piece_uncertainty.values,
                          is2_warren.sit_pmw_uncertainty.values, is2_warren.sit_pmw_sig_uncertainty.values, is2_warren.sit_pmw_piece_uncertainty.values,
                          is2_warren.sit_sm_uncertainty.values, is2_warren.sit_sm_sig_uncertainty.values, is2_warren.sit_sm_piece_uncertainty.values,
                          is2_sm.sit_w99_uncertainty.values, is2_sm.sit_w99_sig_uncertainty.values, is2_sm.sit_w99_piece_uncertainty.values,
                          is2_sm.sit_pmw_uncertainty.values, is2_sm.sit_pmw_sig_uncertainty.values, is2_sm.sit_pmw_piece_uncertainty.values,
                          is2_sm.sit_sm_uncertainty.values, is2_sm.sit_sm_sig_uncertainty.values, is2_sm.sit_sm_piece_uncertainty.values])

is2_yr_mean = np.nanmean(icesat2,axis=2)
is2_sats_mean = np.nanmean(is2_yr_mean,axis=0)

is2_yr_mean_uncer = np.nanmean(icesat2_uncer,axis=2)
is2_sats_mean_uncer = np.nanmean(is2_yr_mean_uncer,axis=0)

is2_sats_std = scipy.stats.sem(is2_yr_mean,axis=0)
      

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree())
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is2_warren.lon, is2_warren.lat, c=is2_sats_mean, cmap='Spectral_r',vmin=0,vmax=2.5,s=1.5,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
fig.colorbar(im, ax=ax, label='sea ice thickness [m]', fraction=0.046, pad=0.04)
plt.title('March mean SIT ICESat-2 (2019-2020) ',fontsize=8)
# plt.savefig('../../figures/figure5b/icesat2_mean.png',dpi=300)

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree())
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is2_warren.lon, is2_warren.lat, c=is2_sats_std, cmap='hot_r',vmin=0,vmax=0.15,s=1.5,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
cbar = fig.colorbar(im, ax=ax, label='standard error [m]', fraction=0.046, pad=0.04)
cbar.ax.locator_params(nbins=5)
plt.title('March std SIT ICESat-2 (2019-2020) ',fontsize=8)
# plt.savefig('../../figures/figure5b/icesat2_std.png',dpi=300)


#%%
#Find east-west location gridcells
shp = shapefile.Reader(docs+'/QGIS/4quad.shp') #open the shapefile
all_shapes = shp.shapes() # get all the polygons
all_records = shp.records()     

#Get east/west location gridcells on 25 km grid (ICESat & Envisat)
len_f = len(is_warren.lat)
eastwest25 = np.zeros(len_f)     
for i in range(len_f):
    pt = (is_warren.lon.values[i],is_warren.lat.values[i])  
    for k in range (len(all_shapes)):             
        boundary = all_shapes[k]
        if Point(pt).within(shape(boundary)): # make a point and see if it's in the polygon
            eastwest25[i] = all_records[k][0]

#Get east/west location gridcells on 12 km grid (ICESat-2 & CryoSat-2)
len_f = len(cs2_warren.lat)
eastwest12 = np.zeros(len_f)     
for i in range(len_f):
    pt = (cs2_warren.lon.values[i],cs2_warren.lat.values[i])  
    for k in range (len(all_shapes)):             
        boundary = all_shapes[k]
        if Point(pt).within(shape(boundary)): # make a point and see if it's in the polygon
            eastwest12[i] = all_records[k][0]
            
#%%    
is_mean_bb = np.nanmean(is_sats_mean[(eastwest25>0)])
envi_mean_bb = np.nanmean(envi_sats_mean[(eastwest25>0)])
cs2_mean_bb = np.nanmean(cs2_sats_mean[(eastwest12>0)])
is2_mean_bb = np.nanmean(is2_sats_mean[(eastwest12>0)])

is_mean_u_bb = np.nanmean(is_sats_mean_uncer[(eastwest25>0)])
envi_mean_u_bb = np.nanmean(envi_sats_mean_uncer[(eastwest25>0)])
cs2_mean_u_bb = np.nanmean(cs2_sats_mean_uncer[(eastwest12>0)])
is2_mean_u_bb = np.nanmean(is2_sats_mean_uncer[(eastwest12>0)])

print(f'March mean (2003-2009) thickness ICESat: {np.round(is_mean_bb,2)}+-{np.round(is_mean_u_bb,2)} m')
print(f'March mean (2003-2012) thickness Envisat: {np.round(envi_mean_bb,2)}+-{np.round(envi_mean_u_bb,2)} m')
print(f'March mean (2011-2020) thickness CryoSat-2: {np.round(cs2_mean_bb,2)}+-{np.round(cs2_mean_u_bb,2)} m')
print(f'March mean (2019-2020) thickness ICESat-2: {np.round(is2_mean_bb,2)}+-{np.round(is2_mean_u_bb,2)} m')
print('----')     

#%%
# Envisat + CS2 mean shorter timeframe

envi_sats_mean = np.nanmean(np.nanmean(envisat[:,:,:7], axis=2),axis=0)
envi_mean_bb = np.nanmean(envi_sats_mean[(eastwest25>0)])
envi_sats_mean_u = np.nanmean(np.nanmean(envisat_uncer[:,:,:7], axis=2),axis=0)
envi_mean_bb_u = np.nanmean(envi_sats_mean_u[(eastwest25>0)])

print(f'March mean (2003-2009) thickness Envisat: {np.round(envi_mean_bb,2)}+-{np.round(envi_mean_bb_u,2)} m')

envi_sats_mean = np.nanmean(np.nanmean(envisat[:,:,8:], axis=2),axis=0)
envi_mean_bb = np.nanmean(envi_sats_mean[(eastwest25>0)])
envi_sats_mean_u = np.nanmean(np.nanmean(envisat_uncer[:,:,8:], axis=2),axis=0)
envi_mean_bb_u = np.nanmean(envi_sats_mean_u[(eastwest25>0)])

print(f'March mean (2011-2012) thickness Envisat: {np.round(envi_mean_bb,2)}+-{np.round(envi_mean_bb_u,2)} m')

cs2_sats_mean = np.nanmean(np.nanmean(cryosat[:,:,:2], axis=2),axis=0)
cs2_mean_bb = np.nanmean(cs2_sats_mean[(eastwest12>0)])
cs2_sats_mean_u = np.nanmean(np.nanmean(cryosat_uncer[:,:,:2], axis=2),axis=0)
cs2_mean_bb_u = np.nanmean(cs2_sats_mean_u[(eastwest12>0)])

print(f'March mean (2011-2012) thickness CryoSat-2: {np.round(cs2_mean_bb,2)}+-{np.round(cs2_mean_bb_u,2)} m')

cs2_sats_mean = np.nanmean(np.nanmean(cryosat[:,:,8:], axis=2),axis=0)
cs2_mean_bb = np.nanmean(cs2_sats_mean[(eastwest12>0)])
cs2_sats_mean_u = np.nanmean(np.nanmean(cryosat_uncer[:,:,8:], axis=2),axis=0)
cs2_mean_bb_u = np.nanmean(cs2_sats_mean_u[(eastwest12>0)])

print(f'March mean (2019-2020) thickness CryoSat-2: {np.round(cs2_mean_bb,2)}+-{np.round(cs2_mean_bb_u,2)} m')

#%%
is_sats_mean_time = np.nanmean(np.nanmean(icesat[:,(eastwest25>0),:], axis=1),axis=0)

cs2_sats_mean_time = np.nanmean(np.nanmean(cryosat[:,(eastwest12>0),:], axis=1),axis=0)

yrs = np.arange(2003,2021)
mean_sit = np.append(is_sats_mean_time, np.nan)
mean_sit = np.append(mean_sit,cs2_sats_mean_time)
mask = ~np.isnan(yrs) & ~np.isnan(mean_sit)

trend, intercept, r_value, p_value, std_err = scipy.stats.linregress(yrs[mask], mean_sit[mask])
print('----------')
print('Trend SIT: ',np.round(trend*100,2),' cm/yr (p=',np.round(p_value,2),')')   

#%%
is_west_mean_time = np.nanmean(np.nanmean(icesat[:,((eastwest25==3)|(eastwest25==4)),:], axis=1),axis=0)
is_east_mean_time = np.nanmean(np.nanmean(icesat[:,((eastwest25==1)|(eastwest25==2)),:], axis=1),axis=0)

cs2_west_mean_time = np.nanmean(np.nanmean(cryosat[:,((eastwest12==3)|(eastwest12==4)),:], axis=1),axis=0)
cs2_east_mean_time = np.nanmean(np.nanmean(cryosat[:,((eastwest12==1)|(eastwest12==2)),:], axis=1),axis=0)

is_asym = is_west_mean_time - is_east_mean_time
cs_asym = cs2_west_mean_time - cs2_east_mean_time

yrs = np.arange(2003,2021)
mean_asym = np.append(is_asym, np.nan)
mean_asym = np.append(mean_asym,cs_asym)
mask = ~np.isnan(yrs) & ~np.isnan(mean_asym)

trend, intercept, r_value, p_value, std_err = scipy.stats.linregress(yrs[mask], mean_asym[mask])
print('----------')
print('Trend asymmetry (west-east): ',np.round(trend*100,2),' cm/yr (p=',np.round(p_value,2),')')   

yrs = np.arange(2003,2021)
mean_asym = np.append(is_west_mean_time, np.nan)
mean_asym = np.append(mean_asym,cs2_west_mean_time)
mask = ~np.isnan(yrs) & ~np.isnan(mean_asym)

trend, intercept, r_value, p_value, std_err = scipy.stats.linregress(yrs[mask], mean_asym[mask])
print('----------')
print('Trend west: ',np.round(trend*100,2),' cm/yr (p=',np.round(p_value,2),')')   


yrs = np.arange(2003,2021)
mean_asym = np.append(is_east_mean_time, np.nan)
mean_asym = np.append(mean_asym,cs2_east_mean_time)
mask = ~np.isnan(yrs) & ~np.isnan(mean_asym)

trend, intercept, r_value, p_value, std_err = scipy.stats.linregress(yrs[mask], mean_asym[mask])
print('----------')
print('Trend east: ',np.round(trend*100,2),' cm/yr (p=',np.round(p_value,2),')')   


