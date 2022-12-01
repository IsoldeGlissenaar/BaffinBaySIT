# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 13:49:39 2021

@author: IsoldeGlissenaar

Plot and print mean sea ice thickness from ICESat for all processing methods.
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import xarray as xr
from os import listdir
from shapely.geometry import Point # Point class
from shapely.geometry import shape # shape() is a function to convert geo objects through the interface
import shapefile
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

#%%
# ## IMPORT SIT DATA
is_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is/SIT_ICESat_March_2003_2009_SnowModelsnowdens.nc')
is_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/is/SIT_ICESat_March_2003_2009_Warrensnowdens.nc')


#%%
is_mean_warren = np.nanmean(is_warren.sit_w99.values, axis=1)
is_mean_warren_eff = np.nanmean(is_warren.sit_w99_sig.values, axis=1)
is_mean_warren_petty = np.nanmean(is_warren.sit_w99_piece.values, axis=1)

is_mean_pmw_petty = np.nanmean(is_warren.sit_pmw_piece.values, axis=1)
is_mean_pmw = np.nanmean(is_warren.sit_pmw.values, axis=1)
is_mean_pmw_eff = np.nanmean(is_warren.sit_pmw_sig.values, axis=1)

is_mean_stroeve = np.nanmean(is_warren.sit_sm.values, axis=1)
is_mean_stroeve_eff = np.nanmean(is_warren.sit_sm_sig.values, axis=1)
is_mean_stroeve_petty = np.nanmean(is_warren.sit_sm_piece.values, axis=1)

#%%
#Plot Warren '99

import cartopy.feature as cfeature
import matplotlib.patches as mpatches

land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
lat_corners = np.array([-100., -100, -100, -100, -100, -100, -100., -90, -80, -70, -60, -50., -50.])
lon_corners = np.array([ 50., 55, 60, 65, 70, 75, 80., 80, 80, 80, 80, 80., 50.]) # offset from gridline for clarity
poly_corners = np.zeros((len(lat_corners), 2), np.float64)
poly_corners[:,0] = lat_corners
poly_corners[:,1] = lon_corners


fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m', linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_warren,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
cbar = fig.colorbar(im, ax=ax, label='[m]', fraction=0.046, pad=0.04)
plt.title('March mean ICESat W99 ',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_warren.png',dpi=300)


fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_warren_eff,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
fig.colorbar(im, ax=ax, label='[m]',fraction=0.046, pad=0.04)
plt.title('March mean ICESat W99 effsd ',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_warren_effsd.png',dpi=300)


fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_warren_petty,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
fig.colorbar(im, ax=ax, label='[m]',fraction=0.046, pad=0.04)
plt.title('March mean ICESat W99 Petty',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_warren_petty.png',dpi=300)

#%%
#Plot pmw

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_pmw,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.title('March mean ICESat PMW ',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_pmw.png',dpi=300)


fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_pmw_eff,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.title('March mean ICESat PMW effsd ',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_pmw_effsd.png',dpi=300)


fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_pmw_petty,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
fig.colorbar(im, ax=ax, label='[m]', fraction=0.046, pad=0.04)
plt.title('March mean ICESat PMW Petty',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_pmw_petty.png',dpi=300)

#%%
#Plot Stroeve

fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_stroeve,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.title('March mean ICESat Stroeve ',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_stroeve.png',dpi=300)


fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_stroeve_eff,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.title('March mean ICESat Stroeve effsd ',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_stroeve_effsd.png',dpi=300)


fig=plt.figure(dpi=200)
# fig.patch.set_facecolor(None)
# fig.patch.set_alpha(0)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
ax.coastlines(resolution='50m',linewidth=0.5)
ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree()) 
ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
im = plt.scatter(is_warren.lon.values, is_warren.lat.values, c=is_mean_stroeve_petty,cmap='Spectral_r',vmin=0,vmax=2.5,s=9,transform=ccrs.PlateCarree())
ax.add_feature(land_50m, facecolor='#eeeeee')
ax.add_patch(mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic()))
fig.colorbar(im, ax=ax, label='[m]', fraction=0.046, pad=0.04)
plt.title('March mean ICESat Stroeve Petty ',fontsize=8)
# plt.savefig('../../figures/figure5/icesat_stroeve_petty_cbar.png',dpi=300)


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


icesat_west = icesat[:,(eastwest25==3)|(eastwest25==4),:]
icesat_east = icesat[:,(eastwest25==1)|(eastwest25==2),:]

is_west_mean = np.nanmean(icesat_west,axis=1)
is_east_mean = np.nanmean(icesat_east,axis=0)


#%%
is_bb = icesat[:,(eastwest25>0),:]
is_bb_mean = np.nanmean(is_bb,axis=1)
sit_allmethods = np.nanmean(is_bb_mean,axis=1)

is_bb_u = icesat_uncer[:,(eastwest25>0),:]
is_bb_u_mean = np.nanmean(is_bb_u,axis=1)
uncer_allmethods = np.nanmean(is_bb_u_mean,axis=1)

print('-------')
print('Warren snowdensity')
print('March mean (2003-2009) thickness W99: ', np.round(sit_allmethods[0],3),'+- ', np.round(uncer_allmethods[0],3),' m')
print('March mean (2003-2009) thickness W99+eff: ', np.round(sit_allmethods[1],3),'+- ', np.round(uncer_allmethods[1],3),' m')
print('March mean (2003-2009) thickness W99+Petty: ', np.round(sit_allmethods[2],3),'+- ', np.round(uncer_allmethods[2],3),' m')
print('March mean (2003-2009) thickness PMW: ', np.round(sit_allmethods[3],3),'+- ', np.round(uncer_allmethods[3],3),' m')
print('March mean (2003-2009) thickness PMW+eff: ', np.round(sit_allmethods[4],3),'+- ', np.round(uncer_allmethods[4],3),' m')
print('March mean (2003-2009) thickness PMW+Petty: ', np.round(sit_allmethods[5],3),'+- ', np.round(uncer_allmethods[5],3),' m')
print('March mean (2003-2009) thickness Stroeve: ', np.round(sit_allmethods[6],3),'+- ', np.round(uncer_allmethods[6],3),' m')
print('March mean (2003-2009) thickness Stroeve+eff: ', np.round(sit_allmethods[7],3),'+- ', np.round(uncer_allmethods[7],3),' m')
print('March mean (2003-2009) thickness Stroeve+Petty: ', np.round(sit_allmethods[8],3),'+- ', np.round(uncer_allmethods[8],3),' m')

print('-------')
print('SnowModel snowdensity')
print('March mean (2003-2009) thickness W99: ', np.round(sit_allmethods[9],3),'+- ', np.round(uncer_allmethods[9],3),' m')
print('March mean (2003-2009) thickness W99+eff: ', np.round(sit_allmethods[10],3),'+- ', np.round(uncer_allmethods[10],3),' m')
print('March mean (2003-2009) thickness W99+Petty: ', np.round(sit_allmethods[11],3),'+- ', np.round(uncer_allmethods[11],3),' m')
print('March mean (2003-2009) thickness PMW: ', np.round(sit_allmethods[12],3),'+- ', np.round(uncer_allmethods[12],3),' m')
print('March mean (2003-2009) thickness PMW+eff: ', np.round(sit_allmethods[13],3),'+- ', np.round(uncer_allmethods[13],3),' m')
print('March mean (2003-2009) thickness PMW+Petty: ', np.round(sit_allmethods[14],3),'+- ', np.round(uncer_allmethods[14],3),' m')
print('March mean (2003-2009) thickness Stroeve: ', np.round(sit_allmethods[15],3),'+- ', np.round(uncer_allmethods[15],3),' m')
print('March mean (2003-2009) thickness Stroeve+eff: ', np.round(sit_allmethods[16],3),'+- ', np.round(uncer_allmethods[16],3),' m')
print('March mean (2003-2009) thickness Stroeve+Petty: ', np.round(sit_allmethods[17],3),'+- ', np.round(uncer_allmethods[17],3),' m')


print('-----')
pmw_vs_w99 = (sit_allmethods[3]-sit_allmethods[0])/sit_allmethods[0]*100
pmw_vs_w99_e = (sit_allmethods[4]-sit_allmethods[1])/sit_allmethods[1]*100 
pmw_vs_w99_p = (sit_allmethods[5]-sit_allmethods[2])/sit_allmethods[2]*100 
sm_vs_w99 = (sit_allmethods[6]-sit_allmethods[0])/sit_allmethods[0]*100
sm_vs_w99_e = (sit_allmethods[7]-sit_allmethods[1])/sit_allmethods[1]*100
sm_vs_w99_p = (sit_allmethods[8]-sit_allmethods[2])/sit_allmethods[2]*100
print('March mean (2003-2009) thickness PMW is ', np.round(pmw_vs_w99,1),'% more than W99')
print('March mean (2003-2009) thickness PMW is ', np.round(pmw_vs_w99_e,1),'% more than W99 effective')
print('March mean (2003-2009) thickness PMW is ', np.round(pmw_vs_w99_p,1),'% more than W99 petty')
print('March mean (2003-2009) thickness SnowModel is ', np.round(sm_vs_w99,1),'% more than W99')
print('March mean (2003-2009) thickness SnowModel is ', np.round(sm_vs_w99_e,1),'% more than W99 effective')
print('March mean (2003-2009) thickness SnowModel is ', np.round(sm_vs_w99_p,1),'% more than W99 petty')

print('---')
pmw_vs_pmw_e = (sit_allmethods[4]-sit_allmethods[3])/sit_allmethods[3]*100 
pmw_vs_pmw_p = (sit_allmethods[5]-sit_allmethods[3])/sit_allmethods[3]*100 
sm_vs_sm_p = (sit_allmethods[8]-sit_allmethods[6])/sit_allmethods[6]*100 
sm_vs_sm_e = (sit_allmethods[7]-sit_allmethods[6])/sit_allmethods[6]*100
w99_vs_w99_p = (sit_allmethods[2]-sit_allmethods[0])/sit_allmethods[0]*100 
w99_vs_w99_e = (sit_allmethods[1]-sit_allmethods[0])/sit_allmethods[0]*100
print('March mean (2003-2009) thickness PMW Eff is ', np.round(pmw_vs_pmw_e,1),'% more than PMW')
print('March mean (2003-2009) thickness PMW Petty is ', np.round(pmw_vs_pmw_p,1),'% more than PMW')
print('March mean (2003-2009) thickness SnowModel Eff is ', np.round(sm_vs_sm_e,1),'% more than SnowModel')
print('March mean (2003-2009) thickness SnowModel Petty is ', np.round(sm_vs_sm_p,1),'% more than SnowModel')
print('March mean (2003-2009) thickness W99 Eff is ', np.round(w99_vs_w99_e,1),'% more than W99')
print('March mean (2003-2009) thickness W99 Petty is ', np.round(w99_vs_w99_p,1),'% more than W99')


print('-------')
mean_sit_w99sdens = np.nanmean(sit_allmethods[0:9])
mean_sit_smsdens = np.nanmean(sit_allmethods[9:18])
w99vssm_dens = (mean_sit_smsdens-mean_sit_w99sdens)/mean_sit_w99sdens*100
print('SnowModel snowdensity SIT is on average ', np.round(w99vssm_dens,1), '% thinner than W99 snowdensity SIT')


print('---')
w99_vs_pmw_p = (sit_allmethods[1]-sit_allmethods[12])/sit_allmethods[12]*100 
print(f'Largest difference ICESat: March mean (2003-2009) thickness W99 eff (warren_dens) is {np.round(w99_vs_pmw_p,1)}% more than PMW petty (sm_dens)')

