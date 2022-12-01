# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 15:10:43 2020

@author: IsoldeGlissenaar

Plot and print mean sea ice thickness from CryoSat-2 for all processing methods.
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
from shapely.geometry import Point # Point class
from shapely.geometry import shape # shape() is a function to convert geo objects through the interface
import shapefile
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

# ## IMPORT SIT DATA
cs2_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2/SIT_CryoSat2_March_2011_2020_SnowModelsnowdens.nc')
cs2_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2/SIT_CryoSat2_March_2011_2020_Warrensnowdens.nc')

cs2_cci_sm = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2_cci/SIT_CryoSat2_March_2011_2017_SnowModelsnowdens.nc')
cs2_cci_warren = xr.open_dataset(docs+'Projects/sea_ice_thickness/data/processed/cs2_cci/SIT_CryoSat2_March_2011_2017_Warrensnowdens.nc')

#%%
cs2_mean_warren = np.nanmean(cs2_warren.sit_w99.values[:,:-3], axis=1)
cs2_mean_pmw = np.nanmean(cs2_warren.sit_pmw.values[:,:-3], axis=1)
cs2_mean_stroeve = np.nanmean(cs2_warren.sit_sm.values[:,:-3], axis=1)

cs2_mean_warren_sm = np.nanmean(cs2_sm.sit_w99.values[:,:-3], axis=1)
cs2_mean_pmw_sm = np.nanmean(cs2_sm.sit_pmw.values[:,:-3], axis=1)
cs2_mean_stroeve_sm = np.nanmean(cs2_sm.sit_sm.values[:,:-3], axis=1)

cs2_cci_mean_warren = np.nanmean(cs2_cci_warren.sit_w99.values, axis=1)
cs2_cci_mean_pmw = np.nanmean(cs2_cci_warren.sit_pmw.values, axis=1)
cs2_cci_mean_stroeve = np.nanmean(cs2_cci_warren.sit_sm.values, axis=1)

cs2_cci_mean_warren_sm = np.nanmean(cs2_cci_sm.sit_w99.values, axis=1)
cs2_cci_mean_pmw_sm = np.nanmean(cs2_cci_sm.sit_pmw.values, axis=1)
cs2_cci_mean_stroeve_sm = np.nanmean(cs2_cci_sm.sit_sm.values, axis=1)

#%%
start_yr = 2011; end_yr = 2020

import cartopy.feature as cfeature
import matplotlib.patches as mpatches

land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
lat_corners = np.array([-100., -100, -100, -100, -100, -100, -100., -90, -80, -70, -60, -50., -50.])
lon_corners = np.array([ 50., 55, 60, 65, 70, 75, 80., 80, 80, 80, 80, 80., 50.]) # offset from gridline for clarity
poly_corners = np.zeros((len(lat_corners), 2), np.float64)
poly_corners[:,0] = lat_corners
poly_corners[:,1] = lon_corners
poly = mpatches.Polygon(poly_corners, closed=True, ec='k', linewidth=0.5, fill=False, zorder=10, transform=ccrs.Geodetic())

def plot_sit(lon, lat, sit, title, save=False, save_name=''):
    fig=plt.figure(dpi=200)
    # fig.patch.set_facecolor(None)
    # fig.patch.set_alpha(0)
    ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
    ax.coastlines(resolution='50m',linewidth=0.5)
    ax.set_extent([-78,-43,58,82],crs=ccrs.PlateCarree())
    ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
    im = plt.scatter(lon, lat, c=sit, cmap='Spectral_r',vmin=0,vmax=2.5,s=1.5,transform=ccrs.PlateCarree())
    ax.add_feature(land_50m, facecolor='#eeeeee')
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.title(title,fontsize=8)
    if save:
        plt.savefig('../../figures/figure5/'+save_name+'.png',dpi=300)
        
#%%
#CS2 with Warren '99 density

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_mean_warren,
         f'March mean SIT - CS2 W99 ({str(start_yr)}-{str(end_yr)}) ',
         save_name='cryosat_warren')

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_mean_pmw,
         f'March mean SIT - CS2 PMW ({str(start_yr)}-{str(end_yr)}) ',
         save_name='cryosat_pmw')
    
plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_mean_stroeve,
         f'March mean SIT - CS2 SnowModel ({str(start_yr)}-{str(end_yr)}) ',
         save_name='cryosat_stroeve')


#%%
#CS2 with SnowModel density

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_mean_warren_sm,
         f'March mean SIT - CS2 W99 (sm_dens) ({str(start_yr)}-{str(end_yr)}) ',
         save_name='sm_dens/cryosat_warren')

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_mean_pmw_sm,
         f'March mean SIT - CS2 PMW (sm_dens) ({str(start_yr)}-{str(end_yr)}) ',
         save_name='sm_dens/cryosat_pmw')

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_mean_stroeve_sm,
         f'March mean SIT - CS2 SnowModel (sm_dens) ({str(start_yr)}-{str(end_yr)}) ',
         save_name='sm_dens/cryosat_stroeve')


#%%
#CS2 CCI with Warren '99 density

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_cci_mean_warren,
         f'March mean SIT - CS2 CCI W99 ({str(start_yr)}-{str(end_yr)}) ',
         save_name='cryosat_warren_cci')

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_cci_mean_pmw,
         f'March mean SIT - CS2 CCI PMW ({str(start_yr)}-{str(end_yr)}) ',
         save_name='cryosat_pmw_cci')
    
plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_cci_mean_stroeve,
         f'March mean SIT - CS2 CCI SnowModel ({str(start_yr)}-{str(end_yr)}) ',
         save_name='cryosat_stroeve_cci')


#%%
#CS2 CCI with SnowModel density

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_cci_mean_warren_sm,
         f'March mean SIT - CS2 CCI W99 (sm_dens) ({str(start_yr)}-{str(end_yr)}) ',
         save_name='sm_dens/cryosat_warren_cci')

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_cci_mean_pmw_sm,
         f'March mean SIT - CS2 CCI PMW (sm_dens) ({str(start_yr)}-{str(end_yr)}) ',
         save_name='sm_dens/cryosat_pmw_cci')

plot_sit(cs2_warren.lon.values, cs2_warren.lat.values,
         cs2_cci_mean_stroeve_sm,
         f'March mean SIT - CS2 CCI SnowModel (sm_dens) ({str(start_yr)}-{str(end_yr)}) ',
         save_name='sm_dens/cryosat_stroeve_cci')

#%%
#Find east-west location gridcells
shp = shapefile.Reader(docs+'/QGIS/4quad.shp') #open the shapefile
all_shapes = shp.shapes() # get all the polygons
all_records = shp.records()   

#Get east/west location gridcells on 12 km grid (ICESat-2 & CryoSat-2)
len_f = len(cs2_warren.lon)
eastwest12 = np.zeros(len_f)     
for i in range(len_f):
    pt = (cs2_warren.lon.values[i],cs2_warren.lat.values[i])  
    for k in range (len(all_shapes)):             
        boundary = all_shapes[k]
        if Point(pt).within(shape(boundary)): # make a point and see if it's in the polygon
            eastwest12[i] = all_records[k][0]

#%%

cryosat = np.array([cs2_warren.sit_w99.values, cs2_warren.sit_pmw.values, cs2_warren.sit_sm.values, 
                    cs2_sm.sit_w99.values, cs2_sm.sit_pmw.values, cs2_sm.sit_sm.values, ])

cryosat_uncer = np.array([cs2_warren.sit_w99_uncertainty.values, cs2_warren.sit_pmw_uncertainty.values, cs2_warren.sit_sm_uncertainty.values, 
                          cs2_sm.sit_w99_uncertainty.values, cs2_sm.sit_pmw_uncertainty.values, cs2_sm.sit_sm_uncertainty.values])

cryosat_cci = np.array([cs2_cci_warren.sit_w99.values, cs2_cci_warren.sit_pmw.values, cs2_cci_warren.sit_sm.values, 
                    cs2_cci_sm.sit_w99.values, cs2_cci_sm.sit_pmw.values, cs2_cci_sm.sit_sm.values, ])

cryosat_cci_uncer = np.array([cs2_cci_warren.sit_w99_uncertainty.values, cs2_cci_warren.sit_pmw_uncertainty.values, cs2_cci_warren.sit_sm_uncertainty.values, 
                          cs2_cci_sm.sit_w99_uncertainty.values, cs2_cci_sm.sit_pmw_uncertainty.values, cs2_cci_sm.sit_sm_uncertainty.values])


#%%
cs2_bb = cryosat[:,(eastwest12>0),:]
cs2_bb_mean = np.nanmean(cs2_bb,axis=1)
sit_allmethods = np.nanmean(cs2_bb_mean[:,:-3],axis=1)

cs2_bb_u = cryosat_uncer[:,(eastwest12>0),:]
cs2_bb_u_mean = np.nanmean(cs2_bb_u,axis=1)
uncer_allmethods = np.nanmean(cs2_bb_u_mean[:,:-3],axis=1)

print('---')
print('LARM')
print('-------')
print('Warren snowdensity')
print('March mean (2011-2017) thickness W99: ', np.round(sit_allmethods[0],3),'+- ', np.round(uncer_allmethods[0],3),' m')
print('March mean (2011-2017) thickness PMW: ', np.round(sit_allmethods[1],3),'+- ', np.round(uncer_allmethods[1],3),' m')
print('March mean (2011-2017) thickness Stroeve: ', np.round(sit_allmethods[2],3),'+- ', np.round(uncer_allmethods[2],3),' m')

print('-------')
print('SnowModel snowdensity')
print('March mean (2011-2017) thickness W99: ', np.round(sit_allmethods[3],3),'+- ', np.round(uncer_allmethods[3],3),' m')
print('March mean (2011-2017) thickness PMW: ', np.round(sit_allmethods[4],3),'+- ', np.round(uncer_allmethods[4],3),' m')
print('March mean (2011-2017) thickness Stroeve: ', np.round(sit_allmethods[5],3),'+- ', np.round(uncer_allmethods[5],3),' m')


print('-----')
pmw_vs_w99 = (sit_allmethods[1]-sit_allmethods[0])/sit_allmethods[0]*100
sm_vs_w99 = (sit_allmethods[2]-sit_allmethods[0])/sit_allmethods[0]*100
print('March mean (2011-2017) thickness PMW is ', np.round(pmw_vs_w99,1),'% more than W99')
print('March mean (2011-2017) thickness SnowModel is ', np.round(sm_vs_w99,1),'% more than W99')


print('-------')
mean_sit_w99sdens = np.nanmean(sit_allmethods[0:3])
mean_sit_smsdens = np.nanmean(sit_allmethods[3:6])
w99vssm_dens = (mean_sit_smsdens-mean_sit_w99sdens)/mean_sit_w99sdens*100
print('SnowModel snowdensity SIT is on average ', np.round(w99vssm_dens,1), '% thinner than W99 snowdensity SIT')


print('---')
# w99_vs_pmw_p = (sit_allmethods[1]-sit_allmethods[12])/sit_allmethods[12]*100 
# print(f'Largest difference ICESat: March mean (2003-2009) thickness W99 eff (warren_dens) is {np.round(w99_vs_pmw_p,1)}% more than PMW petty (sm_dens)')


cs2_cci_bb = cryosat_cci[:,(eastwest12>0),:]
cs2_cci_bb_mean = np.nanmean(cs2_cci_bb,axis=1)
sit_allmethods_cci = np.nanmean(cs2_cci_bb_mean[:,:-3],axis=1)

cs2_cci_bb_u = cryosat_cci_uncer[:,(eastwest12>0),:]
cs2_cci_bb_u_mean = np.nanmean(cs2_cci_bb_u,axis=1)
uncer_allmethods_cci = np.nanmean(cs2_cci_bb_u_mean[:,:-3],axis=1)

print('---')
print('CCI')
print('---')
print('Warren snowdensity')
print('March mean (2011-2017) thickness W99: ', np.round(sit_allmethods_cci[0],3),'+- ', np.round(uncer_allmethods_cci[0],3),' m')
print('March mean (2011-2017) thickness PMW: ', np.round(sit_allmethods_cci[1],3),'+- ', np.round(uncer_allmethods_cci[1],3),' m')
print('March mean (2011-2017) thickness Stroeve: ', np.round(sit_allmethods_cci[2],3),'+- ', np.round(uncer_allmethods_cci[2],3),' m')

print('-------')
print('SnowModel snowdensity')
print('March mean (2011-2017) thickness W99: ', np.round(sit_allmethods_cci[3],3),'+- ', np.round(uncer_allmethods_cci[3],3),' m')
print('March mean (2011-2017) thickness PMW: ', np.round(sit_allmethods_cci[4],3),'+- ', np.round(uncer_allmethods_cci[4],3),' m')
print('March mean (2011-2017) thickness Stroeve: ', np.round(sit_allmethods_cci[5],3),'+- ', np.round(uncer_allmethods_cci[5],3),' m')
print('-----')
cci_w99 = (np.nanmean(cs2_cci_mean_warren)-np.nanmean(cs2_mean_warren))/np.nanmean(cs2_mean_warren)*100
cci_pmw = (np.nanmean(cs2_cci_mean_pmw)-np.nanmean(cs2_mean_pmw))/np.nanmean(cs2_mean_pmw)*100
cci_sm = (np.nanmean(cs2_cci_mean_stroeve)-np.nanmean(cs2_mean_stroeve))/np.nanmean(cs2_mean_stroeve)*100
print('March mean (2011-2017) thickness CCI W99 is ', np.round(cci_w99,1),'% more than LARM W99')
print('March mean (2011-2017) thickness CCI PMW is ', np.round(cci_pmw,1),'% more than LARM PMW')
print('March mean (2011-2017) thickness CCI SnowModel is ', np.round(cci_sm,1),'% more than LARM SnowModel')
print('-----')