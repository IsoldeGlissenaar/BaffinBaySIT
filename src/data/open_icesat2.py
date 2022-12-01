# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 13:41:45 2020

@author: IsoldeGlissenaar

Open and prep ICESat2 data
Uses ICESat-2 L3A ATL10 alongtrack freeboard product. Can be retrieved from
NSIDC: https://nsidc.org/data/icesat-2/data-sets
Prep for icesat2_allmethods.py
"""

import xarray as xr
import numpy as np
import h5py
from os import listdir
import os
import pandas as pd
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

freeboard = np.zeros((0))
freeboard_uncertainty = np.zeros((0))
lats = np.zeros((0))
lons = np.zeros((0))
time = np.zeros((0))

year = 2019
month = 3
# is2_file = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/ICESat-2/ATL10-02_20190310233840_11080201_002_01.h5'
direc = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/ICESat-2/'+str(year)+'_'+str("{:02d}".format(month))+'/v3/'

for filename in listdir(direc):
    if filename.endswith('.h5'):
        is2_file = direc+filename
                   
        with h5py.File(is2_file,'r') as is2:
            is2_groups = list(is2.keys())
            fb = is2['/gt1l/freeboard_beam_segment/beam_freeboard/beam_fb_height'][:]
            fb_unc = is2['/gt1l/freeboard_beam_segment/beam_freeboard/beam_fb_sigma'][:]
            dtime = is2['/gt1l/freeboard_beam_segment/beam_freeboard/delta_time'][:] #seconds since 2018-01-01
            lat = is2['/gt1l/freeboard_beam_segment/beam_freeboard/latitude'][:]
            lon = is2['/gt1l/freeboard_beam_segment/beam_freeboard/longitude'][:]
            
            try:
                fb = np.append(fb, is2['/gt1r/freeboard_beam_segment/beam_freeboard/beam_fb_height'][:])
                fb_unc = np.append(fb_unc, is2['/gt1r/freeboard_beam_segment/beam_freeboard/beam_fb_sigma'][:])
                dtime = np.append(dtime, is2['/gt1r/freeboard_beam_segment/beam_freeboard/delta_time'][:])
                lat = np.append(lat, is2['/gt1r/freeboard_beam_segment/beam_freeboard/latitude'][:])
                lon = np.append(lon, is2['/gt1r/freeboard_beam_segment/beam_freeboard/longitude'][:])
            except KeyError:
                print('no gt1r')
                pass
            
            try:
                fb = np.append(fb, is2['/gt2l/freeboard_beam_segment/beam_freeboard/beam_fb_height'][:])
                fb_unc = np.append(fb_unc, is2['/gt2l/freeboard_beam_segment/beam_freeboard/beam_fb_sigma'][:])
                dtime = np.append(dtime, is2['/gt2l/freeboard_beam_segment/beam_freeboard/delta_time'][:])
                lat = np.append(lat, is2['/gt2l/freeboard_beam_segment/beam_freeboard/latitude'][:])
                lon = np.append(lon, is2['/gt2l/freeboard_beam_segment/beam_freeboard/longitude'][:])
            except KeyError:
                print('no gt2l')
                pass
            
            try:
                fb = np.append(fb, is2['/gt2r/freeboard_beam_segment/beam_freeboard/beam_fb_height'][:])
                fb_unc = np.append(fb_unc, is2['/gt2r/freeboard_beam_segment/beam_freeboard/beam_fb_sigma'][:])
                dtime = np.append(dtime, is2['/gt2r/freeboard_beam_segment/beam_freeboard/delta_time'][:])
                lat = np.append(lat, is2['/gt2r/freeboard_beam_segment/beam_freeboard/latitude'][:])
                lon = np.append(lon, is2['/gt2r/freeboard_beam_segment/beam_freeboard/longitude'][:])
            except KeyError:
                print('no gt2r')
                pass
            
            try:
                fb = np.append(fb, is2['/gt3l/freeboard_beam_segment/beam_freeboard/beam_fb_height'][:])
                fb_unc = np.append(fb_unc, is2['/gt3l/freeboard_beam_segment/beam_freeboard/beam_fb_sigma'][:])
                dtime = np.append(dtime, is2['/gt3l/freeboard_beam_segment/beam_freeboard/delta_time'][:])
                lat = np.append(lat, is2['/gt3l/freeboard_beam_segment/beam_freeboard/latitude'][:])
                lon = np.append(lon, is2['/gt3l/freeboard_beam_segment/beam_freeboard/longitude'][:])
            except KeyError:
                print('no gt3l')
                pass
            
            try:
                fb = np.append(fb, is2['/gt3r/freeboard_beam_segment/beam_freeboard/beam_fb_height'][:])
                fb_unc = np.append(fb_unc, is2['/gt3r/freeboard_beam_segment/beam_freeboard/beam_fb_sigma'][:])
                dtime = np.append(dtime, is2['/gt3r/freeboard_beam_segment/beam_freeboard/delta_time'][:])
                lat = np.append(lat, is2['/gt3r/freeboard_beam_segment/beam_freeboard/latitude'][:])
                lon = np.append(lon, is2['/gt3r/freeboard_beam_segment/beam_freeboard/longitude'][:]) 
            except KeyError:
                print('no gt3r')
                pass
            
        
        #Cut Baffin Bay & remove NaNs
        loc = (lat<80)*(lat>50)*(lon<-40)*(lon>-110)*(fb<100)
        
        freeboard = np.append(freeboard, fb[loc])
        freeboard_uncertainty = np.append(freeboard_uncertainty, fb_unc[loc])
        lats = np.append(lats, lat[loc])
        lons = np.append(lons, lon[loc])
        time = np.append(time, dtime[loc])

        print(filename)


# fig=plt.figure(dpi=200)
# ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-60, central_latitude=60, globe=None))
# ax.coastlines(resolution='50m',linewidth=0.5)
# ax.set_extent([-90,-40,58,85],crs=ccrs.PlateCarree()) 
# ax.gridlines(linewidth=0.3, color='k', alpha=0.5, linestyle=':')
# im = plt.scatter(lons, lats, c=freeboard,marker='o',
#                  cmap='jet',vmin=0,vmax=0.5,s=0.1,transform=ccrs.PlateCarree())
# fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


# intialise data of lists.
data = {'lat':lats,
        'lon':lons,
        'time':time,
        'freeboard':freeboard,
        'freeboard_uncertainty': freeboard_uncertainty}
 
# Create DataFrame
is2_fb = pd.DataFrame(data)
is2_fb.to_csv(docs+'/SatelliteData/ICESat-2/'+str(year)+'_'+str("{:02d}".format(month))+'/dataframe/ATL10_'+str(year)+str("{:02d}".format(month))+'.csv',
             index=False)


