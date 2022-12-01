# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 15:00:23 2020

@author: IsoldeGlissenaar

Open and prep CryoSat-2 CCI data
Uses CryoSat-2 ESA CCI freeboard observations. Can be retrieved from:
https://catalogue.ceda.ac.uk/uuid/5b6033bfb7f241e89132a83fdc3d5364
Prep for cryosat2_cci_allmethods.py
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd
from os import listdir
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

yr = 2011
month = 3    
direc = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/SatelliteData/CryoSat-2_CCI/'+str(yr)+'/'

day_arr = []
lat_arr = []
lon_arr = []
yr_arr = []
month_arr = []
fb_arr = []
fb_uncer_arr = []
for filename in listdir(direc):
    if filename.endswith('.nc'):
        file = Dataset(direc+filename)

        loc = (file['lat'][:]<80)*(file['lat'][:]>50)*(file['lon'][:]<-40)*(file['lon'][:]>-110)

        day = int(filename[50:52])
        lat = file['lat'][:][loc]
        lon = file['lon'][:][loc]
        fb = file['radar_freeboard'][:][loc] #Radar freeboard
        fb_uncer = file['radar_freeboard_uncertainty'][:][loc]
        
        lat_arr = np.append(lat_arr, lat)
        lon_arr = np.append(lon_arr, lon)
        day_arr = np.append(day_arr, np.full(len(lat), day))         
        yr_arr = np.append(yr_arr, np.full(len(lat), yr))        
        month_arr = np.append(month_arr, np.full(len(lat), month))
        fb_arr = np.append(fb_arr, fb)
        fb_uncer_arr = np.append(fb_uncer_arr, fb_uncer)

 
# intialise data of lists.
data = {'lat':lat_arr,
        'lon':lon_arr,
        'yr':yr_arr,
        'month':month_arr,
        'day':day_arr,
        'freeboard':fb_arr,
        'freeboard_uncertainty':fb_uncer_arr}
 
# Create DataFrame
cs_cci = pd.DataFrame(data)
cs_cci.to_csv(docs+'/SatelliteData/CryoSat-2_CCI/dataframe/cs2_cci_'+str(yr)+str("{:02d}".format(month))+'.csv',
             index=False)


