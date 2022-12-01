# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 10:29:48 2019

@author: IsoldeGlissenaar

Calculate alongtrack total freeboard from ICESat measurements (2003-2009) 
according to Landy et al. (2017).
"""

import sys
sys.path.append('../functions')
from func_reproj import reproject
from func_smooth import smooth

import numpy as np
import scipy.io
from scipy import spatial
from netCDF4 import Dataset
from scipy.interpolate import Rbf, PchipInterpolator
from datetime import date
docs = 'C:/Users/zq19140/OneDrive - University of Bristol/Documents/'

sic_folder=docs+'SatelliteData/OSI-SAF_sea_ice_concentration/'
data = scipy.io.loadmat(docs+'Jack_onedrive/Existing Sea Ice Thickness Datasets/Landy et al 2017 ICESat GLAS13/GLAH13_034_ICESat_data_Hudson_Bay_V5')['results']
## [id shot year month day hour min lat lon elev 50km_running_median elevar reflectance geoid_height delta_ellipsoid]

data[:,7],data[:,8] = reproject(data[:,7],data[:,8])


#%%
## Add SIC to icesat data and filter sic

# Identify all icesat points for one day and import SIC
U = np.unique(data[:,2:5],axis=0)
icesat_sic = np.zeros((len(data[:,0])))
for i in range (0,len(U[:,0])):
    idx = np.where((data[:,2]==U[i,0]) & (data[:,3]==U[i,1]) & (data[:,4]==U[i,2]))[0]
    sic_id = '_'+str("{:.0f}".format(U[i,0]))+str("{:02d}".format(int(U[i,1])))+str("{:02d}".format(int(U[i,2])))
    
    sic_filename = sic_folder+str("{:.0f}".format(U[i,0]))+'/'+str("{:02d}".format(int(U[i,1])))+'/ice_conc_nh_ease2-250_cdr-v2p0'+sic_id+'1200.nc'
    sic_lat = Dataset(sic_filename)['lat'][:,:]
    sic_lon = Dataset(sic_filename)['lon'][:,:]
    sic     = Dataset(sic_filename)['ice_conc'][0,:,:]/100
    
    idlat = ((sic_lat>45) & (sic_lat<85))
    idlon = ((sic_lon>-100) & (sic_lon<-40))
    ids = ((idlat==True) & (idlon==True) & (sic>=0))
    sic_lat=sic_lat[ids]; sic_lon=sic_lon[ids]; sic=sic[ids]
    
    sic_lat,sic_lon = reproject(sic_lat,sic_lon)
    
    interp = Rbf(sic_lat, sic_lon ,sic, kind='linear')
    icesat_sic[idx] = interp(data[idx,7],data[idx,8])
    
    #remove icesat_sic values from cells where closest sic data is >25km
    points = np.array(list(zip(sic_lat,sic_lon)))
    tree = spatial.KDTree(points)
    Idx = tree.query_ball_point(data[idx,7:9], 25000)
    empty = np.invert([any(idxi) for idxi in Idx])
    id_empty = np.where(empty)
    icesat_sic[idx[id_empty[0]]] = 0

    print(f"Got to {i} out of {len(U[:,0])}")
    

# Remove points with sic<0.2 (Kwok et al 2007 use sic<0.3, NSIDC use sic<0.2)
data = data[icesat_sic>0.2,:]
icesat_sic = icesat_sic[icesat_sic>0.2]

# Remove days with <100 points
U = np.unique(data[:,2:5],axis=0)
N_samples = np.zeros((len(data[:,0])))
for i in range (len(U[:,0])):
    ids = (data[:,2]==U[i,0]) & (data[:,3]==U[i,1]) & (data[:,4]==U[i,2])
    N_samples[ids] = (np.sum(ids)>100)
data = data[N_samples>0,:] ; icesat_sic = icesat_sic[N_samples>0]

#np.save('Data/data_cell1.npy', data); np.save('Data/icesat_sic_cell1.npy',icesat_sic)

#%%
## Calculate relative elevation and sea level
# icesat_sic = np.load('Data/icesat_sic_cell1.npy')
# data = np.load('Data/data_cell1.npy')

icesat_sic = icesat_sic[data[:,9]<1000]
data = data[data[:,9]<1000]

# Calculate sea level from tiepoints assessed point-to-point,
# including only points with reflectivity >0.2 and elevation variation
# >X below background level [Kwok and Cunningham 2008]
# [results col13] = calibrated reflectance
KandC08_function = scipy.io.loadmat(docs+'Jack_onedrive/Existing Sea Ice Thickness Datasets/Landy et al 2017 ICESat GLAS13/KandC08_function')['KandC08_function']
datnum = np.zeros(len(data[:,0]))
for i in range(len(data[:,0])):
    datnum[i] = date.toordinal(date(int(data[i,2]), int(data[i,3]), int(data[i,4]))) + 366 + (data[i,5])/24. + (data[i,6])/24/60.

neighbourhood=25000  #Size of search neighbourhood in m, default=25km
U = np.unique(data[:,2:5], axis=0)
rel_elev = np.zeros(len(data[:,0]))
sea_level_v1 = np.zeros(len(data[:,0]))
sea_level_v2 = np.zeros(len(data[:,0]))
sea_level_uncertainty = np.zeros(len(data[:,0]))
dist_closest_lead = np.zeros(len(data[:,0]))
# tie_points_arr = np.zeros(len(data[:,0]))

for i in range(0,len(U[:,0])):
    ids = np.nonzero((data[:,2]==U[i,0]) & (data[:,3]==U[i,1]) & (data[:,4]==U[i,2]))[0]
    
    # Identify different segments within day
    q = np.nonzero(np.diff(datnum[ids])>0.025)
    if len(q[0])==0:
        jumps = [0]; jumps.append(len(ids)-1)
    elif q[0][0]==0:
        jumps = [];  jumps.extend(list(q[0]));  jumps.append(len(ids)-1)
    else:
        jumps = [0];  jumps.extend(list(q[0]));  jumps.append(len(ids)-1)
    for j in range(1,len(jumps)):
        idss = ids[(ids>ids[jumps[j-1]]) & (ids<=ids[jumps[j]])]
        #point = np.array(list(zip(data[idss,7],data[idss,8])))
        MdlKDT = spatial.KDTree(data[idss,7:9])
        Idx = MdlKDT.query_ball_point(data[idss,7:9], neighbourhood)
        
        elev = data[idss,9]
        if len(elev)>150:
            elev_smooth = smooth(elev, 147)
        else:
            elev_smooth = np.nanmean(elev)
        elev_segment = np.array(elev-elev_smooth)
        
        #Remove outliers from relative elevation >4std from mean
        outliers = (elev_segment > (np.nanmean(elev_segment) + 4*np.nanstd(elev_segment))) | (
                    elev_segment < (np.nanmean(elev_segment) - 4*np.nanstd(elev_segment)))
        elev_segment[outliers]=0
        
        rel_elev[idss]=elev_segment[:]
        R_segment=data[idss,12]
        
        # Calculate sea level v1 (lowest 5% points within 50 km)
        slv1_idx = MdlKDT.query_ball_point(data[idss,7:9],50000)
        slv1_temp = [np.sort(elev_segment[idxx],axis=0) for idxx in slv1_idx]
        sea_level_v1[idss] = [np.median(x[0:int(np.round(0.05*len(x)))]) for x in slv1_temp]
        
        # Identify leads
        R_mean = np.array([np.mean(R_segment[x]) for x in Idx])
        R_std  = np.array([np.std(R_segment[x],ddof=1) for x in Idx])
        R_threshold = R_mean - 1.5*R_std
        R_bg = [np.mean(R_segment[np.array(x)[np.where(R_segment[x]>R_threshold[x])[0]]]) for x in Idx]
        R_delta = R_bg-R_segment
        
        h_mean = np.array([np.mean(elev_segment[x]) for x in Idx])
        h_std  = np.array([np.std(elev_segment[x],ddof=1) for x in Idx])
        h_threshold = h_mean - 1.5*h_std
        h_bg = np.array([np.mean(elev_segment[np.array(x)[np.where(elev_segment[x]>h_threshold[x])[0]]]) for x in Idx])
        h_d = h_bg-elev_segment

        tie_points = ((R_delta>0.2) & (h_d>0))  #Kwok et al (2007) use Rd>0.3 and hd>0.3
        
        if np.sum(tie_points)>0:
            tp_locs = np.where(tie_points==True)[0]
            R_delta_ties = R_delta[tie_points]
            sd_correction=np.zeros(np.sum(tie_points))
            for k in range(0,np.sum(tie_points)):
                sd_correction[k] = KandC08_function[np.round(R_delta_ties[k],1)==KandC08_function[:,0],1]  #K&C2008 snow depth correction
                
            # Remove outlying leads
            D = np.sqrt((data[idss,7]-data[idss[0],7])**2 + (data[idss,8]-data[idss[0],8])**2)  
            elev_tp = elev_segment[tie_points]+sd_correction
            tp_outliers = (elev_tp>(np.median(elev_tp)+3*np.std(elev_tp))) | (
                           elev_tp<(np.median(elev_tp)-2*np.std(elev_tp)))
            tie_points[tp_locs[tp_outliers]] = 0
            elev_tp = elev_tp[tp_outliers<1]
            
            # Calculate sea level v2
            if tie_points[0]<1 and tie_points[-1]<1:
                tie_points[0]=1; tie_points[-1]=1
                elev_tp2 = [elev_tp[0]]; elev_tp2.extend(list(elev_tp)); elev_tp2.append(elev_tp[-1])
            elif tie_points[0]>0 and tie_points[-1]<1:
                tie_points[0]=1; tie_points[-1]=1
                elev_tp2 = list(elev_tp); elev_tp2.append(elev_tp[-1])
            elif tie_points[0]<1 and tie_points[-1]>0:
                tie_points[0]=1; tie_points[-1]=1
                elev_tp2=[elev_tp[0]]; elev_tp2.extend(list(elev_tp))
            else:
                tie_points[0]=1; tie_points[-1]=1
                elev_tp2 = elev_tp[:]
            
            combo = np.transpose(np.array([D[tie_points],elev_tp2]))
            test = combo[combo[:,0].argsort()]
            sl_interpl = PchipInterpolator(test[:,0], test[:,1])#  interp1d, kind='cubic', fill_value=np.nan) #"extrapolate")
            sl = sl_interpl(D)

            if len(sl)>585:
                sl_smooth = smooth(sl,581)
            else:
                sl_smooth = np.nanmean(sl)
            sea_level_v2[idss] = sl_smooth

            # Calculate uncertainty in sl over 25 km length scale
            sea_level_uncertainty[idss] = np.array([np.std(sl[x]) for x in Idx])
            
            # Calculate distance to closest lead [m]
            dist_closest_lead[idss] = np.array([np.min(np.abs(x-D[tie_points])) for x in D])
            
            # tie_points_arr[idss] = tie_points
    print(f"Got to {i} out of {len(U[:,0])}")
    
            
h_delta = data[:,13]-sea_level_v2  # diff between ICESat sea level and geoid
freeboard = rel_elev-sea_level_v2  # ice freeboard

freeboard[freeboard<0] = 0 # set negative fb's to zero, according to NSIDC  


# np.save('Data/data_cell2.npy', data); np.save('Data/icesat_sic_cell2.npy',icesat_sic)
# np.save('Data/sea_level_v1_cell2.npy', sea_level_v1); np.save('Data/sea_level_v2_cell2.npy', sea_level_v2)
# np.save('Data/sea_level_uncertainty_cell2.npy', sea_level_uncertainty)
# np.save('Data/freeboard_cell2.npy', freeboard); np.save('Data/h_delta_cell2.npy', h_delta)
# np.save('Data/rel_elev_cell2.npy', rel_elev); np.save('Data/dist_closest_lead_cell2.npy', dist_closest_lead)
# np.save('Data/tie_points_cell2.npy', tie_points_arr)
