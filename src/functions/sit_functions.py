# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 13:03:13 2021

@author: zq19140
"""

import numpy as np
from func_easegrid import easegrid
from func_ieasegrid import ieasegrid
from func_reproj import reproject
from math import floor
from scipy import spatial
from os import listdir
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def create_grid(resolution=25):
    # Function create_grid(resolution)
    # Input: resolution
    # Output: reslat_plt, reslon_plt, reslat, reslon
    #
    # creates ease grid for Baffin Bay region
    #
    # given the required resolution (in km), creates the coordinates for
    # easegrid in Baffin Bay. [reslat_plt, reslon_plt] give the coordinates 
    # in the WGS84 coordinate system (for plotting purposes), [reslat, reslot]
    # give the coordinates in EPSG:3143 (NSIDC Sea Ice Polar Stereographic North)
    # Used to create easegrid for the sea ice thickness product.
    
    extent = np.array([[50,0],[50,90],[50,180],[50,270]])
    [thelon, thelat] = easegrid(11,extent[:,0],extent[:,1],6371.228/25)
    x = np.linspace(np.min(thelat), np.max(thelat), floor((np.max(thelat)-np.min(thelat))/(resolution/25)))
    y = np.linspace(np.min(thelon), np.max(thelon), floor((np.max(thelon)-np.min(thelon))/(resolution/25)))
    [x,y] = np.meshgrid(x,y)
    [reslon, reslat] = ieasegrid(11,y,x,6371.228/25)
    reslat=reslat[:]
    reslon=reslon[:]
    # 50<lat<80 & -100<lon<-50
    idl = ((reslat>=50) & (reslat<=80) & (reslon>=-100) & (reslon<=-50))
    reslat_plt = reslat[idl]
    reslon_plt = reslon[idl]
    
    reslat, reslon = reproject(reslat_plt, reslon_plt)
    
    return reslat_plt, reslon_plt, reslat, reslon


def sit_laser(data, ow_density=1023.9, si_density_v1=915.1):
    # Function sit_laser(data, ow_density, si_density_v1)
    # Input: data, ow_density, si_density_v1
    # Output: sea_ice_thickness_v2, si_density_v2, sea_ice_thickness_v1
    #
    # determines sea ice thickness from total freeboard observations
    #
    # given data (DataFrame including total freeboard, snow depth and 
    # snow density), ocean water density (ow_density) and a first estimate
    # of sea ice density (si_density_v1), determines sea ice thickness.
    
    sea_ice_thickness_v1 = (ow_density/(ow_density-si_density_v1))*data.freeboard.values - (
                       (ow_density-data.snowdensity.values)/(ow_density-si_density_v1))*data.snowdepth.values
    si_density_v2 = (0.93633-0.0018*sea_ice_thickness_v1**0.5)*1000  # ice-thickness dependent param, Kovacs [1996]
    si_density_v2[sea_ice_thickness_v1<0]=936.33                     #density when sea ice thickness=0
    sea_ice_thickness_v2 = (ow_density/(ow_density-si_density_v2))*data.freeboard.values - (
                       (ow_density-data.snowdensity.values)/(ow_density-si_density_v2))*data.snowdepth.values
    
    return(sea_ice_thickness_v2,si_density_v2,sea_ice_thickness_v1)


def sit_radar(data, ow_density=1023.9, si_density_v1=915.1):
    # Function sit_radar(data, ow_density, si_density_v1)
    # Input: data, ow_density, si_density_v1
    # Output: sea_ice_thickness_v2, si_density_v2, sea_ice_thickness_v1
    #
    # determines sea ice thickness from radar freeboard observations
    #
    # given data (DataFrame including radar freeboard, snow depth and 
    # snow density), ocean water density (ow_density) and a first estimate
    # of sea ice density (si_density_v1), determines sea ice thickness.
    
    c = 299792458  # speed of light in m/s
    c_snow = c*(1+0.51*(data.snowdensity.values/1000))**(-1.5)
    fb_i = data.freeboard.values + data.snowdepth.values*((c/c_snow)-1) #(1-(c_snow/c))
    
    sea_ice_thickness_v1 = (ow_density/(ow_density-si_density_v1))*fb_i + (
                       (data.snowdensity.values)/(ow_density-si_density_v1))*(data.snowdepth.values)
    si_density_v2 = (0.93633-0.0018*sea_ice_thickness_v1**0.5)*1000  # ice-thickness dependent param, Kovacs [1996]
    si_density_v2[sea_ice_thickness_v1<0]=936.33                     #density when sea ice thickness=0
    sea_ice_thickness_v2 = (ow_density/(ow_density-si_density_v2))*fb_i + (
                           (data.snowdensity.values)/(ow_density-si_density_v2))*(data.snowdepth.values)
    
    return(sea_ice_thickness_v2,si_density_v2,sea_ice_thickness_v1)


def petty(hf, hsl):
    # Function petty(hf, hsl)
    # Input: hf, hsl
    # Output: hs
    #
    # redistributes snow depth according to Petty et al. (2020)
    #
    # given freeboard observations (hf) and large scale snow depth (hsl),
    # returns redistributeed snow along the piecewise function described by 
    # Petty et al. (2020).

    hs = np.zeros(len(hf))
    
    #Petty et al., (in review)
    c1=0.70; c2=0.22; c3=0.16
    c4=1.03; c5=0.01
    hs_thick = c4*hsl + c5
    hf_cutoff = c1*hsl + c2*np.nanmean(hf) + c3
    
    thick = (hf>hf_cutoff)
    hs[thick] = hs_thick
    hs[~thick] = (hf[~thick]/hf_cutoff)*hs_thick
    
    hs[hs>hf] = hf[hs>hf]

    
    diff = hsl-np.nanmean(hs)
    count = 0
    while np.abs(diff)>0.01 and count<10:
        hs_thick = hs_thick + diff
        if hs_thick<0:
            hs_thick=0
        hs[thick] = hs_thick
        hs[~thick] = (hf[~thick]/hf_cutoff)*hs_thick
        hs[hs>hf] = hf[hs>hf]

        diff = hsl-np.nanmean(hs)
        count = count+1

    return(hs)



def gridding(reslat, reslon, data, uncer, search_radius=50, speckle=False):
    # Function gridding(reslat, reslon, data, uncer, search_radius, speckle)
    # Input: reslat, reslon, data, uncer, search_radius, speckle
    # Output: sea_ice_grid, uncertainty_grid
    #
    # grids sea ice thickness and uncertainty
    #
    # given lat and lon of desired grid (reslat, reslon), along-track
    # data including lat, lon, and sea ice thickness (data), the uncertainty 
    # estimate of the along-track data (uncer), a search radius (search_radius,
    # in km), and whether or not speckle uncertainty is included (speckle),
    # grids sea ice thickness on the given grid according to an inverse
    # linear mean filter (weighted by sea level uncertainty and distance)
    
    sea_ice_grid = np.zeros((len(reslat[:])))
    uncertainty_grid = np.zeros((len(reslat[:])))
    coord = np.transpose(np.array([data.lat.values, data.lon.values]))
    MdlKDT = spatial.KDTree(coord)
    D, Idx = MdlKDT.query(np.transpose(np.array([reslat, reslon])), k=10000, distance_upper_bound=search_radius*1e3)
    maxD = np.max(np.where(D<50e3)[1])
    D = D[:,:maxD];    Idx = Idx[:,:maxD].astype(float)
    Idx[np.where(D>50e3)]=np.nan;    D[np.where(D>50e3)]=np.nan
    lenD = np.sum(~np.isnan(D),axis=1)
    
    for j in range(0,len(Idx[:,0])):
        x = Idx[j,(~np.isnan(Idx[j,:]))].astype(int)
        y = D[j,(~np.isnan(D[j,:]))]
        a = np.nanmean(data.sit.values[x]*(1/np.transpose(y))*(1/uncer['total_it_uncertainty'].values[x]))
        b = np.mean(1/y)*np.nanmean(1/uncer['total_it_uncertainty'].values[x])
        sea_ice_grid[j] = a/b
        c_fb = np.nanmean((uncer['fb'].values[x]*(1/y))/np.nanmean(1/y))
        c_sd = np.nanmean((uncer['sd'].values[x]*(1/y))/np.nanmean(1/y))
        c_sdens = np.nanmean((uncer['snow_dens'].values[x]*(1/y))/np.nanmean(1/y))
        c_si = np.nanmean((uncer['si_dens'].values[x]*(1/y))/np.nanmean(1/y))
        c_ow = np.nanmean((uncer['ow_dens'].values[x]*(1/y))/np.nanmean(1/y))
        no_tracks = np.unique(data['track'].values[x]).shape
        d_track = 1/np.sqrt(no_tracks)
        d = 1/np.sqrt(lenD[j])
        fb_uncer_grid = c_fb*d_track
        sd_uncer_grid = c_sd
        sdens_uncer_grid = c_sdens
        si_uncer_grid = c_si*d
        ow_uncer_grid = c_ow*d
        if speckle:
            c_speckle = np.nanmean((uncer['speckle'].values[x]*(1/y))/np.nanmean(1/y))
            speckle_uncer_grid = c_speckle*d
        else:
            speckle_uncer_grid=0
        uncertainty_grid[j] = np.sqrt(fb_uncer_grid**2+speckle_uncer_grid**2+sd_uncer_grid**2+sdens_uncer_grid**2+si_uncer_grid**2+ow_uncer_grid**2)
        
    return (sea_ice_grid, uncertainty_grid)
    
    
def uncertainty(data, ow_density=1023.9, si_density_v1=915.1, icesat2=False, speckle=False):
    # Function uncertainty(data, ow_density, si_density_v1, icesat2, speckle)
    # Input: data, ow_density, si_density_v1, icesat2, speckle
    # Output: uncer
    #
    # determines uncertainty estimate for along-track freeboard obs.
    #
    # given data (including freeboard uncertainty, snow depth, snow density,
    # sea ice thickness and freeboard from along-track observations), ocean 
    # water density (ow_density), and sea ice density (si_density_v1), determines
    # freeboard, snow depth, snow density, sea ice density, ocean water density
    # and total uncertainty estimates, combined in uncer. If freeboard_uncertainty
    # is included in data icesat2=True, otherwise, freeboard uncertainty is
    # determined from sea level uncertainty. If speckle uncertainty is included
    # in data (CryoSat-2), speckle=True.
    
    import pandas as pd
    uncer = pd.DataFrame()
    if icesat2:
        uncer['fb'] = data.freeboard_uncertainty.values
    else:
        uncer['fb'] = (data.sea_level_uncertainty)*(ow_density/(ow_density-si_density_v1))
    uncer['sd'] = (0.06)*(-(ow_density-data.snowdensity.values)/(ow_density-si_density_v1))
    uncer['snow_dens'] = (100)*(data.snowdepth.values/(ow_density-si_density_v1))
    uncer['si_dens'] = (15)*(data.sit.values/(ow_density-si_density_v1))
    uncer['ow_dens'] = (1)*((-si_density_v1*data.freeboard+(
                            si_density_v1-data.snowdensity)*data.snowdepth)/(ow_density-si_density_v1)**2)
    
    if speckle:
        uncer['speckle'] = np.zeros(len(data))
        uncer.speckle = uncer.speckle.where(data.sar_mode != 0 , 0.1)
        uncer.speckle = uncer.speckle.where(data.sar_mode != 1 , 0.14)
        uncer['total_it_uncertainty'] = np.sqrt(uncer['fb']**2+uncer['speckle']**2+uncer['sd']**2+uncer['snow_dens']**2+uncer['si_dens']**2+uncer['ow_dens']**2)
        return uncer
    
    uncer['total_it_uncertainty'] = np.sqrt(uncer['fb']**2+uncer['sd']**2+uncer['snow_dens']**2+uncer['si_dens']**2+uncer['ow_dens']**2)
    return uncer


def freeboard_uncertainty(data, distance_to_lead, ow_density=1023.9, si_density_v1=915.1):
    # Function freeboard_uncertainty(data, distance_to_lead, ow_density, si_density_v1)
    # Input: data, distance_to_lead, ow_density, si_density_v1
    # Output: fb_uncertainty
    #
    # determines uncertainty estimate for along-track freeboard obs.
    #
    # given data (including sea level uncertainty, icesat height and sea level),
    # distance to nearest lead, and ocean water and sea ice density, determines freeboard
    # uncertainty. Used to estimate freeboard uncertainty for ICESat. 
    
    fb_uncertainty = np.zeros(len(data)); fb_uncertainty[:]=np.nan
    under = (distance_to_lead<=25000)
    fb_uncertainty[under] = (data.sea_level_uncertainty[under])*(ow_density/(ow_density-si_density_v1))
    upper = (distance_to_lead>100000)
    fb_uncertainty[upper] = np.abs(data.icesat_height[upper]-data.sea_level[upper])
    middle = (distance_to_lead>25000)*(distance_to_lead<=100000)
    bot = (data.sea_level_uncertainty[middle])*(ow_density/(ow_density-si_density_v1))
    top = np.abs(data.icesat_height[middle]-data.sea_level[middle])
    fb_uncertainty[middle] = (top-bot)/75000.*distance_to_lead[middle]+1.333*bot-0.333*top
    return(fb_uncertainty)
