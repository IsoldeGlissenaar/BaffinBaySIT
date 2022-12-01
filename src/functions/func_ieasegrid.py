# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 11:09:43 2020

@author:
"""
import numpy as np

def ieasegrid(iopt, thelat, thelon, ascale):
    # Function ieasegrid(iopt,thelon,thelat,ascale)
    # Input: iopt, thelon, thelat, ascale
    # Output: alon, alat
    #
	# computes the inverse "ease" grid transform
    #
 	# given the image transformation coordinates (thelon,thelat) and
	# the scale (ascale) the corresponding lon,lat (alon,alat) is computed
	# using the "ease grid" (version 1.0) transformation given in fortran
	# source code supplied by NSIDC
	# iopt is ease type: iopt=11=north, iopt=12=south, iopt=13=cylindrical
    #
	# the radius of the earth used in this projection is imbedded into
	# ascale while the pixel dimension in km is imbedded in bscale
	# the base values are: radius earth= 6371.228 km
	#		               pixel dimen =25.067525 km
	# then, bscale = base_pixel_dimen
	#       ascale = radius_earth/base_pixel_dimen

    pi2 = 1.57079633
    dtr = pi2/90.0
    x1 = thelon
    y1 = thelat
    temp = np.zeros((thelon.shape))
    alat = np.zeros((thelat.shape))
    alon = np.zeros((thelon.shape))
    if iopt==11:    # ease grid north
        alon = np.arctan2(x1,-y1)/dtr
        ind = np.where(np.abs(np.sin(dtr*alon)) > np.abs(np.cos(alon*dtr)))
        temp[ind] = (x1[ind]/np.sin(alon[ind]*dtr))/ascale
        ind = np.where(np.abs(np.sin(dtr*alon)) <= np.abs(np.cos(alon*dtr)))
        temp[ind] = (-y1[ind]/np.cos(alon[ind]*dtr))/ascale
        alat[np.abs(temp)<=1] = 90.0-2*np.arcsin(temp[np.abs(temp) <= 1])/dtr
        alat[np.abs(temp)>1]  = 90.0*np.sign(temp[np.abs(temp)>1])
    elif iopt==12:  # ease grid south
        alon = np.arctan2(x1,y1)/dtr
        alat = np.zeros((alon.shape))
        ind = (np.abs(np.cos(alon*dtr)) > np.abs(np.sin(alon*dtr)))
        temp[ind] = (y1[ind]/np.cos(alon[ind]*dtr))/ascale
        alat[np.abs(temp)<=1] = 90.0-2*np.arccos(temp[np.abs(ind)<=1])/dtr
        alat[np.abs(temp)>1]  = 90.0*np.sign(temp[np.abs(temp)>1])
    elif iopt==13:  # ease cylindrical
        alon = ((x1/ascale)/np.cos(30*dtr))*90/pi2
        temp = (y1*np.cos(30*dtr))/ascale
        alat = np.zeros((alon.shape))
        alat[np.abs(temp)<=1] = np.arcsin(temp[np.abs(temp)<=1])/dtr
        alat[np.abs(temp)>1]  = 90*np.sign(temp[np.abs(temp)>1])
        
    return alon, alat

