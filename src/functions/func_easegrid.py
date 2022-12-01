# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:45:37 2020

@author: 
"""

import numpy as np

def easegrid(iopt, alat, alon, ascale):
    # Function easegrid(iopts, lat, lon, ascale)
    # Input: iopts, lat, lon, ascale
    # Output: thelon, thelat
    #
    # computes the forward 'ease' grid transform
    #
    # given a lat, lon (alat, alon) and the scale (ascale) the image
    # transformation coordinates (thelon, thelat) are computed
    # using the 'ease grid' (version 1.0) transformation given in fortran
    # source code supplied by nsidc.
    #
    # the radius of the earth used in this projection is imbedded into
    # ascale while the pixel dimension in km is imbedded in bscale
    # the base values are: radius earth = 6371.228 km
    #                       pixel dimen = 25.067525 km
    # then, bscale = base_pixel_dimen
    #       ascale = radius_earth/base_pixel_dimen
    #
    # iopt is ease type: iopth=11=north, iopt=12=south, iopt=13=cylindrical
    
    pi2 = 1.57079633
    dtr = pi2/90.0
    
    if iopt==11:        #ease grid north
        thelon = ascale*np.sin(alon*dtr)*np.sin(dtr*(45-0.5*alat))
        thelat = -ascale*np.cos(alon*dtr)*np.sin(dtr*(45-0.5*alat))
    elif iopt==12:      #ease grid south
        thelon = ascale*np.sin(alon*dtr)*np.cos(dtr*(45-0.5*alat))
        thelat = ascale*np.cos(alon*dtr)*np.cos(dtr*(45-0.5*alat))
    elif iopt==13:      #ease cylindrical
        thelon = ascale*pi2*alon*np.cos(30*dtr)/90
        thelat = ascale*np.sin(alat*dtr)/np.cos(30*dtr)
    return thelon, thelat
        