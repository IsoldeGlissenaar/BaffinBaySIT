# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 11:48:27 2020

@author: zq19140
"""

def reproject(x,y):
    # Function reproject(x,y)
    # Input: x, y
    # Output: xnew, ynew
    #
    # reprojects coordinates to another coordinate system
    #
    # given x,y coordinates in the ESPG:4326 coordinate system (WGS84), the
    # coordinates a converted to the EPSG:3413 coordinate system (NSIDC Sea 
    # Ice Polar Stereographic North) using pyproj functions.  
    
    from pyproj import Proj, transform
    inProj = Proj('epsg:4326')
    outProj = Proj('epsg:3413')
    x1,y1 = x, y
    xnew,ynew = transform(inProj,outProj,x1,y1)
    return xnew, ynew


def backproject(x,y):
    # Function backproject(x,y)
    # Input: x, y
    # Output: xnew, ynew
    #
    # reprojects coordinates to another coordinate system
    #
    # given x,y coordinates in the EPSG:3413 coordinate system (NSIDC Sea 
    # Ice Polar Stereographic North), the coordinates a converted to the 
    # ESPG:4326 coordinate system (WGS84) using pyproj functions.  
    
    from pyproj import Proj, transform
    inProj = Proj('epsg:3413')
    outProj = Proj('epsg:4326')
    x1,y1 = x, y
    xnew, ynew = transform(inProj,outProj,x1,y1)
    return xnew, ynew
