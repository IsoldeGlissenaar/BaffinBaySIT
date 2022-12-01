# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 09:51:02 2020

@author: IsoldeGlissenaar
"""

from scipy.interpolate import LinearNDInterpolator
from scipy import spatial
import numpy as np
from func_reproj import reproject

def in_hull(p, hull):
    """
    Test if points in 'p' are in 'hull'

    'p' should be a 'NxK' coordinates of 'N' points in 'K' dimensions
    'hull' is either a scipy.spatial.Delaunay object or the 'MxK' array of the 
    coordinates of 'M' points in 'K'dimensions for which Delaunay triangulation
    will be computed
    By: Juh_
    """
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0


def scatteredInterpolant_CS(field, data):
    # Function scatteredInterpolant_CS(field, data)
    # Input: field, data
    # Output: interpolant
    #
    # interpolates from field to data grid. 
    #
    # given data and coordinates in field[n,3] (lat, lon, values) and
    # coordinates in data, interpolates the values in field to the coordinates
    # in data.
    # Used to interpolate snow depth in func_snowdepth.py
    
    #Find nearest distance between satellite altimetry values and 
    #snow depth data points
    coord = np.transpose(np.array([data.lat.values, data.lon.values]))
    inter = in_hull(coord, field[:,0:2])
    
    interpolant = np.zeros(len(data))
    
    #Liner interpolation inside convex hull
    int_sd = LinearNDInterpolator(points=list(zip(field[:,0], field[:,1])), values=field[:,2])   
    interpolant[inter]= int_sd(coord[inter,0],coord[inter,1])
    
    #Nearest neighbour interpolation outside of convex hull
    tree = spatial.KDTree(field[:,0:2])
    dist, ids = tree.query(coord[~inter])
    interpolant[~inter] = field[ids,2]

    return interpolant
