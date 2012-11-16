#!/usr/bin/env/python

'''Basic Set of Geometrical and Geographical Tools'''
import numpy as np
from pyproj import Geod, Proj
import shapely

def haversine(lon1, lat1, lon2, lat2, radians=False, earth_rad=6371.227):
    """
    Allows to calculate geographical distance
    using the haversine formula.

    :param lon1: longitude of the first set of locations
    :type lon1: numpy.ndarray
    :param lat1: latitude of the frist set of locations
    :type lat1: numpy.ndarray
    :param lon2: longitude of the second set of locations
    :type lon2: numpy.float64
    :param lat2: latitude of the second set of locations
    :type lat2: numpy.float64
    :keyword radians: states if locations are given in terms of radians
    :type radians: bool
    :keyword earth_rad: radius of the earth in km
    :type earth_rad: float
    :returns: geographical distance in km
    :rtype: numpy.ndarray
    """

    if radians == False:
        cfact = np.pi / 180.
        lon1 = cfact * lon1
        lat1 = cfact * lat1
        lon2 = cfact * lon2
        lat2 = cfact * lat2

    # Number of locations in each set of points
    if not np.shape(lon1):
        nlocs1 = 1
        lon1 = np.array([lon1])
        lat1 = np.array([lat1])
    else:
        nlocs1 = np.max(np.shape(lon1))
    if not np.shape(lon2):
        nlocs2 = 1
        lon2 = np.array([lon2])
        lat2 = np.array([lat2])
    else:
        nlocs2 = np.max(np.shape(lon2))
    # Pre-allocate array
    distance = np.zeros((nlocs1, nlocs2))
    i = 0
    while i < nlocs2:
        # Perform distance calculation
        dlat = lat1 - lat2[i]
        dlon = lon1 - lon2[i]
        idl_check = dlon > np.pi
        if np.sum(idl_check) > 0:
            dlon[idl_check] = 2. * np.pi - dlon[idl_check]
        aval = (np.sin(dlat / 2.) ** 2.) + (np.cos(lat1) * np.cos(lat2[i]) *
             (np.sin(dlon / 2.) ** 2.))
        distance[:, i] = (2. * earth_rad * np.arctan2(np.sqrt(aval),
                                                    np.sqrt(1 - aval))).T
        i += 1
    return distance
    
def get_trace_length(longitude, latitude):
    '''Get's the length of the fault trace'''
    num_points = np.shape(longitude)[0]
    if num_points != np.shape(latitude)[0]:
        raise ValueError('Longitudes and Latitudes have different length')
    
    length = 0.
    for iloc in range(1, num_points):
        length = length + haversine(longitude[iloc], latitude[iloc],
                                    longitude[iloc - 1], latitude[iloc - 1])
    return length

def point_in_polygon_across_dateline(long_point, lat_point, long_poly, 
    lat_poly):
    '''Point in polygon calculator - across international dateline'''
    if (max(long_poly) - min(long_poly)):
        return some_sht
        i
