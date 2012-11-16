#!/usr/bin/env/python

'''New basic definition of an MTK Point Source Class'''

import numpy as np
from copy import deepcopy
from nhlib.geo.point import Point
from nhlib.geo.mesh import Mesh
from nhlib.geo._utils import (spherical_to_cartesian,
                              cartesian_to_spherical)

class mtkPointSource(object):
    '''New class to describe the mtkPointsource object'''
    def __init__(self, identifier, name, tectonic_region, aspect_ratio, 
        input_geometry, upper_depth=0.0, lower_depth=None):
        '''Instantiate class with just the basic attributes
        :param identifier: 
            Integer ID code for the source
        :param name:
            Source Name (string)
        :param tectonic_region:
            Tectonic Region Type (String)
        :param aspect_ratio:
            Ratio of along-strike length to down-dip width (float)
        :param input_geometry:
            Point geometry of the source [Long, Lat]
        :param upper_depth:
            Upper seismogenic depth (km) (Non-negative Float) 
        :param lower_depth:
            Lower Seismogenic depth (km) (Non-negative Float)
        '''
        self.typology = 'Point'
        self.source_id = identifier
        self.name = name
        self.tectonic_region_type = tectonic_region
        self.aspect_ratio = aspect_ratio
        self.mfd = None
        self.msr = None
        if upper_depth < 0.0:
            raise ValueError('Upper Depth Must be Non Negative')
        self.upper_depth = upper_depth
        if lower_depth and (lower_depth < upper_depth):
            raise ValueError(
                'Lower Depth must be greater than or equal to upper depth')
        self.lower_depth = lower_depth
        self.nodal_plane_dist = None
        self.hypocentre_dist = None
        # Parse geometry
        self.geometry = Point(input_geometry[0], input_geometry[1],
                              upper_depth)

    
    def find_earthquakes_in_source(self, catalogue, distance, 
        upper_depth=None, lower_depth=None, **kwargs):
        '''Finds earthquakes within a distance of the point
        :param catalogue:
            Catalogue Dictionary (MTKCatalogue.data) 
        :geo_buffer:
            Distance from point to consider
        :upper_depth:
            Upper seismogenic depth (km)
        :lower_depth:
            Lower seismogenic depth (km)
        :**kwargs:
           Possible keys 'distance_type' - 'hypocentral' | '{epicentral}'
                         'geometry_type' - '{circle}' | 'square'
        '''
        self.catalogue = deepcopy(catalogue)
        if 'distance_type' in kwargs.keys():
            distance_type = kwargs['distance_type']
        else:
            distance_type = 'epicentral'
        geometry_type = 'circle'
        if not distance or distance < 0:
            distance = 0.0
        
        if 'geometry_type' in kwargs.keys():
            geometry_type = kwargs['geometry_type']
            if not geometry_type.lower() in ['square', 'circle']:
                raise ValueError('Unrecognised geometry type %s', 
                                 geometry_type)
        source_selector = Selector(self.catalogue.data)
        if geometry_type.lower() is 'square':
            self.catalogue.data, self.catalogue.number_earthquakes = \
                source_selector.select_square_centred_on_point(
                    self.geometry,
                    distance,
                    upper_depth=upper_depth,
                    lower_depth=lower_depth,
                    distance_type=distance_type)
        else:
             self.catalogue.data, self.catalogue.number_earthquakes = \
                source_selector.select_circular_distance_from_point(
                    self.geometry,
                    distance,
                    upper_depth=upper_depth,
                    lower_depth=lower_depth,
                    distance_type=distance_type)           
    
    
    #~ def find_earthquakes_in_source(self, catalogue, geo_buffer = 0.0,        
        #~ upper_depth = None, lower_depth = None, search_geometry=None):
        #~ '''Identifies the locations in the catalogue within a "geo_buffer"
        #~ distance (km) of the source'''
        #~ # Create mesh object from hypocentres
        #~ if catalogue['depth'] is not None:
            #~ hypocentres = Mesh(catalogue['longitude'], 
                               #~ catalogue['latitude'], 
                               #~ None)
        #~ else:
            #~ hypocentres = Mesh(catalogue['longitude'], 
                               #~ catalogue['latitude'], 
                               #~ catalogue['depth'])
        #~ # Find events within distance of point
        #~ if search_geometry == 'square':
            #~ # Geometry defined as a square cell
            #~ limits = _create_limits_of_square_cell(
                #~ self.geometry.longitude,
                #~ self.geometry.latitude,
                #~ geo_buffer)
            #~ valid_long = np.logical_and(
                #~ catalogue['longitude'] >= limits[3],
                #~ catalogue['longitude'] < limits[2])
            #~ valid_dist = np.logical_and(valid_long,
                #~ catalogue['latitude'] >= limits[1],
                #~ catalogue['latitude'] < limits[0])                
        #~ else:
            #~ valid_dist = self.geometry.closer_than(hypocentres, 
                                                   #~ geo_buffer)
        #~ if not upper_depth:
            #~ upper_depth = 0.0
        #~ if not lower_depth:
            #~ lower_depth = np.inf
        #~ valid_depth = np.logical_and(catalogue['depth'] >= upper_depth,
                                     #~ catalogue['depth'] < lower_depth)
        #~ valid_idx = np.logical_and(valid_dist, valid_depth,
                                   #~ catalogue['flagvector'])
        #~ self.catalogue = purge_catalogue(deepcopy(catalogue), valid_idx)
        
#~ def _create_limits_of_square_cell(long, lat, distance):
        #~ '''Function to define the north, south, west and east extent
        #~ of a square geographical cell centred on a point Long, Lat
        #~ '''
        #~ cart_loc = spherical_to_cartesian(long, lat, None)
        #~ west = cartesian_to_spherical(np.array([
            #~ cart_loc[0] - distance, cart_loc[1], cart_loc[2]]))
        #~ east = cartesian_to_spherical(np.array([
            #~ cart_loc[0] + distance, cart_loc[1], cart_loc[2]]))
        #~ north = cartesian_to_spherical(np.array([
            #~ cart_loc[0], cart_loc[1] + distance, cart_loc[2]]))
        #~ south = cartesian_to_spherical(np.array([
            #~ cart_loc[0], cart_loc[1] - distance, cart_loc[2]]))
        #~ return north, south, east, west