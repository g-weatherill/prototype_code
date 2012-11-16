#!/usr/bin/env/python

'''New basic definition of an MTK Area Source Class'''
from math import fabs
from copy import deepcopy
import numpy as np
from nhlib.geo.point import Point
from nhlib.geo.polygon import Polygon
from nhlib.geo._utils import spherical_to_cartesian
from scientific.selector import Selector
from scientific.catalogue_utilities import (events_in_polygon,
                                                purge_catalogue)

class mtkAreaSource(object):
    '''New class to describe the mtkAreaSource object'''
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
            Input geometry as numpy array of [Long., Lat]
        :param upper_depth:
            Upper seismogenic depth (km) (Non-negative Float)
        :param lower_depth:
            Lower Seismogenic depth (km) (Non-negative Float)            
        '''
        self.typology = 'Area'
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
        # Parse Geometry
        self.geometry = []
        for row in input_geometry:
            self.geometry.append(Point(row[0], row[1], 
                                 self.upper_depth))
        self.geometry = Polygon(self.geometry)
        self.catalogue = None
        self.number_events = None
        
        
    def find_earthquakes_in_source(self, input_catalogue, distance,
        upper_depth = None, lower_depth = None, **kwargs):
        '''From a given earthquake catalogue in the form of an MTKCatalogue
        class, retreive only the events inside the source
        :param catalogue:
            Instance of the MTKCatalogue class
        :param geo_buffer:
            Geographical buffer (km) to extend polygon selection
        :param upper_depth:
            Upper seismogenic depth (km)
        :param lower_depth:
            Lower seismogenic depth (km)
        '''
           
        self.catalogue = deepcopy(catalogue)
        source_selector = Selector(catalogue.data)
        self.catalogue.data, self.catalogue.number_events = \
            source_selector.select_within_polygon(
                self.geometry,
                distance=distance,
                upper_depth=upper_depth,
                lower_depth=lower_depth)
        if number_selected < 5:
            raise Warning('Source %i (%s) has fewer than 5 events' 
                          % self.source_id, self.name)
        

         
        
    #~ def find_earthquakes_in_source(self, catalogue, geo_buffer = 0.0, 
        #~ upper_depth = None, lower_depth = None):
        #~ '''For a given input catalogue, selects the earthquakes inside the 
        #~ area polygon'''
        #~ if (isinstance(geo_buffer, float) and 
            #~ (fabs(geo_buffer) > small_number)):
            #~ zone_polygon = self.polygon.dilate(geo_buffer)
        #~ else:
            #~ zone_polygon = self.polygon
        #~ # Convert spherical to cartesian
        #~ number_verts = len(zone_polygon.lons)
        #~ zone_polygon = spherical_to_cartesian(zone_polygon.lons,
            #~ zone_polygon.lats, np.zeros(number_verts, dtype=float))
        #~ zone_polygon = np.column_stack([zone_polygon[:, 0], 
                                        #~ zone_polygon[:, 1]])
        
        #~ if isinstance(catalogue['flagvector'], np.ndarray):
            #~ valid_flag = catalogue['flagvector']
        #~ else:
            #~ valid_flag = None

        #~ if not (isinstance(catalogue['xyz'], np.ndarray) and 
            #~ len(catalogue['xyz'] > 0)):
            #~ catalogue['xyz'] = spherical_to_cartesian(catalogue['longitude'],
                                                      #~ catalogue['latitude'],
                                                      #~ catalogue['depth'])

        #~ catalogue['flagvector'] = events_in_polygon(catalogue['xyz'],
                                                    #~ zone_polygon,
                                                    #~ valid_flag,
                                                    #~ upper_depth,
                                                    #~ lower_depth)
        #~ self.catalogue = purge_catalogue(deepcopy(catalogue), 
                                         #~ catalogue['flagvector'])
