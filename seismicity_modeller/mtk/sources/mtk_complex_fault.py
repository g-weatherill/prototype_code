#!/usr/bin/env/python

'''New basic definition of an MTK Complex Fault Source Class'''

from copy import deepcopy
import numpy as np
from scientific.catalogue_utilities import purge_catalogue
from nhlib.geo.point import Point
from nhlib.geo.polygon import Polygon
from nhlib.geo.line import Line
from nhlib.geo.surface.complex_fault import ComplexFaultSurface
from nhlib.geo.mesh import Mesh, RectangularMesh

class mtkComplexFaultSource(object):
    '''New class to describe the mtk complex fault source object'''
    def __init__(self, identifier, name, tectonic_region, aspect_ratio, 
        edges, mesh_spacing = 1.0):
        '''Instantiate class with just the basic attributes
        :param identifier: 
            Integer ID code for the source
        :param name:
            Source Name (string)
        :param tectonic_region:
            Tectonic Region Type (String)
        :param aspect_ratio:
            Ratio of along-strike length to down-dip width (float)
        :param edges:
            List of edges([Long, Lat, Depth],[Long,Lat, Depth])
         :param mesh_spacing:
            Spacing of mesh (km) (Non-negative Float)       
        '''
        self.typology = 'ComplexFault'
        self.source_id = identifier
        self.name = name
        self.tectonic_region_type = tectonic_region
        self.aspect_ratio = aspect_ratio
        self.mfd = None
        self.msr = None
        
        number_edges = len(edges)
        self.upper_depth = np.inf
        self.lower_depth = -np.inf
        self.traces = []
        for iloc, edge in enumerate(edges):
            # Convert edge to line
            if np.min(edge[:, 2]) < self.upper_depth:
                self.upper_depth = np.copy(np.min(edge[:, 2]))
                self.lower_depth = np.copy(self.upper_depth)
            if np.max(edge[:, 2]) > self.lower_depth:
                self.lower_depth = np.copy(np.max(edge[:, 2]))
            self.traces.append(
                Line([Point(node[0], node[1], node[2]) 
                      for node in edge]))
        self.geometry = ComplexFaultSurface.from_fault_data(
            self.traces, mesh_spacing)
        self.catalogue = None
        
    def find_earthquakes_in_source(self, catalogue, distance, 
        upper_depth = None, lower_depth = None, **kwargs):
        '''Identifies earthquakes within a distance (km) of a fault
        :param catalogue:
            Input catalogue as instance of the MTKCatalogue class
        :param distance:
            Distance in km
        :param upper_depth:
            Upper seismogenic depth (km)
        :param lower_depth:
            Lower seismogenic depth (km)
        :param **kwargs:
            **kwargs can include key 'distance_type' 
              ('rupture' | 'joyner-boore')
        '''
        self.catalogue = deepcopy(catalogue)
        # If distance type is not specified or is set to anything
        # other than 'joyner-boore' then default to 'rupture'
        if 'distance_type' not in kwargs.keys():
            rupture_distance = True
        elif kwargs['distance_type'] is 'joyner-boore':
            rupture_distance = False
        else:
            rupture_distance = True
        
        source_selector = Select(self.catalogue.data)
        if rupture_distance:
            self.catalogue.data, self.catalogue.number_events = \
                source_selector.select_within_rupture_distance(
                    self.geometry,
                    distance,
                    upper_depth=upper_depth,
                    lower_depth=lower_depth)
        else:
            self.catalogue.data, self.catalogue.number_events = \
                source_selector.select_within_joyner_boore_distance(
                    self.geometry,
                    distance,
                    upper_depth=upper_depth,
                    lower_depth=lower_depth)

    #~ def find_earthquakes_in_source(self, catalogue, geo_buffer = 0.0, 
        #~ upper_depth = None, lower_depth = None):
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
            
        #~ rupture_dist = self.geometry.get_min_distance(hypocentres)
        #~ if not upper_depth:
            #~ upper_depth = 0.0
        #~ if not lower_depth:
            #~ lower_depth = np.inf
        #~ valid_depth = np.logical_and(catalogue['depth'] >= upper_depth,
                                     #~ catalogue['depth'] < lower_depth)
        #~ valid_idx = np.logical_and(rupture_dist <= geo_buffer, valid_depth,
                                   #~ catalogue['flagvector'])
        #~ self.catalogue = purge_catalogue(deepcopy(catalogue), valid_idx)
