#!/usr/bin/env/python

'''New basic definition of an MTK Simple Fault Source Class'''

#from math import fabs
from copy import deepcopy
import numpy as np
from scientific.catalogue_utilities import purge_catalogue
from nhlib.geo.point import Point
from nhlib.geo.polygon import Polygon
from nhlib.geo.line import Line
from nhlib.geo.surface.simple_fault import SimpleFaultSurface
from nhlib.geo.mesh import Mesh, RectangularMesh

class mtkSimpleFaultSource(object):
    '''New class to describe the mtk Simple fault source object'''
    def __init__(self, identifier, name, tectonic_region, aspect_ratio,
        fault_trace, dip, upper_depth, lower_depth, mesh_spacing = 1.0):
        '''Instantiate class with just the basic attributes
        :param identifier: 
            Integer ID code for the source
        :param name:
            Source Name (string)
        :param tectonic_region:
            Tectonic Region Type (String)
        :param aspect_ratio:
            Ratio of along-strike length to down-dip width (float)
        :param fault_trace:
            Surface trace of the fault [Long., Lat] as np.ndarray
        :param dip:
            Dip of fault (degrees) (Float between 0 and 90)
        :param upper_depth:
            Upper seismogenic depth (km) (Non-negative Float)
        :param lower_depth:
            Lower Seismogenic depth (km) (Non-negative Float)
        :param mesh_spacing:
            Spacing of mesh (km) (Non-negative Float)
        '''
        self.typology = 'SimpleFault'
        self.source_id = identifier
        self.name = name
        self.tectonic_region_type = tectonic_region
        self.aspect_ratio = aspect_ratio
        self.mfd = None
        self.msr = None
        if upper_depth < 0.0:
            raise ValueError('Upper Depth Must be Non Negative')   
        if lower_depth < 0.0 or (lower_depth < upper_depth):
            raise ValueError('Lower Depth must be Non-Negative and > Upper depth')
        self.trace = []
        if np.shape(fault_trace)[1] < 3:
            depth = np.zeros(np.shape(fault_trace)[0], dtype=float)
        else:
            depth = fault_trace[:, 2]
        # TODO if the trace contains only two points could use a planar surface!
        for iloc, node in enumerate(fault_trace):
            self.trace.append(Point(node[0], node[1], 
                depth[iloc]))
        self.trace = Line(self.trace)
        self.geometry = SimpleFaultSurface.from_fault_data(
            self.trace, upper_depth, lower_depth, dip, mesh_spacing)
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


