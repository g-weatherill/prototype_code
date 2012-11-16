#!/usr/bin/env/python

'''New basic definition of an MTK Point Source Class'''

#from math import fabs
import numpy as np
from nhlib.source.geo.point import Point
from nhlib.source.geo.polygon import Polygon
from nhlib.source.geo.line import Line
from nhlib.geo.surface.simple_fault import SimpleFaultSurface
from nhlib.geo.surface.complex_fault import ComplexFaultSurface

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
            Point geometry of the source [Long, Lat]
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
        
        for iloc, node in enumerate(fault_trace):
            self.trace.append(Point(node[0], node[1], 
                depth[iloc]))
        self.trace = Line(self.trace)
        self.geometry = SimpleFaultSurface.from_fault_data(
            self.trace, self.upper_depth, self.lower_depth,
            self.dip, mesh_spacing)            
            
class mtkComplexFault(object):
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
