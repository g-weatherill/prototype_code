#!/usr/bin/env/python

'''Demo of the fault geometry class - with abstraction'''
import abc
import numpy as np
import math
from geo_utils import get_trace_length
#Eventually import modules from nhlib


class Geometry(object):
    '''Base class for geometry calculations'''
    __metaclass__ = abc.ABCMeta

    
    @abc.abstractmethod
    def get_fault_extent(self, geometry):
        '''Function to extract the trace (as numpy 2d array) and upper and 
        lower depths of the fault'''


    @abc.abstractmethod
    def get_width_from_geometry(self):
        '''Calculates the down-dip width of the fault using the geometry 
        parameters'''

    @abc.abstractmethod
    def get_basic_parameters(self, geometry):
        '''Function to extract the along-strike length, down-dip width and
        fault area'''

    @abc.abstractmethod
    def get_fault_surface(self):
        '''Calculates the fault surface - call library from nhlib (TODO)'''
        
class Simple(Geometry):
    '''Calculates the basic geometrical properties of the fault, assuming
    a simple fault typology'''
    def __init__(self):
        '''Initialise class'''
        self.upper_depth = None
        self.upper_depth_min = None
        self.upper_depth_max = None
        self.upper_depth_quality = None
        self.lower_depth = None
        self.lower_depth_min = None
        self.lower_depth_max = None
        self.lower_depth_quality = None
        self.dip = None
        self.dip_min = None
        self.dip_max = None
        self.dip_quality = None
        self.dip_direction = None
        self.trace = None
        self.length = None
        self.width = None
        self.area = None

    def _determine_upper_depths(self, depth_values):
        '''Returns the upper depth values, including pref, min, max and
        quality if provided as a list (len = 4)'''
        self.upper_depth = depth_values
        if (not isinstance(self.upper_depth, list)) or \
            (len(self.upper_depth) < 4):
            # As thhe preferred value of upper_depth is always listed first
            # need to instantiate list if no value or scalar value is given
            self.upper_depth = self.upper_depth
            self.upper_depth_min = None
            self.upper_depth_max = None
            self.upper_depth_quality = None
        elif len(self.upper_depth) == 4:
            self.upper_depth_min = self.upper_depth[1]
            self.upper_depth_max = self.upper_depth[2]
            self.upper_depth_quality = self.upper_depth[3]
            self.upper_depth = self.upper_depth[0]
    
    def _determine_lower_depths(self, depth_values):
        '''Returns the lower depth values as , including pref, min, max and
        quality if provided as a list (len = 4), or preferred'''
        self.lower_depth = depth_values
        if (not isinstance(self.lower_depth, list)) or \
            (len(self.lower_depth) < 4):
            self.lower_depth = self.lower_depth
            self.lower_depth_min = None
            self.lower_depth_max = None
            self.lower_depth_quality = None
        else:
            self.lower_depth_min = self.lower_depth[1]
            self.lower_depth_max = self.lower_depth[2]
            self.lower_depth_quality = self.lower_depth[3]
            self.lower_depth = self.lower_depth[0]

    def _determine_dip_parameters(self, dip_values):
        '''Returns dip, inclding preferred, min, max and quality if 
        provided as a 4-element list'''
        if (not isinstance(dip_values, list)) or (len(dip_values) < 4):
            self.dip = dip_values
            self.dip_min = None
            self.dip_max = None
            self.dip_quality = None
        else:
            self.dip = dip_values[0]
            self.dip_min = dip_values[1]
            self.dip_max = dip_values[2]
            self.dip_quality = dip_values[3]

    def get_fault_extent(self, geometry):
        '''Initialise'''
        self.trace = np.array(geometry['Fault_Trace'])
        #Re-arrange from list into 2-column array
        self.trace = np.column_stack(
            [self.trace[range(0, len(self.trace), 2)], 
             self.trace[range(1, len(self.trace), 2)]])
        self._determine_upper_depths(geometry['Upper_Depth'])
        self._determine_lower_depths(geometry['Lower_Depth'])
        self._determine_dip_parameters(geometry['dip'])            
        
    def get_width_from_geometry(self):
        '''Calculate Area'''
        self.width = (self.lower_depth - self.upper_depth) /\
            math.sin(self.dip * (math.pi / 180.))
    
    def get_basic_parameters(self, geometry):
        '''Defines the length, width and area of the fault'''
        self.length = \
            get_trace_length(self.trace[:, 0], self.trace[:, 1])[0][0]
        self.dip = geometry['dip']
        self.dip_direction = geometry['dip_direction']
        self.get_width_from_geometry()
        self.area = self.length * self.width

    def get_fault_surface(self):
        '''Defines the fault surface mesh'''
        #TODO - Call nhlib 
        print 'Not implemented yet - wait for nhlib'
