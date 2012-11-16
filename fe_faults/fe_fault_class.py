#!/usr/bin/env/python
'''1st draft of a comprehensive object representation of a Faulted Earth fault
'''
from datetime import datetime
from math import fabs
import numpy as np
import uncertainty as unc
from scipy.stats import truncnorm
from nhlib.geo.point import Point
from nhlib.geo.line import Line
from nhlib.geo.surface.simple_fault import SimpleFaultSurface
from nhlib.geo.surface.complex_fault import ComplexFaultSurface
from sources.mtk_simple_fault import SimpleFaultSource


# Default mesh spacing = 1 km
default_spacing = 1.0

uncertainty_map = {'Unknown': unc.Uncertainty()
                   'Uniform': unc.Uniform(),
                   'Gaussian': unc.Gaussian(),
                   'TruncatedGaussian': unc.TruncatedGaussian()
                   }


def generate_discrete_normal_probability(mean, sigma, upper, lower, step):
    '''Generates a discrete list of probabilities from a Gaussian distribution
    '''
    upper_std = (upper - mean) / sigma
    lower_std = (lower - mean) / sigma
    epsilon = np.arange(lower_std, upper_std, step / 2.)
    probability = truncnorm.cdf(epsilon + (step / 2.), lower_std, upper_std) -\
        truncnorm.cdf(epsilon - (step / 2.), lower_std, upper_std)
    assert(fabs(np.sum(probability) - 1.0) < 1E-5) 
    return mean + (epsilon * sigma), probability



class NeotectonicObject(object):
    '''Class containing attributes and methods common to all definitions of
    a neotectonic fault - i.e. neotectonic section, neotectonic fault, 
    neotectonic fault source'''
    def __init__(self, identifier, name, compiled_by='Unknown', 
        contributed_by='Unknown', last_updated=datetime.utcnow()):
        '''Instantiate basic fault class'''
        self.identifier = identifier
        self.typology = 'simple'
        self.name = name
        #self.is_active = True
        #self.is_episodic = False
        self.length = None
        self.geometry = []
        self.upper_depth = None
        self.lower_depth = None
        self.dip = None
        self.dip_direction = None
        self.slip_rate = None
        self.slip_type = None
        self.aseismic = None
        self.displacement = None
        self.recurrence_interval = None
        self.data_completeness = None
        self.age_last_movement = None
        self.compiler = compiled_by
        self.contributer = contributed_by
        self.last_update = last_updated

    def render_geometry(self, trace, **kwargs):
        '''
        Renders the geometry of the section from the trace
        :param trace: 
            If a simple fault it is an instance of the nhlib.geo.line.Line
            class, for a complex fault it is a list of nhlib.geo.line.Line 
            classes
        :param **kwargs:
        '''
        if (not 'mesh_spacing' in kwargs.keys()) or not kwargs['mesh_spacing']:
            kwargs['mesh_spacing'] = 1.0

        if isinstance(trace, list):
            #Complex fault
            self.typology = 'complex'
            self.geometry = ComplexFaultSurface.from_fault_data(
                trace,
                kwargs['mesh_spacing'])
            self.dip = self.geometry.get_dip()
            self.upper_depth = self.geometry.get_top_edge_depth()
            self.lower_depth = np.max(self.geometry.mesh.depths)
        else:
            self.typology = 'simple'
            self.dip = kwargs['dip']
            self.upper_depth = kwargs['upper_seismogenic_depth']
            self.lower_depth = kwargs['lower_seismogenic_depth']
            self.geometry = SimpleFaultSurface.from_fault_data(
                trace,
                self.dip,
                self.upper_depth,
                self.lower_depth,
                kwargs['mesh_spacing'])


class NeotectonicSection(NeotectonicObject):
    '''Neotectonic Section'''
    def __init__(self, identifier, name, compiled_by='Unknown', 
        contributed_by='Unknown', last_updated=datetime.utcnow()):
        '''Instantiate class'''
        super(NeotectonicSection, self).__init__(identifier, 
                                                 name,
                                                 compiled_by,
                                                 contributed_by,
                                                 last_updated)
        self.vertical_slip_rate = None
        self.horizontal_slip_rate = None
        self.net_slip_rate = None
        self.downthrown_side = None
        self.is_active = True
        self.is_episondic = False


class NeotectonicFault(NeotectonicObject):
    '''Neotectonic Fault'''
    def __init__(self, identifier, name, sections, compiled_by='Unknown',
        contributed_by='Unknown', last_updated=datetime.utcnow()):
        '''Instantiate class
        :param sections:
            List of instances of NeotectonicSection class
        '''
        super(NeotectonicFault, self).__init__(identifier, 
                                               name,
                                               compiled_by,
                                               contributed_by,
                                               last_updated)
        self.number_sections = len(sections)
        self.active_sections = []
        for section in sections:
            self.active_sections.append(section)


class FaultSource(NeotectonicObject):
    '''Fault source object'''
    
    def __init__(self, identifier, data, compiled_by='Unknown',
        contributed_by='Unknown', last_updated=datetime.utcnow()):
        ''''''
        super(FaultSource, self).__init__(identifier,
                                          data['name'],
                                          compiled_by,
                                          contributed_by,
                                          last_updated)
        
        self.tectonic_region = data['tectonic_region']
        
        # Render the geometry
        self.geometry = self.render_geometry(
            data['trace'],
            dip=data['dip'],
            upper_seismogenic_depth=data['upper_seismogenic_depth'],
            lower_seismogenic_depth=data['lower_seismogenic_depth'])
                                             







class FaultTrace(object):
    '''Class to describe a fault trace'''
    def __init__(self, trace_id, locations,  method='Unknown', name='Unknown', 
                 scale=None, notes=None):
       '''Instantiate class parsing the geometry to a line trace'''
       self.identifier = trace_id
       self.method = method
       self.name = name
       self.scale = scale
       self.notes = notes
       self.trace = Line([Point(row[0],row[1]) for row in locations])


class Slip(object):
    '''Defines a class of slip parameters'''




#class FaultedEarthFault(object):
#    '''Master class representation of a Faulted Earth Fault'''
#    def __init__(self, identifier, name, compiled_by='Unknown', 
#        contributed_by='Unkown', last_updated=datetime.utcnow()):
#        '''Instantiate basic fault class'''
#       self.identifier = identifier
#        self.name = name
#        #self.is_active = True
#        #self.is_episodic = False
#        self.length = None
#        self.geometry = []
#        self.slip = None
#        self.aseismic = None
#        self.displacement = None
#        self.recurrence_int = None
#        self.data_completeness = None
#        self.age_last_movement = None
#        self.compiler = compiled_by
#        self.contributer = contributed_by
#        self.last_update = last_updated


#class FaultGeometry(object)
#
#    def __init__(self, fault_trace, upper_depth, lower_depth, dip,
#                    mesh_spacing=default_spacing, dip_direction=None):
#        '''
#        Instantiates the fault geometry class (all distances in km)
#        :param fault_trace: Instance of nhlib.geo.Line class
#        :param upper_depth: Upper Seismogenic Depth (km)
#        :param lower_depth: Lower Seismogenic Depth (km)
#        :param dip: Dip (degrees)
#
#        ''' 
#        self.trace = fault_trace.trace
#        self.upper_seismogenic_depth = upper_depth
#        self.lower_seismogenic_depth = lower_depth
#        self.dip = dip
#        self.dip_direction = dip_direction
#        self.length = fault_trace.get_length()
#        self.width = None
#        self.area = None
#        self.surface = SimpleFaultSurface.from_fault_data(
#            self.trace, 
#            self.upper_seismogenic_depth, 
#            self.lower_seismogenic_depth,
#            self.dip,
#            mesh_spacing)
#        self.width = self.surface.get_width()
#        self.strike = self.surface.get_strike()
#        if not self.dip_direction:
#            self.dip_direction = self.strike + 90.
#            if self.dip_direction  > 360.:
#                self.dip_direction = self.dip_direction - 360.
