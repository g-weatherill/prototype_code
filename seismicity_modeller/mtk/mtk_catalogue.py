#!/usr/bin/env/python
'''Prototype of a 'MTKCatalogue' class'''

import csv
import os
import numpy as np
from copy import deepcopy
from matplotlib.nxutils import points_inside_poly
from shapely import wkt
from shapely.geometry import asPolygon, asMultiPoint
from parsers.catalogue.csv_formats import CatalogueCSVParser
from scientific.selector import Selector
from scientific.declustering import Declustering
from scientific.new_completeness import Completeness
from scientific.new_recurrence import Recurrence
from scientific.catalogue_utilities import events_in_polygon, purge_catalogue
from nhlib.geo._utils import spherical_to_cartesian
from nhlib.geo.polygon import Polygon

TOTAL_ATTRIBUTE_LIST = ['eventID', 'Agency', 'Identifier', 'year', 'month', 
                        'day', 'hour', 'minute', 'second', 'timeError', 
                        'longitude', 'latitude', 'SemiMajor90', 'SemiMinor90', 
                        'ErrorStrike','depth','depthError','magnitude', 
                        'sigmaMagnitude', 'magnitudeType', 'focalMechanism', 
                        'validIndex']

FLOAT_ATTRIBUTE_LIST = ['second', 'timeError', 'longitude', 'latitude', 
                        'SemiMajor90', 'SemiMinor90', 'ErrorStrike', 'depth',
                        'depthError', 'magnitude', 'sigmaMagnitude']

INT_ATTRIBUTE_LIST = ['eventID','year', 'month', 'day', 'hour', 'minute']

STRING_ATTRIBUTE_LIST = ['Agency', 'magnitudeType']

file_parser = {'csv': CatalogueCSVParser()}



class MTKCatalogue(object):
    '''General Catalogue Class'''

    def __init__(self):
        '''Initilise the catalogue dictionary'''
        self.data = {}
        self.processes = {'declustering': None,
                          'completeness': None,
                          'recurrence': None,
                          'Poisson Tests': None}

        for attribute in TOTAL_ATTRIBUTE_LIST:
            if attribute in FLOAT_ATTRIBUTE_LIST:
                self.data[attribute] = np.array([], dtype=float)
            elif attribute in INT_ATTRIBUTE_LIST:
                self.data[attribute] = np.array([], dtype=int)
            else:
                self.data[attribute] = []
        self.data['xyz'] = None
        self.data['flag_vector'] = None
        self.number_earthquakes = 0
        self.default_completeness = None
    
    def read_catalogue(self, input_file, filetype):
        '''Reads the catalogue according to the file type'''
        # Check to see if file type supported
        if filetype in file_parser.keys():
            # Is supported 
            self.data = file_parser[filetype].parse_catalogue(self.data, 
                                                              input_file)
            self.number_earthquakes = len(self.data['eventID'])
            self.default_completeness = np.array([[
                np.min(self.data['year']), np.min(self.data['magnitude'])]])
        else:
            # Throw error
            raise ValueError('Catalogue file type not supported!')


    def write_catalogue(self, output_file, filetype):
        # Nothing here yet!
        raise AttributeError('Not implemented yet!')

    def decluster(self, config_file, purge=False):
        '''Master for declustering algorithm'''
        self.processes['declustering'] = Declustering()
        self.data['cluster_id'], self.data['flag_vector'] =  \
            self.processes['declustering'].run_declustering(self.data, 
                                                            config_file)
        if purge:
            self.data = purge_catalogue(self.data, self.data['flag_vector'])
            self.number_earthquakes = len(self.data['eventID'])
        

    def poisson_test(self, config_file):
        '''Poisson Chi^2 tests -  Nothing here yet!'''
        raise AttributeError('Not implemented yet!')


    def completeness_single(self, config_file, completeness_table=None,
        polygon=None, shallow_limit=None, deep_limit=None, 
        polygon_buffer=None):
        '''Implementes catalogue completeness analysis
        :param config_file:
            Dictionary containing the following configuration settings
            'algorithm': Choice of 'Stepp' | 'AssumedMFD'
            'delta_m': Magnitude window
            'delta_t': Time window
            'tolerance': Tolerance of Stepp method
            'increment_lock': Ensures completeness always reduces {True}

        :param polygon_set:
            List of instances of  nhlib.geo.polygon.Polygon class
        :param shallow_limit:
            Shallowest depth (km) of earthquakes considered in analysis
        :param deep_limit:
            Deepest depth (km) of earthquakes considered in analysis
        '''
        if polygon:
            # TODO Can this be adapted so that you can limit to a depth
            # range without needing a polygon
            self.processes['completeness'] = CatalogueCompleteness(
                self.data
                polygon,
                shallow_limit,
                deep_limit,
                polygon_buffer)
        else:
            self.process['completeness'] = CatalogueCompleteness(self.data)
        
        # If a completeness table is supplied then this overrides the default
        if isinstance(completeness_table, np.ndarray):
            self.default_completeness = completeness_table
            self.process['completeness'].completeness_table = \
                completeness_table
        else:
            self.process['completeness'].get_completeness_table(config_file)


    def recurrence(self, config_file):
        '''Calculates the recurrence parameters'''
        if isinstance(self.processes['completeness'], list):
            # Many different sources were considered - calculate recurrence
            # for each geometry
            self.processes['recurrence'] = []
            recurrence_model = []
            for comp_analysis in self.processes['completeness']:
                recurrence_analysis = Recurrence()
                recurrence_model.append(recurrence_analysis.run_recurrence(
                    comp_analysis.catalogue, 
                    config_file,
                    comp_analysis.completeness_table))
                self.processes['recurrence'].append(recurrence_analysis)
        else:
            self.processes['recurrence'] = Recurrence()
            recurrence_model = self.processes['recurrence'].run_recurrence(
                self.data, 
                config_file,
                self.default_completeness)
        return recurrence_model


    def maximum_magnitude(self, config_file, valid_config=False, **kwargs):
        '''Calculates the maximum magnitude using the instrumental seismicity
        '''
        output = {'Mmax':None, 'Mmax_Uncertainty':None}
        self.process['mmax'] = InstrumentalMmax()
        if not valid_config:
            self.process['mmax'].check_config(config_file)
            valid_config = True
        output['Mmax'], output['Mmax_Uncertainty'] = \
            self.process['mmax'].analyse(self.catalogue, config_file)
        return output, valid_config


#    def completeness(self, config_file, completeness_table=None,
#        polygon_set=None, shallow_limit=None, deep_limit=None, 
#        polygon_buffer=None):
#        '''Implementes catalogue completeness analysis
#        :param config_file:
#            Dictionary containing the following configuration settings
#            'algorithm': Choice of 'Stepp' | 'AssumedMFD'
#            'delta_m': Magnitude window
#            'delta_t': Time window
#            'tolerance': Tolerance of Stepp method
#            'increment_lock': Ensures completeness always reduces {True}
#
#        :param polygon_set:
#            List of instances of  nhlib.geo.polygon.Polygon class
#        :param shallow_limit:
#            Shallowest depth (km) of earthquakes considered in analysis
#        :param deep_limit:
#            Deepest depth (km) of earthquakes considered in analysis
#        '''
#        # If a completeness table is supplied then this overrides the default
#        if isinstance(completeness_table, np.ndarray):
#            has_table = True
#        else:
#            has_table = False
#        valid_config = False
#        if isinstance(polygon_set, list):
#            # Calculates separate completeness tables for each region
#            self.processes['completeness'] = []
#            # TODO This still seems really ugly - maybe need to rethink!
#
#            for polygon in polygon_set:
#                temp_completeness = CatalogueCompleteness(self.data,
#                                                          polygon,
#                                                          shallow_limit,
#                                                          deep_limit,
#                                                          polygon_buffer)
#                                
#                if has_table:
#                    # Completeness table is defined - therefore override
#                    temp_completeness.completeness_table = \
#                        self.default_completeness
#                else:
#                    temp_completeness.get_completeness_table(
#                        config_file,
#                        valid_config)
#                
#                self.processes['completeness'].append(temp_completeness)
#               
#                if not valid_config:
#                    valid_config = True
#
#        else:
#            model = CatalogueCompleteness(self.data)
#            if has_table:
#                model.completeness_table = self.default_completeness
#                #self.process['completeness'] = self.default_completeness
#            else:
#                model.get_completeness_table(config_file, valid_config)
#            if not valid_config:
#                valid_config = True
#            self.processes['completeness'] = [model]


    def _hypocentres_to_xyz(self):
        '''Converts the catalogue hypocentres into a cartesian coordinate
        system'''
        self.data['xyz'] = spherical_to_cartesian(
            self.data['longitude'], self.data['latitude'], 
            self.data['depth'])



def _xyz_to_wkt(xyz):
    '''Converts the xy of the xyz into a wkt MultiPoint format'''
    return asMultiPoint(xyz[:, :-1]).wkt

    
def _polygon_wkt_to_xyz(polygon, as_wkt=False):
    '''Converts a polygon in 'POLYGON' wkt format to xyz. If as_wkt is 
    True, will return the xyz output back as a wkt format, otherwise
    it will return an array'''
    polygon = wkt.loads(polygon)
    polygon = polygon.boundary.xy
    number_nodes = len(polygon[0])
    xyz =  spherical_to_cartesian(polygon[0], polygon[1], 
        np.zeros(number_nodes, dtype=float))
    if as_wkt:
        xyz = asPolygon(xyz[:, :-1]).wkt
    return xyz



class CatalogueCompleteness(object):
    '''A separate object for storing completeness information'''
    def __init__(self, input_catalogue, source_polygon=None, 
        shallow_depth=None, deep_depth=None, polygon_buffer=None):
        '''Initialise completeness class'''
        self.source_geometry = source_polygon
        self.completeness_table = None
        self.upper_depth = shallow_depth
        self.lower_depth = deep_depth
        # Cut catalogue to source geometry
        if isinstance(source_polygon, Polygon):
            catalogue_selector = Selector(input_catalogue)
            self.catalogue = catalogue_selector.select_within_polygon(
                source_polygon,
                distance = polygon_buffer,
                upper_depth = self.upper_depth,
                lower_depth = self.lower_depth)[0]
        else:
            self.catalogue = deepcopy(input_catalogue)

        self.model = Completeness() 
        self.completeness_table = None
        self.valid_config = False


    def get_completeness_table(self, config, is_valid=False):
        '''Applies completeness algorithms'''
        self.valid_config = is_valid
        self.completeness_table = self.model.analysis(self.catalogue, 
                                                      config,
                                                      self.valid_config)
        

#    def select_earthquakes_in_polygon(self, polygon, upper_depth=None,
#        lower_depth = None):
#        '''Creates a subcatalogue containing only the events inside the polygon
#        :param source:
#            polygon is a numpy.ndarray of [Longitude, Latitude]
#        :param upper_depth:
#            Shallowest depth (km) of earthquakes considered in analysis
#        :param lower_depth:
#            Deepest depth (km) of earthquakes considered in analysis
#        :returns: 
#            "sub_catalogue" dictionary of the same format as data, which
#            corresponds to only the events inside the polygon
#        '''
#        # To avoid mis-allocations from application on a sphere convert
#        # Long. and Lat to cartesian
#        if not self.data['xyz']:
#            self._hypocentres_to_xyz()
#        number_nodes = np.shape(polygon)[0]
#        polygon = spherical_to_cartesian(polygon[:, 0], polygon[:, 1],
#                                         np.zeros(number_nodes, dtype=float))
#        selected_id = np.zeros(len(self.data['eventID']), dtype=bool)
#        # Set to True all events inside the depth range, if specified
#        if not upper_depth:
#            upper_depth = 0.0 
#        selected_id[self.data['xyz'][:, 2] >= upper_depth] = True
#        if not lower_depth:
#            lower_depth = np.inf
#        
#        selected_id[self.data['xyz'][:, 2] <= lower_depth] = True
#        valid_locations = np.where(selected_id)[0]
#        
#        in_polygon_id = points_inside_poly(
#            self.data['xyz'][valid_locations, :2], polygon[:, :2])
#        selected_id[np.where(selected_id)[0]] = in_polygon_id
#        
#        return purge_catalogue(self.data, selected_id.astype(int))


class MomentTensor(object):
    '''Creates a general class describing the Moment Tensor of the Event'''
    def __init__(self):
        '''Initialise'''
        self.hypocentre = {}
        self.moment = None
        self.m_w = None
        self.nodal_plane1 = {'strike': None, 'dip': None, 'rake': None}
        self.nodal_plane2 = {'strike': None, 'dip': None, 'rake': None}
        self.npt_axes = {'n_plunge': None, 'n_azimuth': None, 
                         'p_plunge': None, 'p_azimuth': None,
                         't_plunge': None, 't_azimuth': None}
        self.moment_tensor = []
        self.moment_tensor_error = []
