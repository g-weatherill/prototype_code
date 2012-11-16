#!/usr/bin/env/python

'''Class to implement set of functionalities for selecting events from
and earthquake catalogue'''

import numpy as np
from copy import deepcopy
from matplotlib.nxutils import points_inside_poly
from nhlib.geo.geodetic import distance
from nhlib.geo.point import Point
from nhlib.geo._utils import spherical_to_cartesian
#from nhlib.geo.line import Line
from nhlib.geo.polygon import Polygon
from nhlib.geo.mesh import Mesh
from scientific.catalogue_utilities import(decimal_time, purge_catalogue)


def _check_depth_limits(input_dict):
    '''Returns the default upper and lower depth values if not in dictionary
    :param input_dict:
        Dictionary corresponding to the kwargs dictionary of calling function
    :returns:
        'upper_depth': Upper seismogenic depth (float)
        'lower_depth': Lower seismogenic depth (float)
    '''
    if ('upper_depth' in input_dict.keys()) and input_dict['upper_depth']:
        if input_dict['upper_depth'] < 0.:
            raise ValueError('Upper seimogenic depth must be positive')
        else:
            upper_depth = input_dict['upper_depth']
    else:
        upper_depth = 0.0

    if ('lower_depth' in input_dict.keys()) and input_dict['lower_depth']:
        if input_dict['lower_depth'] < upper_depth:
            raise ValueError('Lower depth must take a greater value than'
                             'upper depth!')
        else:
            lower_depth = input_dict['lower_depth']
    else:
        lower_depth = np.inf
    return upper_depth, lower_depth

class Selector(object):
    '''Selector class'''
    def __init__(self, catalogue):
        '''Instantiate
        :param catalogue:
            Instance of MTKCatalogue Class
        '''
        self.catalogue = deepcopy(catalogue)
        self.time_value = decimal_time(catalogue['year'], 
                                       catalogue['month'],
                                       catalogue['day'], 
                                       catalogue['hour'], 
                                       catalogue['minute'],
                                       catalogue['second'])

        self.catalogue_mesh = Mesh(catalogue['longitude'],
                                   catalogue['latitude'],
                                   catalogue['depth'])
        if not isinstance(catalogue['xyz'], np.ndarray):
            self.catalogue['xyz'] = spherical_to_cartesian(
                self.catalogue['longitude'],
                self.catalogue['latitude'],
                self.catalogue['depth'])
        self.number_events = len(catalogue['eventID'])

    def select_within_polygon(self, polygon, distance=None, **kwargs):
        '''Select earthquakes within polygon
        :param point:
            Centre point as instance of nhlib.geo.polygon.Polygon class
        :param distance:
            Buffer distance (km) (can take negative values)
        :returns:
            Selected catalogue (as dictionary) and number of selected events
        
        '''
        if distance:
            zone_polygon = polygon.dilate(distance)
        else:
            zone_polygon = polygon

        zone_polygon = spherical_to_cartesian(
            zone_polygon.lons,
            zone_polygon.lats, 
            np.zeros(len(zone_polygon.lons), dtype=float))
        
        # Initially all points are invalid
        valid_id = np.zeros(self.number_events, dtype=bool)
        
        # Make valid all events inside depth range
        upper_depth, lower_depth = _check_depth_limits(kwargs)
        valid_depth = np.logical_and(self.catalogue['depth'] >= upper_depth,
                                     self.catalogue['depth'] < lower_depth)
        # Events outside polygon returned to invalid assignment
        valid_id[valid_depth] = points_inside_poly(
            self.catalogue['xyz'][valid_depth, :2],
            zone_polygon[:, :2])

        number_selected = np.sum(valid_id)
        valid_id = np.logical_not(valid_id)
        return purge_catalogue(self.catalogue, valid_id.astype(int)), \
            number_selected
                            

    def select_circular_distance_from_point(self, point, distance, **kwargs):
        '''Select earthquakes within a distance from a Point
        :param point:
            Centre point as instance of nhlib.geo.point.Point class
        :param distance:
            Distance (km)
        :returns:
            Selected catalogue (as dictionary) and number of selected events
        '''
        if kwargs['distance_type'] is 'epicentral':
            locations = Mesh(catalogue['longitude'], catalogue['latitude'],
                np.zeros(len(catalogue['longitude']), dtype=float))
            point = Point(point.longitude, point.latitude, 0.0)
        else:
            locations = self.catalogue_mesh

        is_close = (point.closer_than(locations, distance))
        number_selected = np.sum(is_close)
        is_close = np.logical_not(is_close)
        return purge_catalogue(catalogue, is_close.astype(int)), \
            number_selected
        
    def select_square_centred_on_point(self, point, distance, **kwargs):
        '''Select earthquakes from within a square centered on a point
        :param point:
            Centre point as instance of nhlib.geo.point.Point class
        :param distance:
            Distance (km)
        :returns:
            Selected catalogue (as dictionary) and number of selected events
        '''
        point_surface = Point(point.longitude, point.latitude, 0.)
        north_point = point_surface.point_at(distance, 0., 0.)
        east_point = point_surface.point_at(distance, 0., 90.)
        south_point = point_surface.point_at(distance, 0., 180.)
        west_point = point_surface.point_at(distance, 0., 270.)
        is_long = np.logical_and(
            catalogue['longitude'] >= west_point.longitude,
            catalogue['longitude'] < east_point.longitude)
        is_surface = np.logical_and(
            is_long,
            catalogue['latitude'] >= south_point.latitude,
            catalogue['latitude'] < north_point.latitude)
        
        upper_depth, lower_depth = _check_depth_limits(kwargs)
        is_valid = np.logical_and(
            is_surface,
            catalogue['depth'] >= upper_depth,
            catalogue['depth'] < lower_depth)
        number_selected = np.sum(is_valid)
        is_valid = np.logical_not(is_valid)
        return purge_catalogue(catalogue, is_valid.astype(int)), \
            number_selected


    def select_within_joyner_boore_distance(self, surface, distance, **kwargs):
        '''Select events within a Joyner-Boore distance of a fault
        :param surface:
            Fault surface as instance of nhlib.geo.surface.base.BaseSurface
        :param distance:
            Rupture distance (km)
        :returns:
            Selected catalogue (as dictionary) and number of selected events
        '''
        
        upper_depth, lower_depth = _check_depth_limits(kwargs)
        
        rjb = surface.get_joyner_boore_distance(self.catalogue_mesh)
        is_valid = np.logical_and(rjb <= distance,
                                  catalogue['depth'] >= upper_depth,
                                  catalogue['depth'] < lower_depth)

        number_selected = np.sum(is_valid)
        is_valid = np.logical_not(is_valid)
        return purge_catalogue(catalogue, is_valid.astype(int)), \
            number_selected
        

    def select_within_rupture_distance(self, surface, distance,  **kwargs):
        '''Select events within a rupture distance from a fault surface
        :param surface:
            Fault surface as instance of nhlib.geo.surface.base.BaseSurface
        :param distance:
            Rupture distance (km)
        :returns:
            Selected catalogue (as dictionary) and number of selected events
        '''
        # Check for upper and lower depths 
        upper_depth, lower_depth = _check_depth_limits(kwargs)
        
        rrupt = surface.get_min_distance(self.catalogue_mesh)
        is_valid = np.logical_and(rrupt <= distance,
                                  catalogue['depth'] >= upper_depth,
                                  catalogue['depth'] < lower_depth)

        number_selected = np.sum(is_valid)
        is_valid = np.logical_not(is_valid)
        return purge_catalogue(catalogue, is_valid.astype(int)), \
            number_selected

    def select_within_time_period(self, start_time=None, end_time=None):
        '''Select earthquakes occurring within a given time period
        :param start_time:
            Earliest time (as datetime.datetime object)
        :param end_time:
            Latest time (as datetime.datetime object)
        :returns:
            Catalogue (Dict) of earthquakes inside time_interval
        '''
        if not(start_date):
            if not(end_date):
                # No times input, therefore skip everything and return catalog
                return self.catalogue, len(self.catalogue['eventID'])
            else:
                start_time = decimal_time(np.min(catalogue['year']), 
                                      1, 1, 0, 0, 0.)
        else:
            start_time = decimal_time(start_time.year, start_time.month,
                start_time.day, start_time.hour, start_time.minute,
                float(start_time.second))

        if not(end_time):
            end_time = datetime.now()
        end_time = decimal_time(end_time.year, end_time.month,
            end_time.day, end_time.hour, end_time.minute,
            float(end_time.second))

        is_valid = np.logical_and(self.time_value >= start_time,
                                  self.time_value < end_time)
        number_selected = np.sum(is_valid)
        is_valid = np.logical_not(is_valid)
        return purge_catalogue(catalogue, is_valid.astype(int)), \
            number_selected

    def select_within_depth_range(self, **kwargs):
        ''''''
        if (not 'upper_depth' in kwargs.keys()) and (not 'lower_depth' in 
            kwargs.keys()):
            '''No depths have been defined - so return entire catalogue!'''
            return catalogue, self.number_events

        upper_depth, lower_depth = _check_depth_limits(kwargs)
        is_valid = np.logical_and(self.catalogue['depth'] >= upper_depth,
                                  self.catalogue['depth'] < lower_depth)
        number_selected = np.sum(is_valid)
        is_value = np.logical_not(is_valid)
        return purge_catalogue(catalogue, is_valid.astype(int)), \
            number_selected

    def select_within_magnitude_range(self, lower_mag=None, upper_mag=None):
        ''''''
        if not lower_mag:
            lower_mag = -np.inf
        if not upper_mag:
            upper_mag = np.inf
        is_valid = np.logical_and(self.catalogue['magnitude'] >= lower_mag,
                                  self.catalogue['magnitude'] < upper_mag)

        number_selected = np.sum(is_valid)
        is_value = np.logical_not(is_valid)
        return purge_catalogue(catalogue, is_valid.astype(int)), \
            number_selected

