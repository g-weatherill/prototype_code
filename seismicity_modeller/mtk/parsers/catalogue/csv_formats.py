#!/usr/bin/env/python
'''Module csv parsers for catalogue'''

import csv
import numpy as np
from parsers.catalogue.base import BaseCatalogueParser

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


class CatalogueCSVParser(BaseCatalogueParser):
    '''Parses the catalogue from a csv format with named headers'''
    def parse_catalogue(self, catalogue, filename):
        '''Loads the csv file (filename) and parses to pre-formatted catalogue 
        dictionary'''
        filedata = open(filename, 'rt')
        data = csv.DictReader(filedata)
        valid_key_list = None
        for iloc, row in enumerate(data):
            if iloc == 0:
                # Run header checks
                valid_key_list = self._header_check(row.keys(), 
                    TOTAL_ATTRIBUTE_LIST)
            else:
                for key in valid_key_list:
                    if key in FLOAT_ATTRIBUTE_LIST:
                        catalogue[key] = self._float_check(catalogue[key], 
                                                           row[key])
                    elif key in INT_ATTRIBUTE_LIST:
                        catalogue[key] = self._int_check(catalogue[key],
                                                         row[key])
                    else:
                        catalogue[key].append(row[key])
        return catalogue


    def _header_check(self, input_keys, catalogue_keys):
        valid_key_list = []
        for element in input_keys:
            if element in catalogue_keys:
                valid_key_list.append(element)
            else:
                print 'Catalogue Attribute %s is not a recognised catalogue key'\
                    % element
        return valid_key_list

    def _float_check(self, attribute_array, value):
        '''Checks if value is valid float, appends to array is valid, appends
        nan if not'''
        value = value.strip(' ')
        if value:
            attribute_array = np.hstack([attribute_array, float(value)])
        else:    
            attribute_array = np.hstack([attribute_array, np.nan])
        return attribute_array

    def _int_check(self, attribute_array, value):
        '''Checks if value is valid integer, appends to array is valid, appends
        nan if not'''
        value = value.strip(' ')
        if value:
            attribute_array = np.hstack([attribute_array, int(value)])
        else:    
            attribute_array = np.hstack([attribute_array, np.nan])
        return attribute_array
    
