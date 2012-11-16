#!/usr/bin/env/python

'''Tectonic Regionalisation Class
'''

import numpy as np
from copy import deepcopy
from nhlib import scalerel


def _check_weight_sum(name, values, weights, tolerance = 1E-7):
    '''Checks that the number of values and weights are equal. Then
    check that the weights sum to 1, and return weights as numpy array
    if true and raise error if not'''
    if not isinstance(values, list):
        values = [values]
    if not isinstance(weights, list):
        weights = [weights]
    if len(values) != len(weights):
        raise ValueError('Weight array does not correspond to values in %s' 
                         % name)
    if np.abs(np.sum(np.array(weights)) - 1.) > tolerance:
        raise ValueError('Weights do not sum to 1 in %s!' %  name)
    return values, weights

class TectonicRegion(object):
    '''Basic class to define a single tectonic region'''
    def __init__(self, input_dict):
        '''Initialise parameters'''
        self.region_name = input_dict['Region_Type']
        self.identifier = input_dict['Region_Identifier_Code']
        # Parse MSR
        self.msr = input_dict['MagnitudeScalingRelation']
        # Ensure iterable
        if not isinstance(self.msr['Model'], list):
            self.msr['Model'] = [self.msr['Model']]
            self.msr['Weight'] = [self.msr['Weight']]
        
        # If Shear Modulus is included
        if 'Shear_Modulus' in input_dict.keys():
            self.shear_modulus = input_dict['Shear_Modulus']
            self.shear_modulus['Value'], self.shear_modulus['Weight'] = \
                self._check_weight_sum('Shear_Modulus', 
                                  self.shear_modulus['Value'], 
                                  self.shear_modulus['Weight'])
        else:
            self.shear_modulus = {'Value': [30.0], 'Weight': [1.0]}
        
        # If dlr is included
        if 'DisplacementLengthRatio' in input_dict.keys():
            self.dlr = input_dict['DisplacementLengthRatio']
            self.dlr['Value'], self.dlr['Weight'] = \
                self._check_weight_sum('DisplacementLengthRatio',
                                  self.dlr['Value'], 
                                  self.dlr['Weight'])
        else:
            self.dlr = None
        
        # If stress drop is included
        if 'StressDrop' in input_dict.keys():
            self.stress_drop = input_dict['StressDrop']
            self.stress_drop['Value'], self.stress_drop['Weight'] = \
                self._check_weight_sum('ShearModulus', 
                                  self.stress_drop['Value'], 
                                  self.stress_drop['Weight'])
        else:
            self.stress_drop = None

    def _check_weight_sum(self, key, values, weights, tolerance = 1E-7):
        '''Checks that the number of values and weights are equal. Then
        check that the weights sum to 1, and return weights as numpy array
        if true and raise error if not'''
        if not isinstance(values, list):
            values = [values]
        if not isinstance(weights, list):
            weights = [weights]
        if len(values) != len(weights):
            raise ValueError('Weight array does not correspond to values in %s' 
                             % key)
        if np.abs(np.sum(np.array(weights)) - 1.) > tolerance:
            raise ValueError('Weights do not sum to 1 in %s - %s!' 
                             % (self.region_name, key))
        return values, weights



class TectonicRegionalisation(object):
    '''Tectonic Regionalisation'''
    def __init__(self):
        self.regions = []
        self.number_regions = 0
        self.region_keys = []
        
    def parse_region_models_from_yaml(self, regionalisation):
        '''Parse a dictionary or list of dictionaries containing the 
        regionalisation
        '''
        if isinstance(regionalisation, dict):
            regionalisation = [regionalisation]
        if not isinstance(regionalisation, list):
            raise ValueError('Incorrect regionalisation formulation')
        
        for region_definition in regionalisation:
            region = TectonicRegion(region_definition)
            self.regions.append(region)
            self.region_keys.append(region.region_name)
            self.number_regions += 1
