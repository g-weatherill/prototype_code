#!/usr/bin/env/python
'''Prototype unittest code for mmax module'''

import unittest
import numpy as np
from scientific.mmax import (_get_observed_mmax, 
                             _get_magnitude_vector_properties,
                             KijkoSellevolFixedb, KijkoSellevolUncertainb,
                             KijkoNonParametricGaussian, CumulativeMoment)

class MmaxTestCase(unittest.TestCase):
    '''Testing class for Mmax functions'''
    def setUp(self):
        '''Set up the test class'''
        test_data = np.genfromtxt('data/completeness_test_cat.csv',
                                  delimiter=',',
                                  skip_header=1)
        self.catalogue = {'year': test_data[:, 3], 
                          'magnitude': test_data[:, 17],
                          'sigmaMagnitude': test_data[:, 18]}
        self.config = {'algorithm': None,
                       'input_mmax': ,
                       'input_mmax_uncertainty': ,
                       'maximum_iterations': ,
                       'tolerance': ,
                       'input_mmin': , 
                       'b-value': 1.0,
                       'sigma-b': 0.1,
                       'number_samples': 51,
                       'number_earthquakes': 100,
                       'number_bootstraps': }
        #self.completeness = np.array([])

    def test_get_observed_mmax(self):
        '''Tests the function to _get_observed_mmax from data'''

    def test_get_magnitude_vector_properties(self):
        '''Tests the function to retreive mmin and number of earthquakes if 
        required for certain functions'''
    
        
