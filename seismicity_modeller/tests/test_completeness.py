#!/usr/bin/env/python

'''First implementation of test class for completeness tools'''

import unittest
import numpy as np
from scientific.new_completeness import (check_completeness_config_master,
                                         Completeness, 
                                         Stepp1971, 
                                         AssumedMFD)

class CompletenessTestCase(unittest.TestCase):
    '''Main class for completeness tests'''
    def setUp(self):
        '''Instantiates the class with input data'''
        test_data = np.genfromtxt('data/completeness_test_cat.csv', 
                                  delimiter=',', skip_header=1)
        self.test_catalogue = {'year': test_data[:, 3], 
                          'magnitude': test_data[:, 17]}
        #test_config_stepp = {'algorithm': 'Stepp', 'delta_m': 0.1, 
        #                     'delta_t': 1.0, 'tolerance': 0.1, 
        #                     'increment_lock':True}
        self.supported_list = ['Stepp']
    
    def test_master_config_check(self):
        '''Tests the master config checker'''
        supported_list = ['Stepp']
        bad_config_missing_key = {'algorithm': 'Stepp', 'delta_t': 1.0,
                                  'increment_lock': True}
        with self.assertRaises(ValueError):
            check_completeness_config_master(bad_config_missing_key,
                                             self.supported_list)

        bad_config_bad_attribute = {'algorithm': 'Stepp', 'delta_m': -1.2, 
                                   'delta_t': 1.0, 'increment_lock': True}
        
        with self.assertRaises(AssertionError):
            check_completeness_config_master(bad_config_bad_attribute,
                                             self.supported_list)


    def test_Stepp_algorithm(self):
        '''Tests the core Stepp functionalities'''
        expected_table = np.genfromtxt('data/Stepp_Table_1.csv', 
                                       delimiter=',')
        good_Stepp_config =  {'algorithm': 'Stepp', 'delta_m': 0.1, 
                              'delta_t': 1.0, 'tolerance': 0.1, 
                              'increment_lock':True}
        bad_Stepp_config =  {'algorithm': 'Stepp', 'delta_m': 0.1, 
                              'delta_t': 1.0, 'increment_lock':True}
        Stepp = Stepp1971()
        with self.assertRaises(ValueError):
            Stepp.check_completeness_config(bad_Stepp_config,
                                            self.supported_list)
        completeness_table = Stepp.analyse(self.catalogue, good_Stepp_config)
        self.assertTrue(np.allclose(expected_table, completeness_table))


        
