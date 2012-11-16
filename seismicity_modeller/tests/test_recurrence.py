#!/usr/bin/env/python

'''First set of unit tests for recurrence calculators'''

import unittest
import numpy as np
from copy import deepcopy
from scientific.new_recurrence import (recurrence_table, b_max_likelihood
                                       _input_checks, Recurrence, 
                                       bMaximumLikelihood, KijkoSmit, 
                                       Weichert)

class RecurrenceTestCase(unittest.TestCase):
    '''Master unit-testing class for recurrence calculator'''
    def setUp(self):
        '''Set up general data'''
        test_data = np.genfromtxt('data/completeness_test_cat.csv',
                                  delimiter=',', 
                                  skip_header=1)
        self.catalogue = {'year': test_data[:, 3], 
                          'magnitude': test_data[:, 17]}
        self.config = {'algorithm': None,
                       'reference_magnitude': 4.0,
                       'magnitude_interval': 0.1,
                       'Average Type': 'Weighted'}
        self.completeness = np.array([[1990., 4.0], 
                                      [1960., 4.5], 
                                      [1910., 5.0]])

    def _cut_catalogue_to_completeness(self):
        '''Cuts the catalogue outside completeness range'''
        flag = np.ones(np.shape(self.catalogue['magnitude'])[0], dtype=bool)

        if not self.completeness:
            return self.catalogue

        for comp_val in self.completeness:
            id0 = np.logical_and(self.catalogue['year'] < comp_val[0],
                                 self.catalogue['magnitude'] < comp_val[1])
            flag[id0] = False
        cut_catalogue = {}
        for key in self.catalogue.keys():
            cut_catalogue[key] = self.catalogue[key][flag]
        
        return cut_catalogue

    def test_recurrence_table(self, catalogue, config):
        '''Basic recurrence table test'''
        target_data = np.genfromtxt('data/recurrence_table_test_1.csv',
                                    delimiter = ',')
        self.assert(np.allclose(target_data, recurrence_table(
            self.catalogue['magnitude'], self.config['magnitude_interval'],
            self.catalogue['year'])))

    def test_b_max_likelihood(self, catalogue, config, completeness):
        '''Test the basic maximum likelihood function'''
        # Cut the catalogue
        cut_catalogue = self._cut_catalogue_to_completeness()
        rec_table = recurrence_table(cut_catalogue['magnitude'],
                                     self.config['magnitude_interval'],
                                     cut_catalogue['year'])
        bvalue, sigmab = b_max_likelihood(rec_table[:, 0], rec_table[:, 1],
                                          self.config['magnitude_interval']
                                          np.min(cut_catalogue['magnitude']))
        self.assertAlmostEqual(bvalue, 1.02823300686, places=5)
        self.assertAlmostEqual(sigmab, 0.03262690837, places=5)

    def test_input_checks_completeness(self):
        '''Tests the input checks for different completeness inputs'''
        # No completeness input - should take the whole catalogue
        test_cmag, test_ctime, test_refmag, test_dmag = _input_checks(
            self.catalogue, self.config, completeness=None)
        self.assertTrue(np.allclose(test_cmag, np.array(4.0)))
        self.assertTrue(np.allclose(test_ctime, np.array(1910.)))
        self.assertAlmostEqual(test_refmag, 4.0)
        self.assertAlmostEqual(test_dmag, 0.1)
        
        # Single float value as input completeness
        test_cmag, test_ctime, test_refmag, test_dmag = _input_checks(
            self.catalogue, self.config, completeness=4.3)
        np.testing.assert_almost_equal(test_cmag, 4.3)
        np.testing.assert_almost_equal(test_ctime, 1910.)
        
        # Incorrectly configured table
        with self.assertRaises(ValueError) as ae:
            _input_checks(self.catalogue, self.config, 
                          completeness=np.array([4.0, 4.5]))
            _input_checks(self.catalogue, self.config, 
                          completeness=np.empty([3,3]))
        # Correctly configured table
        test_cmag, test_ctime, test_refmag, test_dmag = _input_checks(
            self.catalogue, self.config, self.completeness)
        self.assertTrue(np.allclose(test_cmag, np.array([4.0, 4.5, 5.0])))
        self.assertTrue(np.allclose(test_ctime, 
                                    np.array([1990., 1960., 1910.])))
        # Config file with non specified referenc magnitudes and intervals
        alt_config = deepcopy(self.config)
        alt_config['reference_magnitude'] = None
        alt_config['magnitude_interval'] = None
        test_cmag, test_ctime, test_refmag, test_dmag = _input_checks(
            self.catalogue, self.config, self.completeness)
        self.assertAlmostEqual(test_refmag, 0.0)
        self.assertAlmostEqual(test_dmag, 0.1)


    def test_bmaximum_likelihood_peripherals(self):
        '''Test suite for bMaximumLiklihood peripheral functions'''
        bml = bMaximumLikelihood()
        # Testing the means
        parameters = np.array([[1.3, 0.4, 3.5, 0.6],
                               [0.9, 0.1, 2.7, 0.4],
                               [1.0, 0.2, 3.0, 0.5]])
        neq = np.array([75, 25, 50])
        weighted_arithmetic = np.array([1.1333333, 0.2833333, 3.2, 0.5333333])
        self.assertTrue(np.allclose(weighted_arithmetic, 
                                    bml._weighted_mean(parameters, neq)))
        weighted_harmonic = np.array([1.10725552, 0.2181818, 
                                      3.1675977, 0.52173913])
        self.assertTrue(np.allclose(weighted_harmonic,
                                    bml._weighted_harmonic(parameters, neq)))
        # Testing caller function
        weighted_tuple = bml._average_parameters(parameters, neq, 
                                                 average_type='Weighted')
        self.assertAlmostEqual(weighted_tuple[0], weighted_arithmetic[0])
        self.assertAlmostEqual(weighted_tuple[1], weighted_arithmetic[1])
        self.assertAlmostEqual(weighted_tuple[2], weighted_arithmetic[2])
        self.assertAlmostEqual(weighted_tuple[3], weighted_arithmetic[3])
        
        harmonic_tuple = bml._average_parameters(parameters, neq, 
                                                 average_type='Harmonic')
        self.assertAlmostEqual(harmonic_tuple[0], weighted_harmonic[0])
        self.assertAlmostEqual(harmonic_tuple[1], weighted_harmonic[1])
        self.assertAlmostEqual(harmonic_tuple[2], weighted_harmonic[2])
        self.assertAlmostEqual(harmonic_tuple[3], weighted_harmonic[3])
        
        # Check for error if weighting does not equal parameters
        self.assertRaises(ValueError, bml._average_parameters(parameters,
            np.array([75, 25])))

    def test_bmaximum_likelihood(self):
        '''Tests the bMaximumLikelihood function for different requirements'''
        bml = bMaximumLikelihood()
        alternate_config = deepcopy(self.config)
        alternate_config['reference_magnitude'] = None
        grvars = bml.calculate(self.catalogue,
                               self.config,
                               self.completeness)

        self.assertAlmostEqual(grvars[0], 0.9324271)
        self.assertAlmostEqual(grvars[1], 0.02933169)
        self.assertAlmostEqual(grvars[2], 54.34166455)
        self.assertAlmostEqual(grvars[3], 2.240960457)

        grvars = bml.calculate(self.catalogue,
                               alternate_config,
                               self.completeness)
        self.assertAlmostEqual(grvars[2], 5.6519640)
        self.assertAlmostEqual(grvars[3], 0.0199687)

    def test_kijko_smit(self):
        '''Tests for the Kijko & Smit b-value estimator'''
        # Harmonic mean function
        weighted_harmonic = np.array([1.10725552, 0.2181818, 
                                      3.1675977, 0.52173913])
        bks = KijkoSmit()
        self.assertAlmostEqual(1.02768166, bks._harmonic_mean(
            np.array([1.1, 0.9, 1.0]), np.array([75, 25, 50])))
        # a-value function
        self.assertAlmostEqual(27.358423434, 
            bks._calculate_a_value(0.92, 2501.0, np.array([20., 50., 100.,]),
            np.array([4.0, 4.5, 5.0]), 4.0))
        # Full function - reference magnitude 4.0
        grvars = bks.calculate(self.catalogue,
                               self.config,
                               self.completeness)

        self.assertAlmostEqual(grvars[0], 0.923418065)
        self.assertAlmostEqual(grvars[1], 0.0184646687)
        self.assertAlmostEqual(grvars[2], 27.45136839)
        self.assertAlmostEqual(grvars[3], 0.308825369)
         
        alternate_config = deepcopy(self.config)
        alternate_config['reference_magnitude'] = None

        grvars = bks.calculate(self.catalogue,
                               alternate_config,
                               self.completeness)

        self.assertAlmostEqual(grvars[2], 3.042135567)
        self.assertAlmostEqual(grvars[3], 0.036941254)



