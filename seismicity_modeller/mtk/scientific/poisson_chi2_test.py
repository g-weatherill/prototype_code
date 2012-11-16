#!/usr/bin/env python

'''Prototype script for implementation of a Poisson chi2 test suite.
Implementation based on descriptions found in Luen & Stark (2012)'''

import numpy as np
from scipy.misc import factorial
from scipy.stats import chi2

class PoissonChi2Tests(object):
    '''Class for implementation of Poisson Chi^2 tests'''
    def __init__(self, catalogue, config):
        '''Initialise'''
        self.year = catalogue['year']
        self.month = catalogue['month']
        self.day = catalogue['day']
        self.magnitude = catalogue['magnitude']
        self.longitude = catalogue['longitude']
        self.latitude = catalogue['latitude']
        self.tester_map = {}
        self.start_year = np.min(year)
        self.end_year = np.max(year)
        self.duration = float(self.end_year - self.start_year + 1)
        self.number_events = len(year)

    def multinomial_chi2_test(self, config):
        '''Simple multinomial chi2 test - Based on Gardner & Knopoff (1974)'''
        # Config should include 'K' marker - divide catalogue into intervals
        # of K length
        if config['K'] < 1:
            raise ValueError('K must be greater than or equal to 1')
        start_bin = range(self.start_year, self.end_year, config['K'])
        end_bin = range(self.start_year + config['K'] - 1, self.end_year,
                        config['K'])
        number_ints = len(end_bin)
        time_ints = np.column_stack([np.array(start_bin[:number_ints]),
                                     np.array(end_bin)])
        
        ncount = np.zeros(np.shape(time_ints)[0], dtype=int)

        for iloc, time_bin in enumerate(time_ints):
            ncount[iloc] = np.sum(np.logical_and(self.year >= time_bin[0],
                self.year < time_bin[1]))
        #ncount = ncount.astype(float)
        theoretical_rate = self.number_events / float(config['K'])
        c_value, expected_c = self._get_c_value(config['K'], theoretical_rate)
        observed_c = self._get_obs_c(c_value, ncount)
        chi2m = np.sum((observed_c.astype(float) - 
            (expected_c.astype(float) ** 2.)) / expected_c.astype(float))
        if not(config['dof']):
            config['dof'] = float(c_value[-1] - 2)

        p_value = chi2.sf(chi2m, config['dof'])
        return p_value, chi2m, c_value[-1], config['dof']


    def _get_c_value(self, k_value, rate):
        '''Get C-value satisfying criterion that expected rate for all 
        intervals >= 5'''
        c_value = 2
        breaker = True
        while breaker 
            e_c = self._get_expected_c(c_value, k_value, rate)
            if np.all(e_c >= 5.):
                # Found minimum C for which all Ec >= 5
                breaker = False
            c_value += 1
            if c_value > self.number_events:
                raise ValueError('Expected C reached maximum number of '
                    'iterations. Change parameter of function')

       return np.range(0, c_value), e_c



    def _get_expected_c(self, c_value, k_value, rate):
        '''Calculates the expected c-value (Ec)'''
        k_value = float(k_value)
        trial_c = np.arange(0, cvalue, dtype=int)
        exp_c = k_value * np.exp(rate) * \
            (rate ** trial_c.astype(float) / factorial(trial_c))
        exp_c[-1] = k_value - np.sum(exp_c[:-1])
        return exp_c
    
    def _get_obs_c(self, c_values, ncount):
        '''Returns the Oc values - the number of intervals for which the
        number of events in the interval is equal to C'''
        obs_c = np.zeros(len(c_values), dtype=int)
        for c_value in enumerate(c_values):
            if c_value == np.max(c_value):
                obs_c[c_value] = np.sum(ncount >= c_value)
            else:
                obs_c[c_value] = np.sum(ncount == c_value)
        return obs_c



        
