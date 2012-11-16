"""Module :mod: mfd.characteristic implements class Characteristic"""

import numpy as np
from scipy.stats import truncnorm
from math import exp, log, log10
from mfd.base import BaseMFDfromSlip

class Characteristic(BaseMFDfromSlip):
    '''Calculates Characteristic earthquake model assuming a truncated 
    Gaussian distribution and the method of Wesnousky (1986)'''
    def setUp(self, mfd_conf):
        '''Parse the configuration parameters'''
        self.mfd_type = mfd_conf['Model_Type']
        self.mfd_weight = mfd_conf['Model_Weight']
        self.bin_width = mfd_conf['MFD_spacing'] 
        self.mfd_params = {'Sigma': mfd_conf['Sigma'],
                           'mmax': None,
                           'Lower': mfd_conf['Lower_Bound'],
                           'Upper': mfd_conf['Upper_Bound'],
                           'mmin': None}
        self.occurrence_rate = []

    def get_mmax(self, mfd_conf, msr, rake, area):
        '''Gets the mmax for the fault - reading directly from the config file 
        or using the msr otherwise'''    
        if mfd_conf['Maximum_Magnitude']:
            self.mfd_params['mmax'] = mfd_conf['Maximum_Magnitude']
        else:
            self.mfd_params['mmax'] = msr.get_median_mag(rake, area)

    def get_mfd(self, slip, shear_modulus, area):                                           
        '''Calculates activity rate'''
        self.mfd_params['mmin'] = self.mfd_params['mmax'] + (
            self.mfd_params['Lower'] * self.mfd_params['Sigma'])
        moment_rate = (shear_modulus * 1.E9) * (area * 1.E6) * (slip / 1000.)
        mag_upper = self.mfd_params['mmax'] + (self.mfd_params['Upper'] *
                                               self.mfd_params['Sigma'])
        mag_range = np.arange(self.mfd_params['mmin'], 
                              mag_upper + self.bin_width + 1E-7,
                              self.bin_width)
        #moment_mag = 10. ** (1.5 * mag_range + 9.05)
        moment_mag = 10. ** (1.5 * self.mfd_params['mmax']  + 9.05)
        characteristic_rate = moment_rate / moment_mag

        self.occurrence_rate = characteristic_rate * (truncnorm.cdf(
            mag_range + (self.bin_width / 2.), self.mfd_params['Lower'], 
            self.mfd_params['Upper'], loc=self.mfd_params['mmax'], 
            scale=self.mfd_params['Sigma']) - truncnorm.cdf(
            mag_range - (self.bin_width / 2.), self.mfd_params['Lower'], 
            self.mfd_params['Upper'], loc=self.mfd_params['mmax'], 
            scale=self.mfd_params['Sigma']))
        #number_mags = np.shape(mag_range)[0]
        #prob_mass_func = truncnorm.pdf(np.zeros(number_mags), 
        #                               self.mfd_params['Lower'],
        #                               self.mfd_params['Upper'], loc = 0., 
        #                               scale = 1.)
        
        #prob_mass_func = prob_mass_func / np.sum(prob_mass_func)
        #moment_rate = moment_rate * prob_mass_func
        #self.occurrence_rate = moment_rate / moment
        #self.occurrence_rate[:-1] = self.occurrence_rate[:-1] - \
        #                            self.occurrence_rate[1:]
