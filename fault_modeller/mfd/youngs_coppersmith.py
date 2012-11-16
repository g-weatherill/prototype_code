"""Module :mod: mfd.youngs_coppersmith implements class 
YoungsCoppersmith"""

import numpy as np
from math import exp, log, log10
from mfd.base import BaseMFDfromSlip

class YoungsCoppersmith(BaseMFDfromSlip):
    '''calculates model from Youngs and Coppersmith (1985)'''
    def setUp(self, mfd_conf):
        '''Initialise parameters'''
        self.mfd_type = mfd_conf['Model_Type']
        self.mfd_weight = mfd_conf['Model_Weight']
        self.bin_width = mfd_conf['MFD_spacing']
        self.mfd_params = {'mmax': None,
                           'mmin': mfd_conf['Minimum_Magnitude'],
                           'delta_m1': mfd_conf['delta_m1'],
                           'delta_m2': mfd_conf['delta_m2'],
                           'bvalue': mfd_conf['b_value']}
        self.occurrence_rate = None
        self.m_c = None

    def get_mmax(self, mfd_conf, msr, rake, area):
        '''Gets the mmax for the fault - reading directly from the config file 
        or using the msr otherwise'''    
        if mfd_conf['Maximum_Magnitude']:
            self.mfd_params['mmax'] = mfd_conf['Maximum_Magnitude']
        else:
            self.mfd_params['mmax'] = msr.get_median_mag(rake, area)

        self._delta_defaults()

    def _delta_defaults(self):
        '''If the delta_1 and delta_2 values are not defined in the input
        model then use the defaults defined by Youngs & Coppersmith (1985):
        delta_m2 = 0.5
        delta_m1 = 1.0'''
        if not self.mfd_params['delta_m1']:
            self.mfd_params['delta_m1'] = 1.0

        if not self.mfd_params['delta_m2']:
            self.mfd_params['delta_m2'] = 0.5
        
        self.m_c = self.mfd_params['mmax'] - self.mfd_params['delta_m2']
        if self.m_c < self.mfd_params['mmin']:
            raise ValueError('Characteristic magnitude lower than minimum'\
                'magnitude')

    def get_mfd(self, slip, shear_modulus, area):
        '''Calculate magnitude frequency distribution using the Youngs & 
        Coppersmith (1985) model'''
        print self.mfd_params
        beta = self.mfd_params['bvalue'] * log(10.0)
        magnitudes = np.arange(self.mfd_params['mmin'], 
                               self.mfd_params['mmax'] + self.bin_width, 
                               self.bin_width)
        self.occurrence_rate = np.zeros(len(magnitudes), dtype = float)
        c_value = self._c_factor(beta)
        print magnitudes, self.mfd_params['mmin'], self.m_c, c_value
        iloc = np.where(np.logical_and(magnitudes >= self.mfd_params['mmin'], 
                                    magnitudes < self.m_c))[0]
        self.occurrence_rate[iloc] = ((beta * np.exp(-beta * (
            magnitudes[iloc] - self.mfd_params['mmin']))) /\
            (1. - np.exp(-beta * (self.mfd_params['mmax'] - 
            self.mfd_params['mmin'] - self.mfd_params['delta_m2'])))) * c_value

        iloc = np.where(np.logical_and(magnitudes >= self.mfd_params['mmin'], 
                                    magnitudes >= self.m_c))[0]
        self.occurrence_rate[iloc] = ((beta * np.exp(-beta * (
            self.mfd_params['mmax'] - self.mfd_params['mmin'] - 
            self.mfd_params['delta_m1'] - self.mfd_params['delta_m2']))) /
            (1. - np.exp(-beta * (self.mfd_params['mmax'] -
            self.mfd_params['mmin'] - self.mfd_params['delta_m2'])))) * c_value

        
    def _c_factor(self, beta):
        '''Returns the c factor'''
        c_value = (beta * exp(-beta * (self.mfd_params['mmax'] - 
            self.mfd_params['mmin'] - self.mfd_params['delta_m1'] - 
            self.mfd_params['delta_m2']))) / (1.0 - 
            exp(-beta * (self.mfd_params['mmax'] - self.mfd_params['mmin'] -
            self.mfd_params['delta_m2'])))
        return 1.0 / (1.0 + (c_value * self.mfd_params['delta_m2']))