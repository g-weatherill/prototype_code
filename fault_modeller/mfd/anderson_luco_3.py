"""Module :mod: mfd.anderson_luco_3 implements class AndersonLucoType3"""
import numpy as np
from mfd.base import BaseMFDfromSlip

class AndersonLucoType3(BaseMFDfromSlip):
    '''Class to implement the 3nd fault activity rate calculator                
    of Anderson & Luco (1983)
    In retrospect the only difference between the Anderson & Luco calculators
    is in the very last formula - so this degree of abstraction is overkill!
    Still, you live and you learn :)'''         

    def setUp(self, mfd_conf):
        '''Parse the configuration parameters'''
        self.mfd_type = mfd_conf['Model_Type']
        self.mfd_weight = mfd_conf['Model_Weight']
        self.bin_width = mfd_conf['MFD_spacing'] 
        self.mfd_params = {'mmin': mfd_conf['Minimum_Magnitude'],
                           'mmax': None,
                           'bvalue': mfd_conf['b_value'][0]}
        self.occurrence_rate = []


    def get_mmax(self, mfd_conf, msr, rake, area):
        '''Gets the mmax for the fault - reading directly from the config file 
        or using the msr otherwise'''    
        if mfd_conf['Maximum_Magnitude']:
            self.mfd_params['mmax'] = mfd_conf['Maximum_Magnitude']
        else:
            self.mfd_params['mmax'] = msr.get_median_mag(rake, area)


    def get_mfd(self, slip, disp_length_ratio, shear_modulus, fault_width):                                           
        '''Calculates activity rate'''
        moment_scaling = [16.05, 1.5]
        beta = np.sqrt((disp_length_ratio * (10.0 ** moment_scaling[0])) /             
            ((shear_modulus * 1.0E10) * (fault_width * 1E5)))                 
        dbar = moment_scaling[1] * np.log(10.0)                               
        bbar = self.mfd_params['bvalue']* np.log(10.0)
        mag = np.arange(self.mfd_params['mmin'] - (self.bin_width / 2.), 
                        self.mfd_params['mmax'] + (1.5 * self.bin_width), 
                        self.bin_width)
        print 'Type 3 here'
        for ival in range(0, np.shape(mag)[0] - 2):
            self.occurrence_rate.append(
               self._cumulative_value(slip, self.mfd_params['mmax'],
               mag[ival], bbar, dbar, beta) - 
               self._cumulative_value(slip, self.mfd_params['mmax'],
               mag[ival + 1], bbar, dbar, beta))

    def _cumulative_value(self, slip, mmax, mag_value, bbar, dbar, beta):
        # Calculate N(M > mag_value) using Anderson & Luco Type 3 formula
        # Slip is input in mm/yr but needs to be converted to cm/yr - hence
        # divide by 10.
        return ((dbar * (dbar - bbar) / (bbar)) * ((slip / 10.) / beta) * 
                ((1./ bbar) * (np.exp(bbar * (mmax - mag_value)) - 1.) - 
                (mmax - mag_value)) * np.exp(-(dbar / 2.) * mmax))