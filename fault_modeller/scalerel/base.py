"""
Module :mod:`scalerel.base` defines an abstract base classes
for :class:`ASR <BaseASR>`
"""
import abc

class BaseASR(object):
    '''Quick copy of base ASR class - TODO: Overwrite with nhlib implementation
    '''
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def get_std_dev_mag(self, rake):
        '''Return Std Dev'''
    @abc.abstractmethod
    def get_median_mag(self, area, rake):
        '''Return median area from base ASR'''
   
    def get_mag(self, area, rake, epsilon = 0.0):
        '''Returns the expected are from magnitude'''
        median_mag = get_median_mag(area, rake)
        std_dev = get_std_dev_mag(rake)
        return median_mag + epsilon * std_dev