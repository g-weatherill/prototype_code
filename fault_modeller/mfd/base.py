"""
Module :mod:`mfd.base` defines an abstract base classes
for :class:`BaseMFDfromSlip>`
"""
import abc

class BaseMFDfromSlip(object):
    '''Base class for calculating magnitude frequency distribution
    from a given slip value'''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def setUp(self, mfd_conf):
        '''Initialises the parameters from the mfd type'''

    @abc.abstractmethod
    def get_mmax(self, mfd_conf, msr, rake, area):
        '''Gets the mmax for the fault - reading directly from the config file 
        or using the msr otherwise'''

    @abc.abstractmethod
    def get_mfd(self):
        '''Calculates the magnitude frequency distribution'''