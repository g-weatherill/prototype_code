#!/usr/bin/env/python
'''Example class for the uncertainty module'''
import abc
import numpy as np
import scipy.stats as sts


def _value_in_datalist(value, datalist):
    ''''''
    if value in datalist:
        return value
    else:
        return None

class Uncertainty(object):
    '''Trial uncertainty class'''
    #__metaclass__ = abc.ABCMeta
    
    def __init__(self, pref, **kwargs):
        '''Instantiate the class'''
        self.preferred = pref
        self.distribution = None
        self.datalist = kwargs.keys()
        self.upper_bound = _value_in_datalist('max', self.datalist)
        self.lower_bound = _value_in_datalist('min', self.datalist)
        self.quality = _value_in_datalist('quality', self.datalist)


    def get_samples(self, number_samples):
        '''Generate "number_samples" random samples from the 
        probability distribution'''
        return

class TruncatedGaussian(Uncertainty):
    '''Gaussian distributed random variate'''
    def __init__(self, pref, **kwargs):
        '''Setup the function'''
        super(TruncatedGaussian, self).__init__(pref, **kwargs)
        self.distribution = 'TruncatedGaussian'
        if kwargs['mean']:
            self.mean = kwargs['mean']
        if kwargs['standard_deviation']:
            self.standard_deviation = kwargs['standard_deviation']
        if kwargs['lower_bound']:
            self.lower_bound = kwargs['lower_bound']
        if kwargs['upper_bound']:
            self.upper_bound = kwargs['upper_bound']
    
    def get_samples(self, number_samples):
        '''Generate random sample from Gaussian distribution'''
        standard_samples = sts.truncnorm.rvs(
            (self.lower_bound - self.mean) / self.standard_deviation,
            (self.upper_bound - self.mean) / self.standard_deviation,
            loc = self.mean,
            scale=self.standard_deviation,
            size=number_samples
            )
        return standard_samples

class Gaussian(Uncertainty):
    '''Gaussian distributed random variate'''
    
    def __init__(self, pref, **kwargs):
        ''''''
        super(Gaussian, self).__init__(pref, **kwargs)
        self.distribution = 'Gaussian'
        if kwargs['mean']:
            self.mean = kwargs['mean']
        if kwargs['standard_deviation']:
            self.standard_deviation = kwargs['standard_deviation']
        
        
    def get_samples(self, number_samples):
        '''Generate random samples'''
        if not self.standard_deviation:
            return self.mean * np.ones(number_samples, dtype=float)
        else:
            return np.random.normal(self.mean, self.standard_deviation,
                                    number_samples)

class Uniform(Uncertainty):
    '''Uniformly distributed random variate'''
    def __init__(self, pref, **kwargs):
        ''''''
        super(Uniform, self).__init__(pref, **kwargs)
        self.distribution = 'Uniform'
        if kwargs['lower_bound']:
            self.lower_bound = kwargs['lower_bound']
        else:
            raise ValueError(
                'Lower bound must be defined for Uniform Distribution')
        if kwargs['upper_bound']:
            self.upper_bound = kwargs['upper_bound']
        else:
            raise ValueError(
                'Upper bound must be defined for Uniform Distribution')       
        
    def get_samples(self, number_samples):
        '''Generate random samples'''
        if abs(self.upper_bound - self.lower_bound) < 1E-12:
            return self.upper_bound * np.ones(number_samples, dtype=float)
        else:
            return np.random.uniform(self.lower_bound, 
                                     self.upper_bound, 
                                     number_samples)



