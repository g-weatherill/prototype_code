#!/usr/bin/env/python
'''Example class for the uncertainty module'''
import abc
import numpy as np
import scipy.stats as sts

class Uncertainty(object):
    '''Trial uncertainty class'''
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, **kwargs):
        '''Instantiate the class'''
        self.preferred = None
        self.distribution = None
        self.quality = None
        #self.datalist = kwargs.keys()
    
    @abc.abstractmethod
    def create(self, preferred, **kwargs):
        '''Creates the distribution'''
        self.preferred = preferred
        keylist = kwargs.keys()
        if 'quality' in keylist:
            self.quality = kwargs['quality']
        if 'max' in keylist:
            self.upper_bound = kwargs['max']
        else:
            self.upper_bound = None
        if 'min' in keylist:
            self.lower_bound = kwargs['min']
        else:
            self.lower_bound = None


    
    @abc.abstractmethod
    def get_samples(self, number_samples):
        '''Generate "number_samples" random samples from the 
        probability distribution'''
        

class TruncatedGaussian(Uncertainty):
    '''Gaussian distributed random variate'''
    
    def create(self, preferred, **kwargs):
        '''Set up the distribution'''
        self.distribution = 'TruncatedGaussian'
        self.mean = kwargs['mean']
        if not preferred:
            self.preferred = self.mean
        else:
            self.preferred = preferred
            
        self.standard_deviation = kwargs['standard_deviation']
        if not kwargs['min']:
            self.lower_bound = -np.inf
        else:
            self.lower_bound = kwargs['min']
        if not kwargs['max']:
            self.upper_bound = np.inf
        else:
            self.upper_bound = kwargs['max']
        if self.upper_bound < self.lower_bound:
            raise ValueError('Upper bound of truncated Gaussian must be '
                             'greater than or equal to lower bound!')
        if 'quality' in kwargs.keys():
            self.quality = kwargs['quality']
        

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
    
    def create(self, preferred, **kwargs):
        ''''''
        
        self.distribution = 'Gaussian'
        self.mean = kwargs['mean']
        if not preferred:
            self.preferred = self.mean
        else:
            self.preferred = preferred
        self.standard_deviation = kwargs['standard_deviation']
        if 'quality' in kwargs.keys():
            self.quality = kwargs['quality']
        

    def get_samples(self, number_samples):
        '''Generate random samples'''
        if not self.standard_deviation:
            return self.mean * np.ones(number_samples, dtype=float)
        else:
            return np.random.normal(self.mean, self.standard_deviation,
                                    number_samples)

class Uniform(Uncertainty):
    '''Uniformly distributed random variate'''
        
    def create(self, preferred,  **kwargs):
        ''''''
        #__super__(Uniform, self).__init__()
        self.distribution = 'Uniform'
        self.lower_bound = kwargs['min']
        self.upper_bound = kwargs['max']
        if not preferred:
            # If preferred is not specified - take midpoint of distribution
            self.preferred = (self.lower_bound + self.upper_bound) / 2.
        else:
            self.preferred = preferred
        if 'quality' in kwargs.keys():
            self.quality = kwargs['quality']
        
        
    def get_samples(self, number_samples):
        '''Generate random samples'''
        if abs(self.upper_bound - self.lower_bound) < 1E-12:
            return self.upper_bound * np.ones(number_samples, dtype=float)
        else:
            return np.random.uniform(self.lower_bound, 
                                     self.upper_bound, 
                                     number_samples)
