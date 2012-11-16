#!/usr/bin/env/python

'''Abstract base class for source model parser'''

import abc

class BaseSourceParser(object):
    '''Abstract base class for catalogue parser'''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def parse_source_model(self, filename, config):
        '''Reads in the source model from the file '<filename>' and parses to a list of 
        nhlib sources classess'''
        
