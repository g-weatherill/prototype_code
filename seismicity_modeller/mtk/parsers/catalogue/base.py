#!/usr/bin/env/python

'''Abstract base class for catalogue parsers'''
import abc

class BaseCatalogueParser(object):
    '''Abstract base class for catalogue parser'''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def parse_catalogue(self, filename):
        '''Reads in the filename and parses the catalogue into a dictionary'''
        
