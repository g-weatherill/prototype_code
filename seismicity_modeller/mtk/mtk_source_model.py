#!/usr/bin/env/python

'''Basic Example of a source model class'''

import numpy as np
from parsers.source_model.nrml_minimal_format import nrmlSourceModelParser
from scientific.new_recurrence import Recurrence


file_parser = {'nrml':  nrmlSourceModelParser()}
default_config = {'schema_validation': False,
                  'mesh_spacing': 1.0,
                  'bin_width': 0.1,
                  'area_discretisation': 1.0}

recurrence_config = {'buffer': {'Point': 20., 
                                'Area': None, 
                                'SimpleFault': 30., 
                                'ComplexFault': 30.,
                                'upper_depth': None,
                                'lower_depth': None},
                     'algorithm': 'Weichert',
                     'reference_magnitude': 4.0,
                     'magnitude_interval': 0.1,
                     'bvalue': 1.0,
                     'maxiter': 1000,
                     'itstab': 1E-5}


class mtkSourceModel(object):
    '''Source for a point class'''
    def __init__(self):
        '''Initialise'''
        self.source_model = None
        self.number_sources = None
        self.completeness = None
 
    def parse_source_model(self, filename, filetype, config = default_config):
        '''Read in the source model from nrml 0.4 format (using nrml tools)'''
        if filetype in file_parser.keys():
            self.source_model = \
                file_parser[filetype].read_source_model(filename, config)
        else:
            raise ValueError('File type not supported!')
        self.number_sources = len(self.source_model)



    def get_gutenberg_richter_recurrence(self, catalogue, recurrence_config):
        '''For each source in the source model calculate recurrence parameters
        '''
        recurrence = Recurrence()
        for source in self.source_model
            # Find earthquakes from catalogue
            source.find_earthquakes_in_source(
                catalogue, 
                recurrence_config['buffer'][source.typology],
                recurrence_config['buffer']['upper_depth'],
                recurrence_config['buffer']['lower_depth'])
            source.mfd = source.catalogue.recurrence(recurrence_config)

     
#    def get_gutenberg_richter_recurrence(self, catalogue, recurrence_config, 
#        selection_settings=None):
#        '''Loop throught the sources and calculate a- and b-values for each
#        source'''
#        recurrence = Recurrence()
#        for source in self.source_model:
#            if isinstance(selection_settings, dict):
#                if isinstance(source, mtkAreaSource):
#                    buffer_value = selection_settings['area_buffer']
#                elif isinstance(source, mtkPointSource):
#                    buffer_value = selection_settings['point_buffer']
#                else:
#                    buffer_value = selection_settings['fault_buffer']
#            else:
#                buffer_value = 0.0
#                selection_settings = {'upper_depth': 0.0, 
#                                      'lower_depth': 1.0E100}
#              
#            source.find_earthquakes_in_source(catalogue, buffer_value, 
#                selection_settings['upper_depth'], 
#                selection_settings['lower_depth'])
#            mfd_parameters = recurrence.run_recurrence(working_catalogue, 
#                                                       recurrence_config)    
