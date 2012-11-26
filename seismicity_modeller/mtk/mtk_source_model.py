#!/usr/bin/env/python

'''Basic Example of a source model class'''

import numpy as np
from math import fabs
from parsers.source_model.nrml_minimal_format import nrmlSourceModelParser
#from scientific.new_recurrence import Recurrence


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



    def get_gutenberg_richter_recurrence(self, catalogue, recurrence_config, 
        folder_path=None, file_type=None, file_dpi=None):
        '''For each source in the source model calculate recurrence parameters
        '''
        #recurrence = Recurrence()
        if folder_path:
            if not file_type:
                file_type = 'png'
            if not file_dpi:
                file_dpi = 300
            if not '/' in folder_path[-1]:
                folder_path = folder_path + '/'

        for source in self.source_model:
            # Find earthquakes from catalogue
            source.find_earthquakes_in_source(
                catalogue, 
                recurrence_config['buffer'][source.typology],
                recurrence_config['buffer']['upper_depth'],
                recurrence_config['buffer']['lower_depth'])
            
            if source.catalogue.number_events < 2:
                print 'Source %i, %s contains fewer than 2 events - skipping'\
                    % source.identifier, source.name
                    continue
            source.catalogue.recurrence(recurrence_config)
            source.mfd = source.catalogue.recurrence.mfd_parameters
            if folder_path:
                # Plot figure as output
                # Get filename
                figure_name = folder_path + '_' + str(source.identifier) +\
                    '_' + source.name + '.' + file_type
                mag_values, abs_rates, cum_rates = \
                    source.catalogue.get_completeness_adjusted_rates()
                source.catalogue.plot_completeness_adjusted_rates(
                    figure_name,                    
                    mag_values,
                    abs_rates,
                    cum_rates,
                    fig_dpi = file_dpi,
                    fig_format = file_type)

                                     
        

    def get_instrumental_mmax(self, catalogue, mmax_config):
        '''Calculate the instrumental maximum magnitude for each zone'''
        original_config = deepcopy(mmax_config)
        for source in self.source_model:
            # If a b-value and/or sigma-b is present in the source.mfd then
            # this temporarily overwrites the b-value in mmax_config
            if not source.mfd:
                source.mfd = {'Mmax': None, 'Mmax_Uncertainty': None}

            if 'bvalue' in source.mfd.keys():
                mmax_config['b-value'] = source.mfd['b-value']
            
            if 'sigma_bvalue' in source.mfd.keys():
                mmax_config['sigma b'] = source.mfd['sigma_bvalue']
            output, valid_config = source.catalogue.maximum_magnitude(
                mmax_config, 
                valid_config=False)
            source.mfd['Mmax'] = output['Mmax']
            source.mfd['Mmax_Uncertainty'] = output['Mmax_Uncertainty']


   

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
