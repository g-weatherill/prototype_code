#!/usr/bin/env/python

'''Parser for pre-formatted matlab mat-file read from nrml 0.3'''

import numpy as np
from scipy.io import loadmat
from nhlib.geo.point import Point
from nhlib.geo.line import Line
from nhlib.geo.polygon import Polygon
from sources.mtk_point import mtkPointSource
from sources.mtk_area import mtkAreaSource
from sources.mtk_simple_fault import mtkSimpleFaultSource
from sources.mtk_complex_fault import mtkComplexFaultSource
from base import BaseSourceParser

class matlabNRML03parser(BaseSourceParser):
    '''Parser for nrml 0.3 parsed into a matlab .mat file
    '''
    def read_source_model(self, filename, config):
        '''Main reader'''
        input_source = loadmat(filename)
        input_source = input_source['source_model'][0]
        output_source = []
        for source in input_source:
            source_data = source[0]
            source_model, hypo = self.load_area_source(source_data, config)
            source_model.hypocentre_dist = hypo
            strike, dip, rake, rupture_rate_model = \
                self.get_rupture_rate_model(source_data)
            source_model.nodal_plane_dist = {'strike': strike,
                                            'dip': dip,
                                            'rake': rake}
            preloaded_rrm = False
            for keyval in ['aValue','bValue','minMag','maxMag']:
                if rupture_rate_model[keyval] != 0:
                    preloaded_rrm = True
                source_model.mfd = rupture_rate_model
            output_source.append(source_model)
                        


    def load_area_source(self, source_data, config):
        ''''''
        ident, name, region, rupture_depth, hypo_depth = \
            self.read_ancilliaries(source_data)
        return mtkAreaSource(ident, name, region, config['aspect_ratio'], 
            source_data['polygon']), hypo_depth


    def read_ancilliaries(self, source_dict):
        '''Return ID, name, tectonic region, rupture_depth, hypo_depth'''
        source_id = source_dict['name'][0].encode('utf-8')
        source_name = source_dict['id'][0].encode('utf-8')
        region = source_dict['region'][0].encode('utf-8')
        rupture_depth = source_dict['rupture_depth'][0]
        hypo_depth = source_dict['hypo_depth'][0]
        return source_id, source_name, region, rupture_depth, hypo_depth

   def get_rupture_rate_model(self, source_data)
       ''''''
       rrm = source_data[0][0]
       strike = rrm['strike'][0][0][0]
       dip = rrm['dip'][0][0][0]
       rake = rrm['rake'][0][0][0]
       model_type = rrm['model_type'][0][0].encode('utf-8')
       if 'truncatedGutenbergRichter' in model_type:
           mfd = {'model_type': model_type,
                  'avalue': rrm['aValue'][0][0][0], 
                  'bValue': rrm['bValue'][0][0][0], 
                  'minMag': rrm['minMag'][0][0][0], 
                  'maxMag': rrm['maxMag'][0][0][0]}
       else:
           mfd = None
       return strike, dip, rake, mfd
