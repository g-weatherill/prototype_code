#!/usr/bin/env/python

'''Master script for testing out declustering utilities
'''
import numpy as np
from mtk_catalogue import MTKCatalogue
from parsers.csv_formats import CatalogueCSVParser
from scientific.declustering import Declustering


ifile = 'parsers/emme_test_catalogue.csv'
# Instantiate catalogue
Catalogue = MTKCatalogue()
parser1 = CatalogueCSVParser()

# Parse catalogue
Catalogue.catalogue = parser1.parse_catalogue(Catalogue.catalogue, ifile)
#print Catalogue.catalogue
# Build config
config1 = {'algorithm':'GardnerKnopoffType1', 
           'window_opt': 'GardnerKnopoff', 
           'fs_time_prop': 1.0,
           'purge': True}

# Instantiate Declustering
decluster1 = Declustering()
print config1, decluster1.decluster_master
Catalogue.catalogue, vmainshock = decluster1.run_declustering(
    Catalogue.catalogue, config1)

print Catalogue.catalogue['cluster'], Catalogue.catalogue['flag_vector']
print vmainshock


