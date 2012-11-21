#!/usr/bin/env/python
'''1st Draft of program to implement "SHIFT" for a 
general model'''
import abc
from copy import deepcopy
from matplotlib.nxutils import points_inside_poly
from linecache import getline, getlines

import numpy as np

#T5_column_header =        
# (/"CRB   ","CTF   ","CCB   ","OSRnor","OSRoth","OTFslo","OTFmed","OTFfas","OCB   ","SUB   "/)
#T5_CMT_pure_events =      
# (/   285.9,   198.5,   259.4,   424.3,    77.0,   398.0,   406.9,   376.6,   117.7,  2052.8/)
#T5_CMT_threshold_moment = 
# (/ 1.13E17,  3.5E17,  3.5E17, 1.13E17, 1.13E17,  2.0E17,  2.0E17,  2.0E17,  3.5E17,  3.5E17/)
#T5_beta =                 
# (/    0.65,    0.65,    0.62,    0.92,    0.82,    0.64,    0.65,    0.73,    0.53,    0.64/)
#T5_corner_magnitude =     
# (/    7.64,    8.01,    8.46,    5.86,    7.39,    8.14,    6.55,    6.63,    8.04,    9.58/)
#T5_tGR_model_moment_rate =
# (/ 1.67E12,  3.8E12, 1.06E13,  6.7E11,  1.9E11,  6.7E12,  9.4E11,  9.0E11,  4.6E12, 2.85E14/)
#T5_length_km =            
# (/  18126.,  19375.,  12516.,  61807.,  61807.,  27220.,  10351.,   6331.,  13236.,  38990./)
#T5_mean_velocity_mmpa =   
# (/   18.95,   21.54,   18.16,   46.40,    7.58,   20.68,   57.53,   97.11,   19.22,   61.48/)
#T5_assumed_dip_degrees =   
# (/     55.,     73.,     20.,     55.,     55.,     73.,     73.,     73.,     20.,     14./)
#T5_assumed_mu_GPa =       
# (/    27.7,    27.7,    27.7,    25.7,    25.7,    25.7,    25.7,    25.7,     49.,     49./)
#T5_lineIntegral_Nps =     
# (/   5.5E8,   4.4E8,   6.0E8,   5.0E9,   4.7E8,   5.2E8,   5.3E8,   5.5E8,   1.2E9, 1.58E10/)
#T5_coupledThickness_km =  
# (/     3.0,     8.6,     18.,    0.13,    0.40,     13.,     1.8,     1.6,     3.8,     18./)
#T5_assumed_lithosphere_km =
# (/      6.,     12.,     13.,      8.,      8.,     14.,     14.,     14.,     14.,     26./)
#T5_coupling =             
# (/    0.50,    0.72,    1.00,   0.016,    0.05,    0.93,    0.13,    0.11,    0.27,    0.69/)

CRB_PARAMS = {'CMT_EVENTS': 285.9, 'CMT_moment': 1.13E17,
              'beta': 0.65, 'corner_mag': 7.64, 'tGR_moment_rate': 1.67E12,
              'length': 18126., 'velocity': 18.95, 'assumed_dip': 55.,
              'assumed_mu': 27.7, 'line_integral': 5.5E8,
              'coupled_thickness': 3.0, 'lithosphere': 6.,
              'coupling': 0.50, 'adjustment_factor':1.001}
              
CTF_PARAMS = {'CMT_EVENTS': 198.5, 'CMT_moment': 3.5E17,
              'beta': 0.65, 'corner_mag': 8.01, 'tGR_moment_rate': 3.8E12,
              'length': 19375., 'velocity': 21.54, 'assumed_dip': 73.,
              'assumed_mu': 27.7, 'line_integral': 4.4E8,
              'coupled_thickness': 8.6, 'lithosphere': 12.,
              'coupling': 0.72, 'adjustment_factor':1.001}

CCB_PARAMS = {'CMT_EVENTS': 259.4, 'CMT_moment': 3.5E17,
              'beta': 0.62, 'corner_mag': 8.46, 'tGR_moment_rate': 1.06E13,
              'length': 12516., 'velocity': 18.16, 'assumed_dip': 20.,
              'assumed_mu': 27.7, 'line_integral': 6.0E8,
              'coupled_thickness': 18., 'lithosphere': 13.,
              'coupling': 1.0, 'adjustment_factor':1.001}

OSRnor_PARAMS = {'CMT_EVENTS': 424.3, 'CMT_moment': 1.13E17,
              'beta': 0.92, 'corner_mag': 5.86, 'tGR_moment_rate': 6.7E11,
              'length': 61807., 'velocity': 46.40, 'assumed_dip': 55.,
              'assumed_mu': 25.7, 'line_integral': 5.0E9,
              'coupled_thickness': 0.13, 'lithosphere': 8.,
              'coupling': 0.01, 'adjustment_factor': 1.619}

OSRoth_PARAMS = {'CMT_EVENTS': 77.0, 'CMT_moment': 1.13E17,
              'beta': 0.82, 'corner_mag': 7.39, 'tGR_moment_rate': 1.9E11,
              'length': 61807., 'velocity': 7.58, 'assumed_dip': 55.,
              'assumed_mu': 25.7, 'line_integral': 4.7E8,
              'coupled_thickness': 0.40, 'lithosphere': 8.,
              'coupling': 0.05, 'adjustment_factor': 1.619}

OTFslo_PARAMS = {'CMT_EVENTS': 398.0, 'CMT_moment': 2.0E17,
              'beta': 0.64, 'corner_mag': 8.14, 'tGR_moment_rate': 6.7E12,
              'length': 27220., 'velocity': 20.68, 'assumed_dip': 73.,
              'assumed_mu': 25.7, 'line_integral': 5.2E8,
              'coupled_thickness': 13., 'lithosphere': 14.,
              'coupling': 0.93, 'adjustment_factor': 1.619}

OTFmed_PARAMS = {'CMT_EVENTS': 406.9, 'CMT_moment': 2.0E17,
              'beta': 0.65, 'corner_mag': 6.55, 'tGR_moment_rate': 9.4E11,
              'length': 10351., 'velocity': 57.53, 'assumed_dip': 73.,
              'assumed_mu': 25.7, 'line_integral': 5.3E8,
              'coupled_thickness': 1.8, 'lithosphere': 14.,
              'coupling': 0.13, 'adjustment_factor': 1.619}

OTFfas_PARAMS = {'CMT_EVENTS': 376.6, 'CMT_moment': 2.0E17,
              'beta': 0.73, 'corner_mag': 6.63, 'tGR_moment_rate': 9.0E11,
              'length': 6331., 'velocity': 97.11, 'assumed_dip': 73.,
              'assumed_mu': 25.7, 'line_integral': 5.5E8,
              'coupled_thickness': 1.6, 'lithosphere': 14.,
              'coupling': 0.11, 'adjustment_factor': 1.619}

OCB_PARAMS = {'CMT_EVENTS': 117.7, 'CMT_moment': 3.5E17,
              'beta': 0.53, 'corner_mag': 8.04, 'tGR_moment_rate': 4.6E12,
              'length': 13236., 'velocity': 19.22, 'assumed_dip': 20.,
              'assumed_mu': 49., 'line_integral': 1.2E9,
              'coupled_thickness': 3.8, 'lithosphere': 14.,
              'coupling': 0.27, 'adjustment_factor': 2.000}

SUB_PARAMS = {'CMT_EVENTS': 2052.8, 'CMT_moment': 3.5E17,
              'beta': 0.64, 'corner_mag': 9.58, 'tGR_moment_rate': 2.85E14,
              'length': 38990., 'velocity': 61.48, 'assumed_dip': 14.,
              'assumed_mu': 49., 'line_integral': 1.58E10,
              'coupled_thickness': 18., 'lithosphere': 26.,
              'coupling': 0.69, 'adjustment_factor': 3.434}

IPL_PARAMS = {'CMT_EVENTS': 189.0, 'CMT_moment': 3.47E17,
              'beta': 0.63, 'corner_mag': 9.0, 'tGR_moment_rate': None,
              'length': None, 'velocity': None, 'area': 4.3536E14, 
              'assumed_dip': 14., 'assumed_mu': None, 'line_integral': None,
              'coupled_thickness': None, 'lithosphere': None, 'coupling': None, 
              'adjustment_factor': 1.619}

BIRD_GLOBAL_PARAMETERS = {'CRB': CRB_PARAMS,
                          'CTF': CTF_PARAMS,
                          'CCB': CCB_PARAMS,
                          'OSRnor': OSRnor_PARAMS,
                          'OSRoth': OSRoth_PARAMS,
                          'OTFslo': OTFslo_PARAMS,
                          'OTFmed': OTFmed_PARAMS,
                          'OTFas': OTFfas_PARAMS,
                          'OCB': OCB_PARAMS,
                          'SUB': SUB_PARAMS,
                          'IPL': IPL_PARAMS}

radian_conv = np.pi / 180.
secs_per_year = 365.25 * 24. * 60. * 60.
earth_radius = 6371000.0
# This value of 25.7474 is taken from Bird's analysis - TODO this needs to be
# generalised if integrating w/Modeller
CMT_duration_s = 25.7474 * secs_per_year
def moment_function(magnitude):
    '''Get moment from magnitude using Hanks & Kanamori'''
    return 10. ** ((1.5 * magnitude) + 9.05)

# Apply SI conversion adjustments from Bird (2007)'s code
# TODO This is ugly - reconsider this (maybe require only inputs in SI)
for reg_type in BIRD_GLOBAL_PARAMETERS.keys():
    reg = BIRD_GLOBAL_PARAMETERS[reg_type]
    reg['CMT_pure_event_rate'] = reg['CMT_EVENTS'] / CMT_duration_s
    reg['corner_moment'] = moment_function(reg['corner_mag'])
    if reg != 'IPL':
        reg['length'] = 1000. * reg['length']
        reg['velocity'] = (reg['velocity'] * 0.001) / secs_per_year
        reg['assumed_dip'] = reg['assumed_dip'] * radian_conv
        reg['assumed_mu'] = reg['assumed_mu'] * 1.0E9
        reg['coupled_thickness'] = reg['coupled_thickness'] * 1000.
        reg['lithosphere'] = reg['lithosphere'] * 1000.



def tapered_gutenberg_richter_cdf(moment, moment_threshold, beta, 
    corner_moment):
    '''Tapered Gutenberg Richter Cumulative Density Function'''
    cdf = np.exp((moment_threshold - moment) / corner_moment)
    return ((moment / moment_threshold) ** (-beta)) * cdf

def tapered_gutenberg_richter_pdf(moment, moment_threshold, beta,
    corner_moment):
    '''Tapered Gutenberg-Richter Probability Density Function'''
    return (beta / moment + 1. / corner_moment) * \
        tapered_gutenberg_richter_cdf(moment, moment_threshold, beta, 
        corner_moment)

# Regional data classes
#class RegionalParameters(object):
#    '''Abstract base class'''
#    __metaclass__ = abc.ABCMeta


class Strain(object):
    '''Primary class for SHIFT model'''
    def __init__(self):
        '''Initialise'''
        self.number_points = None
        self.strain_data = {}
        self.regions = None
        self.seismicity_rate = None
        #self.long = None
        #self.lat = None
        #self.exx = None
        #self.eyy = None
        #self.exy = None
        #self.var_exx = None
        #self.var_eyy = None
        #self.var_exy = None
        #self.minimum_magnitude = None
        #self.e1h = None
        #self.e2h = None
        self.second_invariant = None
        #self.dilatation_rate = None

    
    def get_strain_data(self, input_strain_file=None):
        '''Reads in the strain data'''
        if not input_strain_file:
            '''Assumes global data - read global file'''
            input_strain_file = 'global_strain_model.csv'
        
        header_names = getline(input_strain_file, 1).rstrip('\n')
        header_names = header_names.split(',')
        input_data = np.genfromtxt(input_strain_file, delimiter=',', 
                                   skip_header = 1)
        self.number_points = np.shape(input_data)[0]
        for iloc, name in enumerate(header_names):
            self.strain_data[name] = input_data[:, iloc]
        #print self.strain_data
        #breaker = here
        invalid_long = self.strain_data['long'] > 180.
        if np.any(invalid_long):
            # Strain Data is Pacific-centred
            self.strain_data['long'][invalid_long] = \
                self.strain_data['long'][invalid_long] - 360.
        if not 'region' in header_names:
            # Requires regionalisation
            print 'Deriving Kreemer Regionalisation'
            # Pre-allocate region list
            self.strain_data['region'] = \
                [ for xval in range(0, self.number_points)]

            west_limit = np.min(self.strain_data['long'])
            east_limit = np.max(self.strain_data['long'])
            north_limit = np.max(self.strain_data['lat'])
            south_limit = np.min(self.strain_data['lat'])
            regionalisation = self.apply_kreemer_regionalisation(north_limit,
                south_limit, east_limit, west_limit)
            for polygon in regionalisation:
                idlong = np.logical_and(
                    self.strain_data['long'] >= polygon['long_lims'][0],
                    self.strain_data['long'] < polygon['long_lims'][1])
                id0 = np.logical_and(idlong, 
                    self.strain_data['lat'] >= polygon['lat_lims'][0],
                    self.strain_data['lat'] < polygon['lat_lims'][1])
                if len(np.where(id0)[0]) > 0:
                    for iloc in id0:
                        self.strain_data['region'][iloc] = \
                            polygon['region_type']
        else:
            # If 'region' was in the header then the data dictionary will
            # contain a list of NaNs 
            datafile = open(input_strain_file, 'rt')
            for iloc, row in enumerate(csv.DictReader(datafile)):
                self.strain_data['region'][iloc] = row['region']
        self.strain_data['region'] = np.array(self.strain_data['region'])


    def apply_kreemer_regionalisation(self, north=90., south = -90., 
        east = 180., west = -180.):
        '''Applies the regionalisation of Kreemer (2003)'''
        input_data = getlines('data/tectonic_areas.dat.txt')
        kreemer_polygons = []
        temp_dict = {'cell': None, 'region_type':'', 'long_lims':None,
                     'lat_lims': None}
        
        for line_loc, line in enumerate(input_data):
            if '>' in line[0]:
                polygon_dict = deepcopy(temp_dict)
                polygon_dict['region_type'] = line[2].rstrip('\n')
                #print input_data, line_loc
                polygon_dict['cell'] = self._build_kreemer_cell(input_data, 
                                                                line_loc)
                polygon_dict['long_lims'] = np.array([
                    np.min(polygon_dict['cell'][:, 0]), 
                    np.max(polygon_dict['cell'][:, 0])])
                polygon_dict['lat_lims'] = np.array([
                    np.min(polygon_dict['cell'][:, 1]),
                    np.max(polygon_dict['cell'][:, 1])])
                polygon_dict['cell'] = None
                #print polygon_dict
                if polygon_dict['long_lims'][1] > 180:
                    polygon_dict['long_lims'] = \
                        polygon_dict['long_lims'] - 360.0
                if (polygon_dict['long_lims'][0] >= west) and\
                        (polygon_dict['long_lims'][1] <= east) and\
                        (polygon_dict['lat_lims'][0] >= south) and\
                        (polygon_dict['lat_lims'][1] <= north):
                    kreemer_polygons.append(polygon_dict)
        return kreemer_polygons        

        # Apply point in polygon
        
    def _build_kreemer_cell(self, data, loc):
        temp_poly = np.empty([4,2], dtype = float)
        for ival in range(1,5):
            value = data[loc + ival].rstrip('\n')
            value = value.lstrip(' ')
            value = np.array((value.split(' ',1))).astype(float)
            #print ival, value
            temp_poly[ival - 1,:] = value.flatten()
        return temp_poly


    def shift(self, minimum_magnitude): 
            #regionalisation = BIRD_GLOBAL_PARAMETERS):
        '''Applies SHIFT model - python implementation of the 
        SHIFT_continuum_seismicity subroutine in SHIFT_GSRM.f90'''
        threshold_moment = moment_function(minimum_magnitude)
        self.second_invariant = np.sqrt((self.strain_data['exx'] ** 2.) +
                                        (self.strain_data['eyy'] ** 2.) + 
                                        2.0 * (self.strain_data['exy'] ** 2.))
        dilatation_rate = self.strain_data['exx'] + self.strain_data['eyy']
        err_value = -1. * dilatation_rate
        center_normal_rate = (self.strain_data['exx'] + 
                              self.strain_data['eyy']) / 2.
        radius_rate = sqrt((self.strain_data['exx'] - center_normal_rate) ** 2.
                            + (self.strain_data['exy'] ** 2.))
        e1h = center_normal_rate - radius_rate
        e2h = center_normal_rate + radius_rate

        self.seismicity_rate = np.zeros(self.number_points, dtype=float)

        for region in self.regionalisation.keys():
            id0 = self.strain_data['region'] == region
            self.seismicity_rate[id0] = self.continuum_seismicity(
                threshold_moment, e1h[id0], e2h[id0], err_value[id0], 
                self.regionalisation[region])


    def continuum_seismicity(self, threshold_moment, e1h, e2h, err, 
        region_params):
        '''Function to implement the continuum seismicity calculation given
        vectors of input rates e1h, e2h [np.ndarray] and a dictionary of 
        the corresponding regionalisation params
        returns a vector of the corresponding seismicity rates'''
        
        strain_values = np.column_stack([e1h, e2h, err])
        e1_rate = np.min(strain_data, dimension = 1)
        e3_rate = np.max(strain_data, dimension = 1)
        e2_rate = 0. - e1_rate - e3_rate
        # Pre-allocate seismicity rate with zeros
        seismicity_rate = np.zeros(
            [np.shape(strain_values)[0], len(threshold_moment)], 
            dtype=float)
        # Calculate moment rate per unit area
        temp_e_rate = 2.0 * (-e1_rate)
        id0 = np.where(e2_rate < 0.0)[0]
        temp_e_rate[id0] = 2.0 * e3_rate[id0]
        M_persec_per_m2 = region_params['assumed_mu'] * temp_e_rate *\
                region_params['coupled_thickness']

        # Calculate seismicity rate at the threshold moment of the CMT
        # catalogue - Eq 6 in Bird et al (2010)
        seismicity_at_cmt_threshold = region_params['CMT_pure_event_rate'] * \
            (M_persec_per_m2 / region_params['tGR_moment_rate'])
        # Adjust forecast rate to desired rate using tapered G-R model
        # Taken from Eq 7 (Bird et al. 2010) and Eq 9 (Bird & Kagan, 2004)
        for iloc, moment_thresh in enumerate(threshold_moment):
            arg_check = (region_params['CMT_moment'] - moment_thresh) /\
		         region_params['corner_moment']
            id0 = arg_check >= -100.0
            seismicity_rate[id0, iloc] = ((moment_thresh / 
                region_params['CMT_moment']) ** region_pars['beta']) * \
                np.exp(arg_check[id0]) * seismicity_at_cmt_threshold * \
                region_params['adjustment_factor']
        return seismicity_rate
