#!/usr/bin/env/python

'''Demo prototype of a fault calculation class'''
import abc
import numpy as np
import yaml
from math import log, log10, sin, exp
from geo_utils import haversine, get_trace_length
from fault_geometry import Simple
from mfd.characteristic import Characteristic
from mfd.youngs_coppersmith import YoungsCoppersmith
from mfd.anderson_luco_1 import AndersonLucoType1
from mfd.anderson_luco_2 import AndersonLucoType2
from mfd.anderson_luco_3 import AndersonLucoType3
from scalerel.wc1994 import WC1994

class FaultModel:
    '''Fault class'''
    def __init__(self):
        '''Initialise'''
        self.tectonic_model = None
        self.fault_model = None
        self.number_faults = 0
        self.mfd_spacing = 0.1
        self.tectonic_regionalisation = None
        self.faults = None

    def read_fault_input_model(self, input_filename):
        '''Read the configuration file'''
        input_data = yaml.load(open(input_filename, 'rt'))
        self.tectonic_model = input_data['tectonic_regionalisation']
        self.fault_model = input_data['Fault_Model']
        self.number_faults = len(self.fault_model)
        self.tectonic_regionalisation = TectonicRegion()
        self.tectonic_regionalisation.parse_region_models(self.tectonic_model)

    def get_activity_rate(self, output_filename):
        '''Calculate the activity rates on each fault'''
        file_id = open(output_filename, 'wt')
        if isinstance(self.fault_model, list):
            for fault_settings in self.fault_model:
                fault = Fault(fault_settings)
                fault.create_fault_definition(fault_settings,
                                              self.tectonic_regionalisation)
                fault.get_recurrence_model(fault_settings)
                self._write_to_output_file(fault, file_id)
        else:
            fault = Fault(self.fault_model)
            fault.create_fault_definition(self.fault_model,
                                          self.tectonic_regionalisation)
            fault.get_recurrence_model(self.fault_model)
            self._write_to_output_file(fault, file_id)
        file_id.close()

    def _write_to_output_file(self, fault, file_id):
        '''Writes output to simple ascii'''
        print >> file_id, 'Fault Name: %s' % fault.fault_name
        print >> file_id, 'Fault ID: %10d' % fault.fault_id
        print >> file_id, 'Fault Region: %s' % fault.fault_region
        print >> file_id, '%s Fault Typology' % fault.typology
        print >> file_id, 'Trace:'
        for row in fault.fault_geometry.trace:
            print >> file_id, '%12.5f, %12.5f' %(row[0], row[1])
        print >> file_id, 'Upper Seismogenic Depth: %12.5f' \
            % fault.fault_geometry.upper_depth
        print >> file_id, 'Lower Seismogenic Depth: %12.5f' \
            % fault.fault_geometry.lower_depth
        print >> file_id, 'Dip: %8.1f' % fault.fault_geometry.dip
        print >> file_id, 'Rake: %8.1f' % fault.rake
        for rec_mod in fault.recurrence_model:
            print >> file_id, 'Occurrence Rate Mmin: %8.2f, binSize: %8.3f' \
                % (rec_mod.mfd_params['mmin'], rec_mod.bin_width)
            occ_rate_string = ''
            occ_rate_string = occ_rate_string.join(
                [str(value).format('%15.8f') + ', ' 
                for value in rec_mod.occurrence_rate])
            print >> file_id, '%s' % occ_rate_string[:-2]
            

class TectonicRegion:
    '''Creates instances of tectonic region class'''
    def __init__(self):

        self.tect_reg = []
        self.number_regions = None
        self.key_list = []


    def parse_region_models(self, tect_reg_list):
        
        if isinstance(tect_reg_list, list):
            # Many region types
            self.number_regions = len(tect_reg_list)
            self.tect_reg = dict()
            for region in tect_reg_list:
                self.tect_reg[region['Region_Type']] = \
                    self._rearrange_tect_region_dictionary(region)
                #self.tect_reg.append(tr_dict)
                self.key_list.append(region['Region_Type'])
        elif isinstance(tect_reg_list, dict):
            # Only one region type
            self.tect_reg = dict()
            self.tect_reg[tect_reg_list['Region_Type']] = \
            self._rearrange_tect_region_dictionary(tect_reg_list)
            self.tect_reg = [self.tect_reg]
            self.key_list.append(region['Region_Type'])
        else:
            #No regionalisation or incorrect formulation
            print 'No tectonic region model - or incorrect formulation'''
            print 'If any fault model is missing these attributes an error'
            print 'will be raised shortly - you were warned bitches!'

    
    def _check_weight_sum(self, values, weights, tolerance = 1E-7):
        '''Checks that the number of values and weights are equal. Then
        check that the weights sum to 1, and return weights as numpy array
        if true and raise error if not'''
        if isinstance(values, list) and isinstance(weights, list) and\
            len(values) != len(weights):
            raise ValueError('Weight array does not correspond to values')
        if np.abs(np.sum(np.array(weights)) - 1) > tolerance:
            raise ValueError('Weights do not sum to 1!')
        return weights

    def _rearrange_tect_region_dictionary(self, tr_dict):
        '''Re-arrange the dictionary structure of the tectonic region classes'''
        new_dict = dict()
        new_dict['Strain_Drop_Values'] = \
            tr_dict['DisplacementLengthRatio']['Value']
        new_dict['Strain_Drop_Weight'] = self._check_weight_sum(
            new_dict['Strain_Drop_Values'],
            tr_dict['DisplacementLengthRatio']['Weight'])
        
        new_dict['MSR_Models'] = \
            tr_dict['MagnitudeScalingRelation']['Model']
        new_dict['MSR_Models_Weight'] = self._check_weight_sum(
            new_dict['MSR_Models'],
            tr_dict['MagnitudeScalingRelation']['Weight'])

        new_dict['Shear_Modulus'] = \
            tr_dict['Shear_Modulus']['Value']
        new_dict['Shear_Modulus_Weight'] = self._check_weight_sum(
            new_dict['Shear_Modulus'],
            tr_dict['Shear_Modulus']['Weight'])
        new_dict['ID'] = tr_dict['Region_Identifier_Code']
        return new_dict


class Fault(object):
    '''Creates the fault definition'''
    def __init__(self, fault_dict):
        '''Parses the fault dictionary representation into the corresponding
        fault class attributes'''
        self.fault_id = fault_dict['ID']
        self.fault_name = fault_dict['Fault_Name']
        self.fault_region = fault_dict['tectonic_region']
        self.typology = fault_dict['Fault_Geometry']['Fault_Typology']
        self.fault_geometry = Simple()
        self.slip_type = fault_dict['slip_type']
        self.slip_completeness = fault_dict['slip_completeness_factor']
        self.rake = fault_dict['rake']
        self.MFD = []
        self.slip = None
        self.slip_min = None
        self.slip_max = None
        self.slip_quality = None
        self.aseismic = None
        self.aseismic_completeness = None
        self.megazone = fault_dict['Megazone']
        self.shear_modulus = fault_dict['Shear_Modulus']
        self.msr = fault_dict['MagnitudeScalingRelation']
        self.aspect_ratio = fault_dict['AspectRatio']
        self.disp_length_ratio = fault_dict['DisplacementLengthRatio']
        self.recurrence_model = []

    def create_fault_definition(self, fault_dict, tect_regions):
        '''Filles out the general fault definitions'''
        # Slip model
        self._parse_fault_slip(fault_dict['slip']) 
        # Aseismic Data
        self._parse_aseismic_slip(fault_dict['aseismic'])
        # Define the fault geometries
        self.fault_geometry.get_fault_extent(fault_dict['Fault_Geometry'])
        self.fault_geometry.get_basic_parameters(fault_dict['Fault_Geometry'])
        # Fill in missing attributes - check if fault has a corresponing
        if not self.fault_region in tect_regions.key_list:
            print self.fault_region
            raise ValueError(
                'Unspecified tectonic region - I told you so bitches!')
        
        if not self.shear_modulus:
            self.shear_modulus = \
                tect_regions.tect_reg[self.fault_region]['Shear_Modulus']
            # Ultimately we would like to support multiple values in a logic 
            # tree (TODO)- but for now just take highest weighted value
            if isinstance(self.shear_modulus, list) and \
                len(self.shear_modulus) > 1:
                max_loc = np.argmax(np.array(
                    tect_regions.tect_reg[self.fault_region]\
                    ['Shear_Modulus_Weight']))
                self.shear_modulus = self.shear_modulus[max_loc]

        if not self.disp_length_ratio:
            self.disp_length_ratio = \
                tect_regions.tect_reg[self.fault_region]['Strain_Drop_Values']
            # Ultimately we would like to support multiple values in a logic 
            # tree (TODO)- but for now just take highest weighted value
            if isinstance(self.disp_length_ratio, list) and \
                len(self.disp_length_ratio) > 1:
                max_loc = np.argmax(np.array(
                    tect_regions.tect_reg[self.fault_region]\
                    ['Strain_Drop_Weight']))
                self.disp_length_ratio = self.disp_length_ratio[max_loc]
            
        if not self.msr:
            self.msr = \
                tect_regions.tect_reg[self.fault_region]\
                ['MSR_Models']
            # Ultimately we would like to support multiple values in a logic 
            # tree (TODO)- but for now just take highest weighted value
            if isinstance(self.msr, list) and len(self.msr) > 1:
                max_loc = np.argmax(np.array(
                    tect_regions.tect_reg[self.fault_region]\
                    ['MSR_Models_Weight']))
                self.msr = self.msr[max_loc]
                   
 
    def _parse_aseismic_slip(self, aseismic_input):
        '''Parses the aseismic slip information'''
        if isinstance(aseismic_input, list) and len(aseismic_input) == 2:
            self.aseismic = aseismic_input[0]
            self.aseismic_completeness = aseismic_input[1]
        else:
            self.aseismic = float(aseismic_input)
        

    def _parse_fault_slip(self, slip_input):
        '''Parses the fault slip data into the correct format'''
        if isinstance(slip_input, list) and \
                len(slip_input) == 4:
            self.slip = slip_input[0]
            self.slip_min = slip_input[1]
            self.slip_max = slip_input[2]
            self.slip_quality = slip_input[3]
        else:
            self.slip = slip_input

    def get_recurrence_model(self, fault_dict):
        '''Calculates MFD for the fault'''
        #Temporary: Fix msr to WC1994()
        self.msr = WC1994()

        # Parse recurrence_params
        self.mfd_model = {'AndersonLucoType1': AndersonLucoType1(),
                          'AndersonLucoType2': AndersonLucoType2(),
                          'AndersonLucoType3': AndersonLucoType3(),
                          'Characteristic': Characteristic(),
                          'YoungsCoppersmith': YoungsCoppersmith()}

        if __name__ == "__main__":
            self.mfd_choice = self.mfd_type

        # Reduce slip rate by using proportion of aseismic slip
        seismic_slip = self.slip * (1.0 - self.aseismic) 
        for mfd_config in fault_dict['MFD_Model']:
           mfd = self.mfd_model[mfd_config['Model_Type']]
           mfd.setUp(mfd_config)
           mfd.get_mmax(mfd_config, self.msr, self.rake, 
                        self.fault_geometry.area)

           if 'AndersonLuco' in mfd_config['Model_Type']:
               mfd.get_mfd(seismic_slip, self.disp_length_ratio, 
                           self.shear_modulus, self.fault_geometry.width)
           else:
               mfd.get_mfd(seismic_slip, self.shear_modulus, 
                           self.fault_geometry.area)

           self.recurrence_model.append(mfd)

