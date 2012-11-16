#!/usr/bin/env/python

'''Wrapper file to call the fault calculator activity class'''

import os, sys
import new_fault_class_4 as nfc

input_file = sys.argv[1]
output_file = sys.argv[2]

fault_model = nfc.FaultModel()
fault_model.read_fault_input_model(input_file)
fault_model.get_activity_rate(output_file)
