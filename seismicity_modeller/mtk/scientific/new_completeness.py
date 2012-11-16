# -*- coding: utf-8 -*-

# Copyright (c) 2010-2012, GEM Foundation.
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.


"""
The purpose of this module is to provide functions
which implement completeness algorithms. Implemented
algorithms are:

* Stepp
* AssumedMFD
"""


import numpy as np
import abc
#import logging

#LOGGER = logging.getLogger('mt_logger')

def check_completeness_config_master(config, supported_algorithms):
    '''Checks whether the config file contains all the keys and valid
    parameters common to all completeness algorithms'''

    expected_list = ['algorithm', 'delta_m', 'delta_t', 'increment_lock']

    key_list = config.keys()
    key_list.sort()
    
    for key in key_list:
        if key in expected_list:
            if key is 'algorithm':
                assert(isinstance(config[key], str) and 
                       config[key] in supported_algorithms)
            if (key is 'delta_m') or (key is 'delta_t'):
                assert(isinstance(config[key], float) and 
                       config[key] > 0.0)
            if key is 'increment_lock':
                assert (isinstance(config[key], bool))
            expected_list.remove(key)

    if len(expected_list) > 0:
        print 'Expected keys not defined:'
        print expected_list
        raise ValueError('Not all required keys defined in config')


class Completeness(object):
    '''Class for completeness analysis'''
    def __init__(self):
        '''Initialise'''
        self.completeness_master = {'Stepp': Stepp1971()}

    def analysis(self, catalogue, config, valid_config=False):
        '''Performs the completeness analysis'''
        if not valid_config:
            self.completeness_master[config['algorithm']].\
            check_completeness_config(config, self.completeness_master.keys())

        completeness_table = \
            self.completeness_master[config['algorithm']].analyse(catalogue, 
                                                                 config)
        return completeness_table


class CatalogueCompleteness(object):
    '''Abstract base class for completeness algorithms'''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def analyse(self, catalogue, config):
        '''Base class for completeness catalogue
        :param catalogue: Catalogue dictionary originating from the 
                          mtk_catalogue class
        :type catalogue: Dictionary
        :param config: Configuraition dictionary containing the following keys- 
                       'delta_m' - magnitude window
                       'delta_t' - time window (years)
                       'tolerance' - the sensitivity parameter
                       'increment_lock'
        :type config: Dictionary'''
        return

class Stepp1971(CatalogueCompleteness):
    '''Class to implement the Stepp (1971) algorithm'''

    def analyse(self, catalogue, config):
        '''Implements the completeness analysis'''
        # Round off the magnitudes to 2 d.p
        m_w = np.around(catalogue['magnitude'], 2)
        lowm = np.floor(10. * np.min(m_w)) / 10.
        highm = np.ceil(10. * np.max(m_w)) / 10.
        # Determine magnitude bins
        mbin = np.arange(lowm, highm + config['delta_m'], config['delta_m'])
        ntb = len(mbin)
        # Determine time bins
        end_time = np.max(catalogue['year'])
        start_time = np.min(catalogue['year'])
        time_range = np.arange(config['delta_t'], end_time - start_time + 2,
                               config['delta_t'])
        n_t = len(time_range)
        t_upper_bound = end_time * np.ones(n_t)
        t_lower_bound = t_upper_bound - time_range
        t_rate = 1. / np.sqrt(time_range)  # Poisson rate

        number_obs = np.zeros((n_t, ntb - 1))
        lamda = np.zeros((n_t, ntb - 1))
        siglam = np.zeros((n_t, ntb - 1))
        ii = 0
        # count number of events catalogue and magnitude windows
        for iloc, lower_year in enumerate(t_lower_bound):
        # Select earthquakes later than or in Year[ii]
            yrchk = catalogue['year'] > lower_year
            mtmp = m_w[yrchk]
            jj = 0
            for jloc, lower_mag in enumerate(mbin[:-1]):
                #Count earthquakes in magnitude bin
                if jloc == ntb - 1:
                    number_obs[iloc, jloc] = np.sum(mtmp >= lower_mag)
                else:
                    number_obs[iloc, jloc] = np.sum(np.logical_and(
                        mtmp >= lower_mag, mtmp < mbin[jloc + 1]))

        time_diff = (np.log10(t_rate[1:]) - np.log10(t_rate[:-1]))
        time_diff = time_diff / (
            np.log10(time_range[1:]) - np.log10(time_range[:-1]))
        comp_length = np.zeros((ntb - 1, 1))
        tloc = np.zeros((ntb - 1, 1), dtype=int)
        ii = 0
        for iloc in range(0, ntb - 1):
            lamda[:, iloc] = number_obs[:, iloc] / time_range
            siglam[:, iloc] = np.sqrt(lamda[:, iloc] / time_range)
            zero_find = siglam[:, iloc] < 1E-14   # To avoid divide by zero
            siglam[zero_find, iloc] = 1E-14
            grad1 = (np.log10(siglam[1:, iloc]) - np.log10(siglam[:-1, iloc]))
            grad1 = grad1 / (np.log10(time_range[1:]) - 
                             np.log10(time_range[:-1]))
            resid1 = grad1 - time_diff
            test1 = np.abs(resid1[1:] - resid1[:-1])
            tloct = np.nonzero(test1 > config['tolerance'])[0]
            if not(np.any(tloct)):
                tloct = -1
            else:
                tloct = tloct[-1]
            if tloct < 0:
                # No location passes test
                if iloc > 0:
                    # Use previous value
                    tloc[iloc] = tloc[iloc - 1]
                else:
                    # Print warning
                    error(
                        "Fitting tolerance removed all data - change parameter")
            else:
                tloc[iloc] = tloct
            if tloct > np.max(np.shape(time_range)):
                tloc[iloc] = np.max(np.shape(time_range))

            if iloc > 0:
                # If the increasing completeness is option is set
                # and the completeness is lower than the previous value
                # then fix at previous value
                if config['increment_lock'] and (tloc[iloc] < tloc[iloc - 1]):
                    tloc[iloc] = tloc[iloc - 1]
            comp_length[iloc] = time_range[tloc[iloc]]

        completeness_table = np.column_stack(
            [end_time - comp_length, mbin[:-1].T])

        return completeness_table


    def check_completeness_config(self, config, supported_algorithms):
        '''Performs set of checks on the completeness file specific to the
        Stepp (1971) algorithm'''
        # Check the core master values
        check_completeness_config_master(config, supported_algorithms)

        # The only additional parameter needed by Stepp1971 is "tolerance"
        if not 'tolerance' in config.keys():
            raise ValueError('Tolerance parameter needed in config for Stepp')

        if config['tolerance'] < 1E-12:
            raise ValueError('Stepp Tolerance must be > 0')
        


class AssumedMFD(CatalogueCompleteness):
    '''Class to implement the set of completeness analysis methods
    requiring the assumption of an exponential magnitude frequency 
    distribution'''
    def analyse(self, catalogue, config):
        '''Implements algorithm
        :param config: Configuration dictionary containing the following 
                       fields: 
                       'mfd_function' - The name of the function to be used
                       'time_bin_size' - The size of the time window
                       'magnitude_bin_size' - The size of the magnitude window
                       'time_step'
                       'min_earthquakes' - The minimum number of earthquakes
                       'increment_lock'
        :type config: Dictionary
        
        '''
        # Initialise base class
        analysis_tools = {'Maximum Curvature': MaximumCurvature(),
                          'b-value Stability': bValueStability(),
                          'Goodness-of-fit': GoodnessOfFit()}

        time_lims = _get_input_completeness_limits(catalogue['year'], 
                                                   config['time_bin_size'], 
                                                   config['min_earthquakes'])
        n_int = np.shape(time_lims)[0]
        completeness_table = np.zeros([n_int, 2], dtype=float)
        if 'fitness_level' in config.keys():
            gof_level = config['fitness_level']
        else:
            gof_level = 90.0

        iloc = 0
        while iloc < n_int:
        #for iloc in range(0, n_int):
            isel = np.logical_and(year <= time_lims[iloc, 1], 
                                  year >= time_lims[iloc, 0])
            completeness_table[iloc, 1] = time_lims[iloc, 0]

            if np.sum(isel) < 3:
                # Array is empty
                print 'Fewer than 3 events in time period %8.0f to %8.0f' % \
                        (time_lims[iloc, 0], time_lims[iloc, 1])
                completeness_table[iloc, 0] = completeness_table[iloc - 1, 0]
                continue
                
            rec_tab = recurrence_table(catalogue['magnitude'][isel], 
                                       config['magnitude_bin_size'], 
                                       catalogue['year'][isel])
            mc = analysis_tools[config['algorithm']].get_mc(
                rec_tab, 
                delta_m = config['magnitude_bin_size'], 
                fitness_level=gof_level)
            if len(mc) > 1:
                completeness_table[iloc, 0] = mc[0]
            else:
                completeness_table[iloc, 0] = mc
            
            # See if failures have occured
            if np.isnan(completeness_table[iloc, 0]):
                if iloc == 0:
                    print 'Completeness table failed on 1st window'
                    print 'Merging windows'
                    time_lims, n_int, completeness_table = \
                        _cut_completeness_first_row(time_lims, n_int,
                        completeness_table)
                    continue
                else:
                    completeness_table[iloc, 0] = completeness_table[iloc - 1, 
                                                                     0]

            if config['increment_lock'] and (iloc > 0) and \
                (completeness_table[iloc, 0] < 
                 completeness_table[iloc - 1, 0]):
                # Mc for the current increment is lower than for the previous
                # increment. The increment_lock option has been selected so 
                #this is not permitted; hence the Mc is rendered
                completeness_table[iloc, 0] = completeness_table[iloc - 1, 0]

            iloc += 1
        return completeness_table


def _get_input_completeness_limits(year, dtime=False, mineq=False):
    '''Function to define the completeness table, given the choice of input
    settings
    :param year: Year of earthquake
    :type year: numpy.ndarray
    :keyword dtime: Tuple containing time window size (years), time increment
    :type dtime: Tuple
    :keyword mineq: Minimum number of events for searching (ignored if False)
    :type mineq: integer
    :returns time_lims: Time earliest and latest times for each completeness
                        window
    :rtype time_lims: numpy.ndarray   
     '''
    
    if dtime:
        # Use fixed time window dtime (years)
        time_lims = _get_time_window_fixed(year, dtime)
    
    elif mineq:
        # Use minimum number of earthquakes
        time_lims = _get_time_window_mineq(year, mineq)[0]
        
    else:
        # Use whole catalogue
        time_lims = np.hstack([np.min(year), np.max(year)])

    return time_lims
    
def _cut_completeness_first_row(time_lims, completeness_table):
    '''Function to remove a row from the completeness table and
     revise the time limit'''
    time_lims[i + 1, 1] = time_limes[i, 1]
    time_lims = time_lims[i+1,:]
    nint = np.shape(time_lims)[0]
    completeness_table = completeness_table[i+1,:]
    return time_lims, nint, completeness_table

def _get_time_window_fixed(year, dtime):
    '''Function to define time window vector defined by moving
    a fixed time window over an increment
    :param year: Year of event
    :type year: numpy.ndarray
    :param dtime: Tuple containing time window size (years), time increment
    :type dtime: Tuple
    :returns time_lims: Vector of start and end years for each time window
    :rtype time_lims: numpy.ndarray
    '''

    time_window = dtime[0]
    time_inc = dtime[1]
    endt = np.max(year)
    startt = endt - time_window + 1
    time_lims = np.zeros((1, 2))
    time_lims[0, 0] = startt
    time_lims[0, 1] = endt
    while np.min(time_lims[:, 0]) >= np.min(year):
        time_lims = np.vstack([time_lims, time_lims[-1, :] - time_inc])
    
    time_lims = time_lims[:-2, :]
    time_lims[-1, 0] = np.min(year)
    
    return time_lims
    
def _get_time_window_mineq(year, mineq):
    '''Function to determine time windows based in a minimum number 
    of events
    :param year: Year of event
    :type year: numpy.ndarray
    :param mineq: Minimum number of events for time window
    :type mineq: int
    :returns: **time_lims** Vector of start and end years for each time window
               and **neq** the number of events in each time window
    :rtype: numpy.ndarray
    '''
    endt = np.max(year)
    mint = np.min(year)
    year_range = np.arange(endt, mint - 1, -1)
    nyear = np.shape(year_range)[0]
    t_1 = endt
    i = 0
    time_lims = np.array([], dtype = float)
    while i < nyear:
        for j in year_range[i:]:
            tsel = np.logical_and(year <= t_1, year >= j)
            if np.sum(tsel) >= mineq:
                # Range has at least minEQ events
                if np.shape(time_lims)[0] == 0:
                    time_lims = np.hstack([j, endt])
                    neq = np.array(np.sum(tsel))
                else:
                    time_lims = np.vstack([time_lims, 
                                            np.hstack([j, t_1])])
                    neq = np.hstack([neq, np.sum(tsel)])
                t_1 = j - 1
                i += 1
                break
            else:
                i += 1
    
    return time_lims, neq

class AssumedMFDMethods(object):
    '''Abstract base class for each of the methods in the AssumedMFD class'''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_mc(self, rec_table, **kwargs):
        '''Calculates the completeness magnitude based on an input catalogue
        table'''
        return

class MaximumCurvature(AssumedMFDMethods):
    '''Implments maximum curvature method of estimating completeness - as
    described in Woessner & Wiemer (2005)'''

    def get_mc(self, rec_table, **kwargs):
        """Function to calculate the completeness of a set of data using the
        maximum curvature method
        :param rec_table: Table counting numbers of observed magnitudes (output
                           from recurrence.recurrence_table
        :type rec_table: numpy.ndarray
        :returns mc: Completeness Magnitude
        :rtype mc: Float
        """
        
        nobs = rec_table[:, 1]
        mc = rec_table[np.argmax(nobs), 0]
         
        return mc
 

class bValueStability(AssumedMFDMethods):
    '''Implments the b-value stability  method of estimating completeness - as
    described in Wiemer & Wyss (2000) and Woessner & Wiemer (2005)'''

    def get_mc(self, rec_table, **kwargs):
        """Function to calculate completeness magnitude using the b-value
        stability method (Wiemer & Wyss, 2000)
        :param rec_table: Table counting numbers of observed magnitudes (output
                           from recurrence.recurrence_table
        :type rec_table: numpy.ndarray
        :returns mc: Completeness Magnitude
        :rtype mc: Float
        """
        nmag = np.shape(rec_table)[0]
        
        dm = rec_table[1, 0] - rec_table[0, 0]
        if nmag < 5:
            print "Magnitude range too small to determine Mc via b-value stability"
            mc = np.nan
        else:
            # Get b-values
            i = 0
            bvalue = np.zeros((nmag - 2, 1))
            sigma_b = np.zeros((nmag - 2, 1))

            for iloc in range(0, nmag - 2):
                bvalue[iloc], sigma_b[iloc] = b_max_likelihood(
                     rec_table[iloc:, 0], 
                     rec_table[iloc:, 1], 
                     kwargs['delta_m'], 
                     rec_table[iloc, 0])
            
            mc_check = 0
            for iloc in range(0, nmag - 5):
                b_ave = np.sum(bvalue[iloc:iloc  + 4]) / 5
                delta_b = np.abs(b_ave - bvalue[iloc])
                if np.abs(delta_b < sigma_b[iloc]):
                    mc_check = 1
                    break
                
            
            if mc_check == 0:
                print "b-value did not stabilise:"
                print "completeness magnitude not estimated"
                mc = np.nan
            else:
                mc = rec_table[iloc, 0]
        return mc
        
    
class GoodnessOfFit(AssumedMFDMethods):
    
    '''Implements Goodness of Fit algorithm for estimation of completeness 
    magnitude. As described in Wiemer & Wyss (2000)'''
    def get_mc(self, rec_table, **kwargs):
        '''Function to calculate Mc via the Goodness of Fit method 
        (Wiemer & Wyss, 2000)
        :param rec_table: Table counting numbers of observed magnitudes (output
                           from recurrence.recurrence_table
        :type rec_table: numpy.ndarray
        :keyword fit_level: Confidence level for the Goodness of Fit test
        :type fit_level: Float
        :returns mc: Completeness Magnitude
        :rtype mc: Float
        '''
        if not kwargs['fitness_level']:
            fit_level = 90.0
        else:
            fit_level = kwargs['fitness_level']

        nmags = np.shape(rec_table)[0]
        temp_mc = np.zeros(nmags - 2, dtype=float)
        d_m = kwargs['delta_m']
        rfit = np.zeros(nmags - 2, dtype=float)

        for iloc in range(0, nmags - 2):
            ref_mags = rec_table[iloc:, 0]
            ref_n = rec_table[iloc:, 1]
            ref_nc = rec_table[iloc:, 2]
            # Calculate b-value for observations
            temp_mc[iloc] = rec_table[iloc, 0]
            bval = b_max_likelihood(ref_mags, ref_n, d_m, temp_mc[iloc])[0]
            # Calculate rate >= mc
            if ref_nc[0] == 0:
                rfit[iloc] = -np.inf
            else:
                aval = np.log10(ref_nc[0]) + bval * temp_mc[iloc]
                # Calculate expected GR distribution
                lognc_exp = aval - bval * ref_mags
                nc_exp = np.floor(10. ** (lognc_exp))
                # Calculate Goodness of fit term (R)
                rfit[iloc] = np.sum(np.abs(ref_nc - nc_exp)) / np.sum(ref_nc)
                rfit[iloc] = 100. - 100. * rfit[iloc]
        
        # Find the lowest value of mc for which rfit >= fit level
        ifit = np.nonzero(rfit > fit_level)[0]

        if np.shape(ifit)[0] == 0:
            # No value of mc reaches this fit level
            mc = np.nan
            iloc = np.argmax(rfit)
            best_tuple = np.hstack([temp_mc[iloc], rfit[iloc]])
        else:
            #
            iloc = ifit[0]
            mc = temp_mc[iloc]
            best_tuple = np.hstack([mc, rfit[iloc]])
        return [mc, best_tuple]
