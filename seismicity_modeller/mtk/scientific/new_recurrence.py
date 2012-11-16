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
which implement recurrence algorithms. Implemented
algorithms are:

* Recurrence analysis (Weichert, MLE)
"""

import abc
import numpy as np
from math import sqrt

def recurrence_table(mag, dmag, year):
    """
    Table of recurrence statistics for each magnitude
    [Magnitude, Number of Observations, Cumulative Number
    of Observations >= M, Number of Observations
    (normalised to annual value), Cumulative Number of
    Observations (normalised to annual value)]
    Counts number and cumulative number of occurrences of
    each magnitude in catalogue

    :param mag: catalog matrix magnitude column
    :type mag: numpy.ndarray
    :param dmag: magnitude interval
    :type dmag: numpy.ndarray
    :param year: catalog matrix year column
    :type year: numpy.ndarray
    :returns: recurrence table
    :rtype: numpy.ndarray
    """

    # Define magnitude vectors
    num_year = np.max(year) - np.min(year) + 1.
    upper_m = np.max(np.ceil(10.0 * mag) / 10.0)
    lower_m = np.min(np.floor(10.0 * mag) / 10.0)
    mag_range = np.arange(lower_m, upper_m + (1.5 * dmag), dmag)
    mval = mag_range[:-1] + (dmag / 2.0)
    # Find number of earthquakes inside range
    number_obs = np.histogram(mag, mag_range)[0]
    number_rows = np.shape(number_obs)[0]
    # Cumulative number of events
    n_c = np.zeros((number_rows, 1))
    i = 0
    while i < number_rows:
        n_c[i] = np.sum(number_obs[i:], axis=0)
        i += 1

    # Normalise to Annual Rate
    number_obs_annual = number_obs / num_year
    n_c_annual = n_c / num_year

    rec_table = np.column_stack([mval, number_obs, n_c, number_obs_annual,
                                 n_c_annual])

    return rec_table


def b_max_likelihood(mval, number_obs, dmag=0.1, m_c=0.0):
    """
    Calculation of b-value and its uncertainty for a given catalogue,
    using the maximum likelihood method of Aki (1965), with a correction
    for discrete bin width (Bender, 1983).

    :param mval: array of reference magnitudes
                 (column 0 from recurrence table)
    :type mval: numpy.ndarray
    :param number_obs: number of observations in magnitude bin
                       (column 1 from recurrence table)
    :type number_obs: numpy.ndarray
    :keyword dmag: magnitude interval
    :type dmag: positive float
    :keyword m_c: completeness magnitude
    :type m_c: float
    :returns: bvalue and sigma_b
    :rtype: float
    """

    # Exclude data below Mc
    id0 = mval >= m_c
    mval = mval[id0]
    number_obs = number_obs[id0]
    # Get Number of events, minimum magnitude and mean magnitude
    neq = np.sum(number_obs)
    m_min = np.min(mval)
    m_ave = np.sum(mval * number_obs) / neq
    # Calculate b-value
    bval = np.log10(np.exp(1.0)) / (m_ave - m_min + (dmag / 2.))
    # Calculate sigma b from Bender estimator
    sigma_b = np.sum(number_obs * ((mval - m_ave) ** 2.0)) / (neq * (neq - 1))
    sigma_b = np.log(10.) * (bval ** 2.0) * np.sqrt(sigma_b)
    return bval, sigma_b

def _input_checks(catalogue, config, completeness):
    '''Performs the basic set of input checks on the data'''

    if isinstance(completeness, np.ndarray):
        # completeness table is a numpy array (i.e. [year, magnitude])
        if np.shape(completeness)[1] != 2:
            raise ValueError('Completeness Table incorrectly configured')
        else:
            cmag = completeness[:, 1]
            ctime = completeness[:, 0]
    elif isinstance(completeness, float):
        # Completeness corresponds to a single magnitude (i.e. applies to
        # the entire catalogue)
        cmag = np.array(completeness)
        ctime = np.array(np.min(catalogue['year']))
    else:
        # Everything is valid - i.e. no completeness magnitude
        cmag = np.array(np.min(catalogue['magnitude']))
        ctime = np.array(np.min(catalogue['year']))
     
    #Set reference magnitude - if not in config then default to M = 0.
    if not config['reference_magnitude']:
        ref_mag = 0.0
    else:
        ref_mag = config['reference_magnitude']

    if not config['magnitude_interval']:
        dmag = 0.1
    else:
        dmag = config['magnitude_interval']
    
    return cmag, ctime, ref_mag, dmag


class Recurrence(object):
    '''Master recurrence class'''
    def __init__(self):
        '''Initialise'''
        self.recurrence_master = {'bMaximumLikelihood': bMaximumLikelihood(), 
                                  'KijkoSmit': KijkoSmit(),
                                  'Weichert': Weichert()}
        self.mfd_parameters = {'rate': None, 'sigma_rate': None, 
                               'bvalue': None, 'sigma_bvalue': None, 
                               'corner_magnitude': None, 'sigma_mc': None}

    def run_recurrence(self, catalogue, config, completeness=None):
        '''Execute the recurrence calculations'''
        bvalue, sigma_bvalue, rate, sigma_rate = self.recurrence_master[
            config['algorithm']].calculate(catalogue, config, completeness)
        self.mfd_parameters['rate'] = rate
        self.mfd_parameters['sigma_rate'] = sigma_rate
        self.mfd_parameters['bvalue'] = bvalue
        self.mfd_parameters['sigma_bvalue'] = sigma_bvalue
        
        return self.mfd_parameters


class SeismicityRecurrence(object):
    '''Implements recurrence calculations for instrumental seismicity'''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def calculate(self, catalogue, config, completeness=None):
        '''Implements recurrence calculation'''
        return

class bMaximumLikelihood(SeismicityRecurrence):
    '''Implements maximum likelihood calculations taking into account time
    variation in completeness'''
   
    def calculate(self, catalogue, config, completeness=None):
        '''Calculates recurrence parameters a_value and b_value, and their 
        respective uncertainties
        :param catalogue: Earthquake Catalogue (dictionary)
        :type catalogue: dictionary
        :param config: Configuration file
        :type config: dictionary
        :param completeness: Completeness table 
        :type completeness: List/numpy.ndarray or float
        
        '''

        # Input checks
        cmag, ctime, ref_mag, dmag = _input_checks(catalogue, config,
                                                    completeness)

        if not config['Average Type'] in ['Weighted','Harmonic']:
            raise ValueError('Average type not recognised in bMaxLiklihood!')

        ival = 0
        mag_eq_tolerance = 1E-5
        while ival < np.shape(ctime)[0]:
            id0 = np.abs(ctime - ctime[ival]) < mag_eq_tolerance
            m_c = np.min(cmag[id0])
            # Find events later than cut-off year, and with magnitude
            # greater than or equal to the corresponding completeness magnitude.
            # m_c - mag_eq_tolerance is required to correct floating point
            # differences.
            id1 = np.logical_and(catalogue['year'] >= ctime[ival],
                catalogue['magnitude'] >= (m_c - mag_eq_tolerance))
            nyr = np.float(np.max(catalogue['year'][id1])) - ctime[ival] + 1.
            # Get a- and b- value for the selected events
            temp_rec_table = recurrence_table(catalogue['magnitude'][id1], 
                                              dmag, 
                                              catalogue['year'][id1])
            bval, sigma_b = b_max_likelihood(temp_rec_table[:, 0],
                                             temp_rec_table[:, 1], dmag, m_c)

            aval = np.log10(np.float(np.sum(id1)) / nyr) + bval * m_c
            sigma_a = np.abs(np.log10(np.float(np.sum(id1)) / nyr) +
                (bval + sigma_b) * ref_mag - aval)

            # Calculate reference rate
            rate = 10.0 ** (aval - bval * ref_mag)
            sigrate = 10.0 ** ((aval + sigma_a) - (bval * ref_mag) -
                np.log10(rate))
            if ival == 0:
                gr_pars = np.array([np.hstack([bval, sigma_b, rate, sigrate])])
                neq = np.sum(id1)  # Number of events
            else:
                gr_pars = np.vstack([gr_pars, np.hstack([bval, sigma_b, rate,
                                                         sigrate])])
                neq = np.hstack([neq, np.sum(id1)])
            ival = ival + np.sum(id0)
        bval, sigma_b, aval, sigma_a = self._average_parameters(gr_pars, neq, 
            config['Average Type'])
        if not config['reference_magnitude']:
            #print aval, sigma_a, bval, sigma_b
            d_aval = aval - sigma_a
            aval = np.log10(aval)
            sigma_a = aval - np.log10(d_aval)
            #print aval, sigma_a
        return bval, sigma_b, aval, sigma_a
    
    def _average_parameters(self, gr_params, neq, average_type='Weighted'):
        '''Calculates the average of a set of Gutenberg-Richter parameters
        depending on the average type
        :param gr_params: Gutenberg-Richter parameters [b, sigma_b, a, sigma_a]
        :type gr_params: numpy.ndarray
        '''
        if np.shape(gr_params)[0] != np.shape(neq)[0]:
            raise ValueError('Number of weights does not correspond'
                             ' to number of parameters')
        
        #neq = neq.astype(float)
        if 'Harmonic' in average_type:
            average_parameters = self._harmonic_mean(gr_params, neq)
        else:
            average_parameters = self._weighted_mean(gr_params, neq)
        bval = average_parameters[0]
        sigma_b = average_parameters[1]
        aval = average_parameters[2]
        sigma_a = average_parameters[3]
        return bval, sigma_b, aval, sigma_a 
        

    def _weighted_mean(self, parameters, neq):
        '''Simple weighted mean'''
        weight = neq.astype(float) / np.sum(neq)
        if np.shape(parameters)[0] != np.shape(weight)[0]:
            raise ValueError('Parameter vector not same shape as weights')
        else:
            average_value = np.zeros(np.shape(parameters)[1], dtype=float)
            for iloc in range(0, np.shape(parameters)[1]):
                average_value[iloc] = np.sum(parameters[:, iloc] * weight)
        return average_value

    def _harmonic_mean(self, parameters, neq):
        '''Harmonic mean'''
        weight = neq.astype(float) / np.sum(neq)
        if np.shape(parameters)[0] != np.shape(weight)[0]:
                raise ValueError('Parameter vector not same shape as weights')
        else:
            average_value = np.zeros(np.shape(parameters)[1], dtype=float)
            for iloc in range(0, np.shape(parameters)[1]):
                average_value[iloc] = 1. / np.sum(
                    (weight * (1. / parameters[:, iloc])))
        return average_value


class KijkoSmit(SeismicityRecurrence):
    '''Class to Implement the Kijko & Smit (2012) algorithm for estimation
    of a- and b-value'''
    def calculate(self, catalogue, config, completeness=None):
        '''Main function to calculate the a- and b-value'''
        # Input checks
        cmag, ctime, ref_mag, dmag = _input_checks(catalogue, config, 
                                                   completeness)
        ival = 0
        mag_eq_tolerance = 1E-5
        number_intervals = np.shape(ctime)[0]
        b_est = np.zeros(number_intervals, dtype=float)
        #sigma_b_est = np.zeros(number_intervals, dtype=float)
        neq = np.zeros(number_intervals, dtype=float)
        nyr = np.zeros(number_intervals, dtype=float)
        while ival < number_intervals:
            id0 = np.abs(ctime - ctime[ival]) < mag_eq_tolerance
            m_c = np.min(cmag[id0])
            # Find events later than cut-off year, and with magnitude
            # greater than or equal to the corresponding completeness magnitude.
            # m_c - mag_eq_tolerance is required to correct floating point
            # differences.
            id1 = np.logical_and(catalogue['year'] >= ctime[ival],
                catalogue['magnitude'] >= (m_c - mag_eq_tolerance))
            nyr[ival] = np.float(np.max(catalogue['year'][id1]) -
                                 np.min(catalogue['year'][id1]) + 1)
            neq[ival] = np.sum(id1)
            # Get a- and b- value for the selected events
            temp_rec_table = recurrence_table(catalogue['magnitude'][id1], 
                                              dmag, 
                                              catalogue['year'][id1])
            b_est[ival]= b_max_likelihood(temp_rec_table[:, 0], 
                                          temp_rec_table[:, 1], 
                                          dmag, m_c)[0]
            ival += 1

        total_neq = np.float(np.sum(neq))
        #print b_est, neq
        bval = self._harmonic_mean(b_est, neq)
        #sigma_b_est = self.harmonic_mean(sigma_b_est, neq)
        sigma_b = bval / sqrt(total_neq)
        #print bval, total_neq, nyr, cmag, ref_mag
        aval = self._calculate_a_value(bval, total_neq, nyr, cmag, ref_mag)
        sigma_a = self._calculate_a_value(bval + sigma_b, total_neq, nyr, 
                                          cmag, ref_mag)
        if not config['reference_magnitude']:
            #print aval, sigma_a, bval, sigma_b
            #d_aval = sigma_a - aval
            aval = np.log10(aval)
            sigma_a = np.log10(sigma_a) - aval
            #print aval, sigma_a
        else: 
            sigma_a = sigma_a - aval
        return bval, sigma_b, aval, sigma_a 
        
    def _harmonic_mean(self, parameters, neq):
        '''Harmonic mean'''
        weight = neq.astype(float) / np.sum(neq)
        if np.shape(parameters)[0] != np.shape(weight)[0]:
            raise ValueError('Parameter vector not same shape as weights')
        else:
            average_value = np.zeros(np.shape(parameters)[0], dtype=float)
            average_value = 1. / np.sum(weight / parameters)
        return average_value

    def _calculate_a_value(self, bval, nvalue, nyr, cmag, ref_mag):
        '''Calculates the rate of events >= ref_mag using the b-value estimator
        and Eq. 10 of Kijko & Smit'''

        denominator = np.sum(nyr * np.exp(-bval * (cmag - ref_mag)))
        return nvalue / denominator

class Weichert(SeismicityRecurrence):
    '''Class to Implement Weichert Algorithm'''

    def calculate(self, catalogue, config, completeness=None):
        '''Calculates recurrence using the Weichert (1980) method'''
        # Input checks
        cmag, ctime, ref_mag, dmag = _input_checks(catalogue, config,
                                                   completeness)
        # Apply Weichert preparation
        cent_mag, t_per, n_obs = self._weichert_prep(catalogue['year'],
                                                     catalogue['magnitude'],
                                                     ctime, cmag, dmag)

        # A few more Weichert checks
        key_list = config.keys()
        if (not 'bvalue'  in key_list) or (not config['bvalue']):
            config['bvalue'] = 1.0
        if (not 'itstab' in key_list) or (not config['itstab']):
            config['itstab'] = 1E-5
        if (not 'maxiter' in key_list) or (not config['maxiter']):
            config['maxiter'] = 1000

        bval, sigma_b, aval, sigma_a = self.weichert_algorithm(t_per,
            cent_mag, n_obs, ref_mag, config['bvalue'], config['itstab'], 
            config['maxiter'])

        if not config['reference_magnitude']:
            aval = np.log10(aval)
            sigma_a = np.log10(sigma_a)
        return bval, sigma_b, aval, sigma_a

    def _weichert_prep(self, year, magnitude, ctime, cmag, dmag, dtime=1.0):
        """
        Allows to prepare table input for Weichert algorithm

        :param year: catalog matrix year column
        :type year: numpy.ndarray
        :param magnitude: catalog matrix magnitude column
        :type magnitude: numpy.ndarray
        :param ctime: year of completeness for each period
        :type ctime: numpy.ndarray
        :param cmag: completeness magnitude for each period
        :type cmag: numpy.ndarray
        :param dmag: magnitude bin size (config file)
        :type dmag: positive float
        :param dtime: time bin size from config file)
        
        :type dtime: float
        :returns: central magnitude, tper length of observation period,
                  n_obs number of events in magnitude increment
        """

        # In the case that the user defines a single value for ctime or cmag
        # that is not an array
        if not(isinstance(ctime, np.ndarray)) and not(isinstance(ctime, list)):
            ctime = np.array([ctime])
        if not(isinstance(cmag, np.ndarray)) and not(isinstance(cmag, list)):
            cmag = np.array([cmag])
        valid_events = np.ones(np.shape(year)[0], dtype = bool)
        # Remove events from catalogue below completeness intervals
        mag_eq_tolerance = dmag / 1.E7
        time_tolerance = dtime / 1.E7

        for iloc, mag in enumerate(cmag):
            index0 = np.logical_and(magnitude < (mag - mag_eq_tolerance), 
                                    year < (ctime[iloc] - time_tolerance))
            valid_events[index0] = False
      
           
        year = year[valid_events]
        magnitude = magnitude[valid_events]
        for iloc, yr in enumerate(year):
            dum = np.hstack([iloc, yr, magnitude[iloc]])

        mag_range = np.arange(np.min(magnitude) - dmag / 2., 
                              np.max(magnitude) + (2.0 * dmag), dmag)
        time_range = np.arange(np.min(year) - dtime / 2.,
                               np.max(year) + (2.0 * dtime), 
                               dtime)
        
        # Histogram data
        fullcount1 = np.histogram2d(year, magnitude, 
                                    bins = [time_range, mag_range])[0]
        n_y = np.shape(fullcount1)[1] - 1
        cent_mag = ((mag_range[:-1] + mag_range[1:]) / 2.)[:-1]
        
        n_obs = np.sum(fullcount1, axis=0)[:-1]
        t_per = np.zeros(n_y)
        for iloc, mag in enumerate(cmag):
            index0 = cent_mag > (mag - (dmag / 2. - mag_eq_tolerance))
            t_per[index0] = np.max(year) - ctime[iloc] + 1

        #cut off magnitudes below the lowest magnitude of completeness
        valid_location = np.nonzero(t_per)[0][0]
        cent_mag = cent_mag[valid_location:]
        t_per = t_per[valid_location:]
        n_obs = n_obs[valid_location:]
        
        return cent_mag, t_per, n_obs

    def weichert_algorithm(self, tper, fmag, nobs, mrate=0.0, bval=1.0, 
                           itstab=1E-5, maxiter=1000):
        """
        Weichert algorithm

        :param tper: length of observation period corresponding to magnitude
        :type tper: numpy.ndarray (float)
        :param fmag: central magnitude
        :type fmag: numpy.ndarray (float)
        :param nobs: number of events in magnitude increment
        :type nobs: numpy.ndarray (int)
        :keyword mrate: reference magnitude
        :type mrate: float
        :keyword bval: initial value for b-value
        :type beta: float
        :keyword itstab: stabilisation tolerance
        :type itstab: float
        :keyword maxiter: Maximum number of iterations
        :type maxiter: Int
        :returns: b-value, sigma_b, a-value, sigma_a
        :rtype: float
        """
        beta = bval * np.log(10.)
        d_m = fmag[1] - fmag[0]
        itbreak = 0
        snm = np.sum(nobs * fmag)
        nkount = np.sum(nobs)
        iteration = 1
        while (itbreak != 1):
            beta_exp = np.exp(-beta * fmag)
            tjexp = tper * beta_exp
            tmexp = tjexp * fmag
            sumexp = np.sum(beta_exp)
            stmex = np.sum(tmexp)
            sumtex = np.sum(tjexp)
            stm2x = np.sum(fmag * tmexp)
            dldb = stmex / sumtex
            if np.isnan(stmex) or np.isnan(sumtex):
                raise ValueError('NaN occers in Weichert iteration')

            d2ldb2 = nkount * ((dldb ** 2.0) - (stm2x / sumtex))
            dldb = (dldb * nkount) - snm
            betl = np.copy(beta)
            beta = beta - (dldb / d2ldb2)
            sigbeta = np.sqrt(-1. / d2ldb2)

            if np.abs(beta - betl) <= itstab:
                # Iteration has reached convergence
                bval = beta / np.log(10.0)
                sigb = sigbeta / np.log(10.)
                fngtm0 = nkount * (sumexp / sumtex)
                fn0 = fngtm0 * np.exp((-beta) * (fmag[0] - (d_m / 2.0)))
                stdfn0 = fn0 / np.sqrt(nkount)
                if mrate == 0.:
                    a_m = fngtm0
                    siga_m = stdfn0
                else:
                    a_m = fngtm0 * np.exp((-beta) * (mrate -
                                                    (fmag[0] - (d_m / 2.0))))
                    siga_m = a_m / np.sqrt(nkount)
                itbreak = 1
            else:
                iteration += 1
                if iteration > maxiter:
                    raise RuntimeError('Maximum Number of Iterations reached')
                continue
        return bval, sigb, a_m, siga_m
