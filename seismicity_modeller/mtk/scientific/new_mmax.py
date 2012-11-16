#!/usr/bin/env python

'''Module to implement statistically interpreted maximum magnitude
estimation algorithms'''

import abc
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy.stats.mstats import mquantiles
from scipy.integrate import quadrature



def _get_observed_mmax(catalogue, config):
    '''Check see if observed mmax values are input, if not then take
    from the catalogue'''

    if not config['input_mmax']:
        # If maxmag is False then maxmag is maximum from magnitude list
        max_location = np.argmax(catalogue['magnitude'])
        obsmax = catalogue['magnitude'][max_location]
        obsmaxsig = catalogue['sigmaMagnitude'][max_location]
    else:
        obsmaxsig = config['input_mmax_uncertainty']
        obsmax = config['input_mmax']
    return obsmax, obsmaxsig

def _get_magnitude_vector_properties(catalogue, config):
    '''If an input minimum magnitude is given then consider catalogue
    only above the minimum magnitude - returns corresponding properties'''
    
    if config['input_mmin']:
        neq = np.float(np.sum(catalogue['magnitude'] >= 
                              config['input_mmin'] - 1.E-7))
        mmin = config['input_mmin']
    else:
        neq = np.float(len(catalogue['magnitude'])
        mmin = np.min(catalogue['magnitude'])
    return neq, mmin


class InstrumentalMMax(object):
    '''Class for estimation of mmax from instrumental seismicity'''
    def __init__(self):
        self.mmax_master = {'KijkoSellevol': KijkoSellevolFixedb(),
                            'KijkoSellevolBayes': KijkoSellevolUncertainb(), 
                            'KijkoNPG': KijkoNonParametricGaussian(),
                            'CumulativeMoment': CumulativeMoment()}

    def analyse(self, catalogue, config):
        '''Calculate mmax'''
        mmax, sigma_mmax = self.mmax_master[config['algorithm']].get_mmax(
            catalogue,
            config)
        return mmax, sigma_mmax


class MaximumMagnitude(object):
    '''Base class for implementation of maximum magnitude calculators'''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_mmax(self, catalogue, config, **kwargs)
        '''Calculates maximum magnitude from an input catalogue
        :param catalogue: Catalogue dictionary (as defined in mtk_catalogue)
        :type catalogue: dictionary
        :param config: configuration settings containing the following - 
                       'algorithm' - Choice of algorithm - 
                                     'Kijko Sellevol Fixed b'
                                     'Kijko Sellevol Uncertain b'
                                     'Kijko Non-Parametric'
                                     'Cumulative Moment'
                       'input_mmax' - Input maximum magnitude (if false will
                                      take from catalogue'
                       'input_mmax_uncertainty' - As above for mmax uncertainty
                       'input mmin' - Input mininum magnitude
                       'maximum_iterations' - Maximum number of iterations
                       'tolerance' - Parameter for stabilising iteration
                       'number_samples' - Number of samples for Kijko NPG
                       'number_bootstraps' - Number of bootstraps for 
                                             Cumulative Moment
                       'b-value' - Input b-value
                       'sigma b' - Input b-value uncertainty'
                       'number_earthquakes' - Number of earthquakes (for 
                                              Kijko Non-Parametric)
                                             '''

class KijkoSellevolFixedb(MaximumMagnitude):
    '''Implements Kijko and Sellevol estimator with a fixed b-value'''
    def get_mmax(self, catalogue, config)
        '''Calculate mmax
        
        :return: **mmax** Maximum magnitude and **mmax_sig** corresponding
                    uncertainty
        :rtype: Float
        '''

        obsmax, obsmaxsig = _get_observed_mmax(catalogue, config)
        
        mmin = config['input_mmin']
        beta = config['b-value'] * np.log(10.)
       
        neq, mmin = _get_magnitude_vector_properties(catalogue, config)
        
        mmax = deepcopy(obsmax) 
        d_t = np.inf
        iterator = 0
        while d_t > config['tolerance']:
            delta = quadrature(_ks_intfunc, mmin, mmax, 
                               args = [neq, mmax, mmin, beta])[0]
            #print mmin, neq, delta, mmax
            tmmax = obsmax + delta
            d_t = np.abs(tmmax - mmax)
            mmax = deepcopy(tmmax)
            iterator += 1
            if iterator > config['maximum_iterations']:
                print 'Kijko-Sellevol estimator reached maximum # of iterations'
                d_t = -np.inf
        return mmax, np.sqrt(obsmaxsig ** 2. + delta ** 2.)      
    
    def check_config(self, config):
        '''Check config file inputs'''
        essential_keys = ['input_mmin', 'b-value']
        for key in essential_keys:
            if not key in config.keys():
                raise ValueError('For KijkoSellevol the key %s needs to be set'
                                 'in the configuation' % key)
        if 'tolerance' not in config.keys() or not config['tolerance']:
            config['tolerance'] = 1E-5

        if 'maximum_iterations' not in config.keys() \
            or not config['maximum_iterations']:
            config['maximum_iterations'] = 1000
        return config

    def _ks_intfunc(mval, neq, mmax, mmin, beta):
        '''Integral function inside Kijko-Sellevol estimator 
        (Eq. 6 in Kijko, 2004)
        :param mval: 
        :type mval:
        :param neq: Number of earthquakes
        :type neq:
        :param mmax: Maximum Magnitude
        :type mmax: Float
        :param mmin: Minimum Magnitude
        :type mmin: Float
        :param beta: Beta-value of the distribution
        :type beta: Float
        :return func1: Integrand of Kijko-Sellevol estimator
        :type func1: Float
        '''
        func1 = 1. - np.exp(-beta * (mval - mmin)) 
        func1 = (func1 / (1. - np.exp(-beta * (mmax - mmin)))) ** neq
        return func1


class KijkoSellevolUncertainb(MaximumMagnitude):
    '''Class to implement Kijko & Sellevol Bayesian estimator of Mmax, with
    uncertain b-value'''

    def get_mmax(self, catalogue, config):
        '''Calculate maximum magnitude
        :return: **mmax** Maximum magnitude and **mmax_sig** corresponding
                    uncertainty
        :rtype: Float
        '''
        obsmax, obsmaxsig = _get_observed_mmax(catalogue, config)
        
        beta = config['b-value'] * np.log(10.)
        sigbeta = config['sigma-b'] * np.log(10.)

        neq, mmin = _get_magnitude_vector_properties(catalogue, config)

        pval = beta / (sigbeta ** 2.)
        qval = (beta / sigbeta) ** 2.
        
        mmax = deepcopy(obsmax) 
        d_t = np.inf
        iterator = 0
        while d_t > config['tolerance']:
            rval = pval / (pval + mmax - mmin)
            ldelt = (1. / (1. - (rval ** qval))) ** neq
            delta = ldelt * quadrature(_ksb_intfunc, mmin, mmax, 
                                       args = [neq, mmin, pval, qval])[0]

            tmmax = obsmax + delta
            d_t = np.abs(tmmax - mmax)
            mmax = deepcopy(tmmax)
            iterator += 1
            if iterator > config['max_iterations']:
                print 'Kijko-Sellevol-Bayes estimator reached'
                print 'maximum # of iterations'
                d_t = -np.inf
        
        return mmax, np.sqrt(obsmaxsig ** 2. + delta ** 2.)        

    def _ksb_intfunc(mval, neq, mmin, pval, qval):
        '''Integral function inside Kijko-Sellevol-Bayes estimator
        (part of Eq. 10 in Kijko, 2004 - section 3.2)
        :param mval: Magnitude
        :type mval: Float
        :param neq: Number of Earthquakes
        :type neq: Float
        :param mmin: Minimum Magnitude
        :type mmin: Float
        :param pval: p-value (see Kijko, 2004 - section 3.2)
        :type pval: Float
        :param qval: q-value (see Kijki, 2004 - section 3.2)
        :type qval: Float
        :return func1: Output of function integrand
        :rtype func1: Float
        '''
        func1 = (1. - ((pval / (pval + mval - mmin)) ** qval)) ** neq
        return func1

    def check_config(self, config):
        '''Check config file inputs'''
        essential_keys = ['input_mmin', 'b-value', 'sigma-b']
        for key in essential_keys:
            if not key in config.keys():
                raise ValueError('For KijkoSellevolBayes the key %s needs to ' 
                    'be set in the configuation' % key)
        if 'tolerance' not in config.keys() or not config['tolerance']:
            config['tolerance'] = 1E-5

        if 'maximum_iterations' not in config.keys() \
            or not config['maximum_iterations']:
            config['maximum_iterations'] = 1000
        return config

class KijkoNonParametricGaussian(MaximumMagnitude):
    '''Class to implement non-parametric Gaussian methodology of Kijko (2004)
    '''

    def get_mmax(self, catalogue, config):
        '''Calculates maximum magnitude
        
        '''
        obsmax, obsmaxsig = _get_observed_mmax(catalogue, config)
        
        # Find number_eqs largest events
        if np.shape(catalogue['magnitude'])[0] <= config['number_earthquakes']:
            # Catalogue smaller than number of required events
            mag = deepcopy(catalogue['magnitude'])
            neq = np.float(np.shape(mag)[0])
        else:
            # Select number_eqs largest events
            mag = np.sort(catalogue['magnitude'], kind='quicksort')
            mag = mag[-config['number_earthquakes']:]
        
        mmin = np.min(mag)   
        # Get smoothing factor
        hfact = self.h_smooth(mag)
        mmax = deepcopy(obsmax)
        d_t = np.inf
        iterator = 0
        while d_t > config['tolerance']:
            # Generate exponentially spaced samples
            magval = np.log(np.hstack([np.exp(mmin) +  np.arange(0., 
                config['number_samples'] - 1., 1.) * ((np.exp(mmax) - 
                np.exp(mmin)) / (config['number_samples'] - 1.)), 
                np.exp(mmax)]))

            # Evaluate integral function using Simpson's method
            delta = _kijko_npg_intfunc_simps(magval, mag, hfact, neq)
            
            tmmax = obsmax + delta
            d_t  = np.abs(tmmax - mmax)
            mmax = deepcopy(tmmax)
            iterator += 1
            if iterator > config['maximum iterations']:
                print 'Kijko-Non-Parametric Gaussian estimator reached'
                print 'maximum # of iterations'
                d_t = -np.inf
        return mmax, np.sqrt(obsmaxsig ** 2. + delta ** 2.)
    
    
    def check_config(self, config):
        '''Check config file inputs'''
        essential_keys = ['input_mmin', 'input_mmax', 'number_earthquakes']
        for key in essential_keys:
            if not key in config.keys():
                raise ValueError('For KijkoSellevolBayes the key %s needs to ' 
                    'be set in the configuation' % key)
        if 'tolerance' not in config.keys() or not config['tolerance']:
            config['tolerance'] = 1E-5

        if 'maximum_iterations' not in config.keys() \
            or not config['maximum_iterations']:
            config['maximum_iterations'] = 1000


        if 'number_samples' not in config.keys() \
            or not config['number_samples']:
            config['number_samples'] = 1000

        return config


    def h_smooth(self, mag):
        '''Function to calculate smoothing coefficient (h) for Gaussian
        Kernel estimation - based on Silverman (1986) formula
        :param mag: Magnitude vector
        :type mag: numpy.ndarray
        :return hfact: Smoothing coefficient (h)
        :rtype hfact: Float
        '''
        neq = np.float(np.shape(mag)[0])
        
        # Calculate inter-quartile range
        qtiles = mquantiles(mag, prob = [0.25, 0.75])
        iqr = qtiles[1] - qtiles[0]
        hfact = 0.9 * np.min([np.std, iqr / 1.34]) * (neq ** (-1. / 5.))
        # Round h to 2 dp
        hfact = np.round(100. * hfact) / 100.
        return hfact

    def _gauss_cdf_hastings(xval, barx = 0.0, sigx = 1.0):
        '''Function to implement Hasting's approximation of the normalised
        cumulative normal function - this is taken from Kijko's own code
        so I don't really know why this is here!!!!!
        :param xval: x variate
        :type xval: Float or numpy.ndarray
        :param barx: Mean of the distribution
        :type barx: Float
        :param sigx: Standard Deviation
        :type sigx: Float
        :return yval: Gaussian Cumulative Distribution
        :rtype yval: Float
        '''
        x_norm = (xval - barx) / sigx

        # Fixed distribution co-efficients
        a_1 = 0.196854
        a_2 = -0.115194
        a_3 = 0.000344
        a_4 = 0.019527
        x_a = np.abs(x_norm)
        yval = 1.0 - 0.5 * (1. + a_1 * x_a + (a_2 * (x_a ** 2.)) + 
                            (a_3 * (x_a ** 3.)) + (a_4 * (x_a ** 4.))) ** (-4.)
        # To deal with precision errors for tail ends
        yval[x_norm < -5.] = 0.
        yval[x_norm > 5.] = 1.
        # Finally to normalise
        yval[x_norm < 0.] = 1. - yval[x_norm < 0.]
        return yval
        
    def _kijko_npg_intfunc_simps(mval, mag, hfact, neq):
        '''Integral function for non-parametric Gaussuan assuming that
        Simpson's rule has been invoked for exponentially spaced samples
        :param mval: Target Magnitude
        :type mval: Float
        :param mag: Observed Magnitude values
        :type mag: numpy.ndarray
        :param hfact: Smoothing coefficient (output of h_smooth)
        :type hfact: Float
        :param neq: Number of earthquakes
        :type neq: Float
        :return intfunc: Integral of non-Parametric Gaussian function
        :rtype intfunc: Float
        '''
        nmval = np.shape(mval)[0]
        mmin = np.min(mval)
        mmax = np.max(mval)

        cdf_func = np.zeros(nmval)
        for ival, target_mag in enumerate(mval):
            # Calculate normalised magnitudes
            p_min = _gauss_cdf_hastings((mmin - mag) / hfact)
            p_max = _gauss_cdf_hastings((mmax - mag) / hfact)
            p_mag = _gauss_cdf_hastings((target_mag - mag) / hfact)
            cdf_func[ival] = ((np.sum(p_mag - p_min)) / 
                              (np.sum(p_max - p_min))) ** neq
        # Now to perform integration via mid-point rule
        intfunc = 0.5 * cdf_func[0] * (mval[1] - mval[0])
        for iloc in range(0,  nmval-1):
            intfunc = intfunc + (0.5 * cdf_func[iloc] * (mval[iloc + 1] -
                                                         mval[iloc - 1]))

        intfunc = intfunc + (0.5 * cdf_func[-1] * (mval[-1] - mval[-2]))

        return intfunc



class CumulativeMoment(MaximumMagnitude):
    '''Class to implement the bootstrapped cumulative moment estimator of
    maximum magnitude. Adapted by G. Weatherill from the Cumulative Strain 
    Energy approach originally suggested by Makropoulos & Burton (1983)'''

    def get_mmax(self, catalogue, config):
        '''Calculates Mmax with bootstrap wrapper
        '''
        mag = catalogue['magnitude']
        sigma_m = catalogue['sigmaMagnitude']
        neq = np.shape(mag)[0]
        mmax_samp = np.zeros(config['number_bootstraps', dtype=float)
        for iloc in range(0, config['number_bootstraps']):
            mw_sample = mag + sigma_m * np.random.normal(0, 1, neq)
            mmax_samp[iloc] = cumulative_moment(catalogue['year'], mw_sample)
        
        # Return mean and standard deviation of sample
        return np.mean(mmax_samp), np.std(mmax_samp, ddof = 1)


    def cumulative_moment(self, year, mag, iplot = False):
        '''Calculation of Mmax using aCumulative Moment approach, adapted from
        the cumulative strain energy method of Makropoulos & Burton (1983)
        :param year: Year of Earthquake
        :type year: numpy.ndarray
        :param mag: Magnitude of Earthquake
        :type mag: numpy.ndarray
        :keyword iplot: Include cumulative moment plot
        :type iplot: Boolean
        :return mmax: Returns Maximum Magnitude
        :rtype mmax: Float
        '''
        # Calculate seismic moment
        m_o =  10. ** (9.05 + 1.5 * mag)
        year_range = np.arange(np.min(year), np.max(year) + 1, 1)
        nyr = np.float(np.shape(year_range)[0])
        morate = np.zeros(nyr, dtype=float)
        # Get moment release per year
        for loc, tyr in enumerate(year_range):
            idx = np.abs(year - tyr) < 1E-5
            if np.sum(idx) > 0:
                # Some moment release in that year
                morate[loc] = np.sum(m_o[idx])
        ave_morate = np.sum(morate) / nyr

        # Average moment rate vector
        exp_morate = np.cumsum(ave_morate * np.ones(nyr))
        modiff = np.abs(np.max(np.cumsum(morate) - exp_morate)) + \
                        np.abs(np.min(np.cumsum(morate) - exp_morate))
        # Return back to Mw
        mmax = (2./ 3.) * (np.log10(modiff) - 9.05)
        if iplot:
            plt.step(year_range, np.cumsum(morate), 'b-', linewidth = 2)
            plt.plot(year_range, exp_morate, 'r-', linewidth = 2)
            # Get offsets
            upper_morate = exp_morate + (np.max(np.cumsum(morate) - exp_morate))
            lower_morate = exp_morate + (np.min(np.cumsum(morate) - exp_morate))
            plt.plot(year_range, upper_morate, 'r--', linewidth = 1)
            plt.plot(year_range, lower_morate, 'r--', linewidth = 1)
            plt.axis([np.min(year), np.max(year), 0.0, np.sum(morate)])
            plt.show()
        return mmax
    
    def cum_mo_uncertainty(year, mag, sigma_m, nbootstrap = 1000):
        '''Function to calculate mmax wth uncertainty using the cumulative 
        moment formulation (adapted from Makropoulos & Burton, 1983) with 
        bootstrapping
        :param year: Year of earthquake
        :type year: numpy.ndarray
        :param mag: Magnitude of earthquake
        :type mag: numpy.ndarray
        :param sigma_m: Uncertainties on the magnitudes
        :type sigma_m: numpy.ndarray
        :keyword nbootstrap: Number of samples for bootstrapping
        :type nbootstrap: Integer
        :return: **mmax** Maximum magnitude and **mmax_sig** corresponding
                 uncertainty
        :rtype: Float          
        '''
        
        neq = np.shape(mag)[0]
        mmax_samp = np.zeros(nbootstrap)
        for i in range(0, nbootstrap):
            mw_sample = mag + sigma_m * np.random.normal(0, 1, neq)
            mmax_samp[i] = cumulative_moment(year, mw_sample)
        
        # Return mean and standard deviation of sample
        return np.mean(mmax_samp), np.std(mmax_samp, ddof = 1)
