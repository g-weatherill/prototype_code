#!/usr/bin/env python

'''Module to implement statistically interpreted maximum magnitude
estimation algorithms'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles
from scipy.integrate import quadrature


def ks_intfunc(mval, neq, mmax, mmin, beta):
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
    
def kijko_sellevol(mag, mag_sig, bval, tol = 1.0E-5, maxiter = 1E3,  
                   obsmax=False):
    '''Function to implement Kijko & Sellevol estimator of Mmax, assuming
    a fixed b-value. 
    :param mag: Earthquake magnitudes
    :type mag: numpy.ndarray
    :param mag_sig: Earthquake Magnitude Uncertainty
    :type mag_sig: numpy.ndarray
    :param bval: b-value
    :type bval: float
    :keyword tol: Intergral tolerance
    :type tol: Float
    :keyword maxiter: Maximum number of Iterations
    :type maxiter: Int
    :keyword obsmax: Maximum Observed Magnitude (if not in magnitude array)
                     and its corresponding uncertainty (sigma)
    :type obsmax: Tuple (float) or Boolean
    :return: **mmax** Maximum magnitude and **mmax_sig** corresponding
             uncertainty
    :rtype: Float
    '''
    if not(obsmax):
        # If maxmag is False then maxmag is maximum from magnitude list
        del(obsmax)
        obsmax = np.max(mag)
        obsmaxsig = mag_sig[np.argmax(mag)]
    else:
        obsmaxsig = obsmax[1]
        obsmax = obsmax[0]

    
    beta = bval * np.log(10.)
    neq = np.float(np.shape(mag)[0])
    
    mmax = np.copy(obsmax) 
    d_t = 1.0E8
    
    #mmin = np.min(mag)
    mmin = 4.0
    j = 1
    while d_t > tol:
        delta = quadrature(ks_intfunc, mmin, mmax, 
                           args = [neq, mmax, mmin, beta])[0]
        print mmin, neq, delta, mmax
        tmmax = obsmax + delta
        d_t = np.abs(tmmax - mmax)
        mmax = np.copy(tmmax)
        j += 1
        if j > maxiter:
            print 'Kijko-Sellevol estimator reached maximum # of iterations'
            d_t = 0.5 * tol
    mmax_sig = np.sqrt(obsmaxsig ** 2. + delta ** 2.)
    return mmax, mmax_sig        
        
    
def ksb_intfunc(mval, neq, mmin, pval, qval):
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

    
def kijko_sellevol_bayes(mag, mag_sig, bval, sigbval, tol = 1.0E-5, 
                         maxiter = 1E3, obsmax=False):
    '''Kijko-Sellevol-Bayes estimator of MMax and its associated uncertainty
    (Eq. 12 - Kijko, 2004). 
    :param mag: Earthquake magnitudes
    :type mag: numpy.ndarray
    :param mag_sig: Earthquake Magnitude Uncertainty
    :type mag_sig: numpy.ndarray
    :param bval: b-value
    :type bval: float
    :param sigbval: b-value uncertainty
    :type sigbval: float
    :keyword tol: Intergral tolerance
    :type tol: Float
    :keyword maxiter: Maximum number of Iterations
    :type maxiter: Int
    :keyword obsmax: Maximum Observed Magnitude (if not in magnitude array)
                     and its corresponding uncertainty (sigma)
    :type obsmax: Tuple (float) or Boolean
    :return: **mmax** Maximum magnitude and **mmax_sig** corresponding
             uncertainty
    :rtype: Float    
    '''
    
    if not(obsmax):
        # If maxmag is False then maxmag is maximum from magnitude list
        del(obsmax)
        obsmax = np.max(mag)
        obsmaxsig = mag_sig[np.argmax(mag)]
    else:
        obsmaxsig = obsmax[1]
        obsmax = obsmax[0]
    
    neq = np.float(np.shape(mag)[0])
    mmin = np.min(mag)
    mmax = np.copy(obsmax)
    beta = bval * np.log(10.)
    sigbeta = sigbval * np.log(10.)
    pval = beta / (sigbeta ** 2.)
    qval = (beta / sigbeta) ** 2.
    
    d_t = 1.E8
    j = 0
    while d_t > tol:
        rval = pval / (pval + mmax - mmin)
        ldelt = (1. / (1. - (rval ** qval))) ** neq
        delta = ldelt * quadrature(ksb_intfunc, mmin, mmax, 
                                   args = [neq, mmin, pval, qval])[0]

        tmmax = obsmax + delta
        d_t = np.abs(tmmax - mmax)
        mmax = np.copy(tmmax)
        j += 1
        if j > maxiter:
            print 'Kijko-Sellevol-Bayes estimator reached'
            print 'maximum # of iterations'
            d_t = 0.5 * tol
    mmax_sig = np.sqrt(obsmaxsig ** 2. + delta ** 2.)
    return mmax, mmax_sig


def h_smooth(mag):
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

def gauss_cdf_hastings(xval, barx = 0.0, sigx = 1.0):
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
    
def kijko_npg_intfunc_simps(mval, mag, hfact, neq):
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
        p_min = gauss_cdf_hastings((mmin - mag) / hfact)
        p_max = gauss_cdf_hastings((mmax - mag) / hfact)
        p_mag = gauss_cdf_hastings((target_mag - mag) / hfact)
        cdf_func[ival] = ((np.sum(p_mag - p_min)) / 
                          (np.sum(p_max - p_min))) ** neq
    # Now to perform integration via mid-point rule
    intfunc = 0.5 * cdf_func[0] * (mval[1] - mval[0])
    i = 1
    while i < nmval-1:
        intfunc = intfunc + (0.5 * cdf_func[i] * (mval[i + 1] - mval[i - 1]))
        i += 1
    intfunc = intfunc + (0.5 * cdf_func[-1] * (mval[-1] - mval[-2]))

    return intfunc


def kijko_nonparametric_gauss(mag, mag_sig, neq=100,
        nsamp = 51, tol=1.0E-3, maxiter=1E3, obsmax=False):
    '''Function to implement Kijko (2004) Nonparametric Gaussian method
    for estimation of Mmax
    :param mag: Observed magnitudes
    :type mag: numpy.ndarray
    :param mag_sig: Uncertainties on observed magnitudes
    :type mag_sig: numpy.ndarray
    :keyword nsamp: Number of sampling points of integral function
    :type nsamp: Integer
    :keyword tol: Intergral tolerance
    :type tol: Float
    :keyword maxiter: Maximum number of Iterations
    :type maxiter: Int    
    :keyword obsmax: Maximum Observed Magnitude (if not in magnitude array)
                     and its corresponding uncertainty (sigma)
    :type obsmax: Tuple (float) or Boolean
    :return: **mmax** Maximum magnitude and **mmax_sig** corresponding
             uncertainty
    :rtype: Float      
    '''
    
    if not(obsmax):
        # If maxmag is False then maxmag is maximum from magnitude list
        del(obsmax)
        obsmax = np.max(mag)
        obsmaxsig = mag_sig[np.argmax(mag)]

    else:
        obsmaxsig = obsmax[1]
        obsmax = obsmax[0]
    
    #neq = np.shape(mag)[0]
    # Find number_eqs largest events
    if np.shape(mag)[0] <= neq:
        # Catalogue smaller than number of required events
        neq = np.float(np.shape(mag)[0])
    else:
        # Select number_eqs largest events
        mag = np.sort(mag, kind='quicksort')
        mag = mag[-neq:]
    
    mmin = np.min(mag)   
    print mmin, neq, mag
    # Get smoothing factor
    hfact = h_smooth(mag)
    mmax = np.copy(obsmax)
    d_t = 1.E8
    j = 0
    while d_t > tol:
        # Generate exponentially spaced samples
        magval = np.log(np.hstack([np.exp(mmin) + 
            np.arange(0., nsamp - 1., 1.) * ((np.exp(mmax) - np.exp(mmin)) /
            (nsamp - 1.)), np.exp(mmax)]))
        # Evaluate integral function using Simpson's method
        delta = kijko_npg_intfunc_simps(magval, mag, hfact, neq)
        
        tmmax = obsmax + delta
        d_t  = np.abs(tmmax - mmax)
        mmax = np.copy(tmmax)
        j += 1
        if j > maxiter:
            print 'Kijko-Non-Parametric Gaussian estimator reached'
            print 'maximum # of iterations'
            d_t = 0.5 * tol
    mmax_sig = np.sqrt(obsmaxsig ** 2. + delta ** 2.)
    return mmax, mmax_sig

def cumulative_moment(year, mag, iplot = False):
    '''Calculation of Mmax using aCumulative Moment approach, adapeted from
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
    '''Function to calculate mmax wth uncertainty using the cumulative moment
    formulation (adapted from Makropoulos & Burton, 1983) with bootstrapping
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
