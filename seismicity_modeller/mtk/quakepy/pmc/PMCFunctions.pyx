# -*- coding: utf-8 -*-
#
# quakepy/pmc/PMCFunctions.pyx
# $Id: PMCFunctions.pyx 335 2012-06-08 12:08:36Z tyrone $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2007-2009 by Fabian Euchner and Danijel Schorlemmer     #
#    fabian@fabian-euchner.de                                              #
#                                                                          #
#    This program is free software; you can redistribute it and#or modify  #
#    it under the terms of the GNU General Public License as published by  #
#    the Free Software Foundation; either version 2 of the License, or     #
#    (at your option) any later version.                                   #
#                                                                          #
#    This program is distributed in the hope that it will be useful,       #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#    GNU General Public License for more details.                          #
#                                                                          #
#    You should have received a copy of the GNU General Public License     #
#    along with this program; if not, write to the                         #
#    Free Software Foundation, Inc.,                                       #
#    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
############################################################################

"""
The QuakePy package
http://www.quakepy.org
"""

__version__  = '$Id: PMCFunctions.pyx 335 2012-06-08 12:08:36Z tyrone $'
__revision__ = '$Revision: 335 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import numpy as np

## Cython stuff

cimport  numpy as np
ctypedef np.float_t DTYPE_t

def computeCombinedProbabilities( np.ndarray probabilities not None, 
    int stationsRequiredForTriggering ):
    """
    compute combined probability for a vector of station probabilities
    using convolution.
    Idea and reference implementation by Matthias Holschneider
    (May 2012, Universität Potsdam)

    input parameter probabilities has to be numpy array of propabilities
    of involved stations (not sorted)

    replaces method in PMC class
    """

    if probabilities.shape[0] < stationsRequiredForTriggering:
        return 0.0

    # cdefs of local helpers:
    cdef DTYPE_t u
    cdef DTYPE_t v
    cdef int i
    cdef int j

    cdef np.ndarray res = np.zeros( (stationsRequiredForTriggering + 1), dtype=np.float )

    # start:
    res[0] = probabilities[0,1]
    res[1] = probabilities[0,0]

    for i in xrange(1, len(probabilities[:,0])):
        u = res[0]
        res[0] = res[0] * probabilities[i,1]
        for j in xrange(1, (stationsRequiredForTriggering)):
            v = probabilities[i,0] * u + probabilities[i,1] * res[j]
            u = res[j]
            res[j] = v

    return 1.0 - res.sum()


def calcRawProbability( np.ndarray pickInfo not None, float magnitude, 
    float distance, calcMagnitudeFromDistance ):
    """
    calc raw probability for given magnitude and distance

    was intended to replace method in PMCInventory class, but not 
    used because it is SLOWER
    """

    cdef int     numSample    = 0

    cdef DTYPE_t probability  = 0.0
    cdef DTYPE_t logDistance  = 0.0

    # if pickInfo has no rows, return probability 0.0
    if pickInfo.shape[0] == 0:
        return ( probability, numSample )

    # create numpy arrays with 4 columns and rows as in pickInfo
    cdef np.ndarray values = np.zeros( ( pickInfo.shape[0], 4 ), 
        dtype=np.float )
    cdef np.ndarray sample = np.zeros( ( pickInfo.shape[0], 4 ), 
        dtype=np.float )
    cdef np.ndarray nonSample = np.zeros( ( pickInfo.shape[0], 4 ), 
        dtype=np.float )
    cdef np.ndarray newSample = np.zeros( ( pickInfo.shape[0], 4 ), 
        dtype=np.float )
    cdef np.ndarray selection = np.zeros( ( pickInfo.shape[0], 4 ), 
        dtype=np.float )
    cdef np.ndarray sortedSelection = np.zeros( ( pickInfo.shape[0], 4 ), 
        dtype=np.float )

    cdef np.ndarray sel = np.zeros( pickInfo.shape[0], dtype=np.int )
    cdef np.ndarray indices = np.zeros( pickInfo.shape[0], dtype=np.int )

    # 1st column: picked or not picked?
    values[:,0] = pickInfo[:,2]

    # 2nd column: magnitude differences
    values[:,1] = pickInfo[:,1] - magnitude

    # 3rd column: translate difference in distances into a magnitude difference
    logDistance = calcMagnitudeFromDistance( distance )
    values[:,2] = pickInfo[:,4] - logDistance

    # 4rd column: calculate total difference in "magnitude units"
    values[:,3] = np.sqrt( values[:,1]**2 + values[:,2]**2 )

    # select only picks with a "magnitude difference" of less than 0.1
    sel = ( values[:,3] < 0.1 )             # selector: TODO use magnitude threshold

    # check if selector contains at least one 'True'
    if sel.sum() != 0.0:
        sample    = values[sel.T,:]
        numSample = sample.shape[0]
    else:
        numSample   = 0

    if ( numSample < 10 ):                                    # not enough samples available
        nonSample = values[np.logical_not(sel),:]             # invert selector

        sel       = np.logical_and( ( nonSample[:,1] <= 0 ),
                                    ( nonSample[:,2] >= 0 ) ) # select points w/ larger distance, smaller magnitude

        selection = nonSample[sel.T,:]
        indices   = selection[:, 3].argsort( axis=0 )          # sort for distance

        if indices.shape[0] < (10 - numSample):               # enough points to fill array up to 10?

            probability = 0.0                                 # no, give up (return zero)
        else:
            sortedSelection = selection[indices]

            if numSample > 0:
                newSample = np.concatenate( 
                    ( sample, sortedSelection[0:(10-numSample),:] ), axis=0 )
            else:
                newSample = sortedSelection[0:10, :]

            probability = newSample[:, 0].sum() / 10.0
            numSample = 10
    else:
        probability = sample[:, 0].sum() / numSample

    return ( probability, numSample )