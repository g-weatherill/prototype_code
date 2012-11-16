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
which implement declustering algorithms. Implemented
algorithms are:

* GardnerKnopoff
* Afteran
"""

import abc
import numpy as np
#import logging

from catalogue_utilities import (decimal_year, haversine, purge_catalogue)


# Time dist window objects

class Window(object):
    """
    Defines the space and time windows,
    within which an event is identified
    as a cluster.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def calc(self, magnitude):
        """
        Allows to calculate distance and time windows (sw_space, sw_time)
        see reference: `Van Stiphout et al (2010)`.

        :param magnitude: magnitude
        :type magnitude: float
        :returns: distance and time windows
        :rtype: numpy.ndarray
        """
        return


class GardnerKnopoffWindow(Window):
    """
    Gardner Knopoff method for calculating
    distance and time windows
    """

    def calc(self, magnitude):
        """
        >>> import numpy as np
        >>> gw = GardnerKnopoffWindow()
        >>> mag = np.array([4.0, 4.7, 8.9])
        >>> sw_space, sw_time = gw.calc(mag)
        >>> np.allclose(\
                sw_space, np.array([30.07460971, 36.71639218, 121.56820382]))
        True
        >>> np.allclose(\
                sw_time, np.array([0.11332015, 0.27097994, 2.89339106]))
        True
        """

        sw_space = np.power(10.0, 0.1238 * magnitude + 0.983)
        sw_time = np.power(10.0, 0.032 * magnitude + 2.7389) / 365.
        sw_time[magnitude < 6.5] = np.power(
            10.0, 0.5409 * magnitude[magnitude < 6.5] - 0.547) / 365.

        return sw_space, sw_time


class GruenthalWindow(Window):
    """
    Gruenthal method for calculating
    distance and time windows
    """

    def calc(self, magnitude):
        """
        >>> import numpy as np
        >>> gw = GruenthalWindow()
        >>> mag = np.array([4.0, 4.7, 8.9])
        >>> sw_space, sw_time = gw.calc(mag)
        >>> np.allclose(\
                sw_space, np.array([44.65825539, 52.87621383, 120.19384661]))
        True
        >>> np.allclose(\
                sw_time, np.array([0.22553603, 0.45240063, 2.82687845]))
        True
        """

        sw_space = np.exp(1.77 + np.sqrt(0.037 + 1.02 * magnitude))
        sw_time = np.abs(
            (np.exp(-3.95 + np.sqrt(0.62 + 17.32 * magnitude))) / 365.)
        sw_time[magnitude >= 6.5] = np.power(
            10, 2.8 + 0.024 * magnitude[magnitude >= 6.5]) / 365.

        return sw_space, sw_time


class UhrhammerWindow(Window):
    """
    Uhrhammer method for calculating
    distance and time windows
    """

    def calc(self, magnitude):
        """
        >>> import numpy as np
        >>> uw = UhrhammerWindow()
        >>> mag = np.array([4.0, 4.7, 8.9])
        >>> sw_space, sw_time = uw.calc(mag)
        >>> np.allclose(\
                sw_space, np.array([8.95310142, 15.71789701, 460.17184693]))
        True
        >>> np.allclose(\
                sw_time, np.array([0.05747152, 0.0576078, 0.05843231]))
        True
        """

        sw_space = np.exp(-1.024 + 0.804 * magnitude)
        sw_time = np.exp(-2.87 + 1.235 * magnitude / 365.)

        return sw_space, sw_time
#LOGGER = logging.getLogger('mt_logger')

TDW_GARDNERKNOPOFF = 'GardnerKnopoff'
TDW_GRUENTHAL = 'Gruenthal'
TDW_UHRHAMMER = 'Uhrhammer'


time_dist_windows = {TDW_GARDNERKNOPOFF: GardnerKnopoffWindow(),
                     TDW_GRUENTHAL: GruenthalWindow(),
                     TDW_UHRHAMMER: UhrhammerWindow()}


class Declustering(object):
    '''Master class for declustering algorithms'''
    def __init__(self):
        '''Initialise'''
        self.decluster_master = {'GardnerKnopoffType1': GardnerKnopoffType1(),
                                 'Afteran': Afteran()}

    def run_declustering(self, catalogue, config):
        '''Apply declustering algorithm'''
        cluster_index, flag_vector = \
            self.decluster_master[config['algorithm']].decluster(catalogue, 
            config)
        #print cluster_index, flag_vector
        return cluster_index, flag_vector

class CatalogueDecluster(object):
    '''Abstract base class for implementation of declustering algorithms'''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def decluster(self, catalogue, config):
       '''Implements declustering algorithms
       :param catalogue: Catalogue of earthquakes
       :type catalogue: Dictionary
       :param config: Configation parameters
       :type config: Dictionary
       '''
       return


class GardnerKnopoffType1(CatalogueDecluster):
    '''Gardner Knopoff algorithm'''
    
    def decluster(self, catalogue, config):
        """
        :param catalogue: Catalogue of earthquakes
        :type catalogue: Dictionary
        :param config: Configation parameters
        :type config: Dictionary

        :returns: **vcl vector** indicating cluster number, 
              **vmain_shock catalog** containing non-clustered events, 
              **flagvector** indicating which eq events belong to a cluster
        :rtype: numpy.ndarray
        """

        # Get relevent parameters
        neq = len(catalogue['magnitude'])  # Number of earthquakes
        # Get decimal year (needed for time windows)
        year_dec = decimal_year(
            catalogue['year'], catalogue['month'], catalogue['day'])
        # Get space and time windows corresponding to each event
        sw_space, sw_time = \
            time_dist_windows[config['window_opt']].\
            calc(catalogue['magnitude'])
        eqid = np.arange(0, neq, 1)  # Initial Position Identifier

        # Pre-allocate cluster index vectors
        vcl = np.zeros(neq, dtype=int)
        #print catalogue['magnitude'], year_dec
        # Sort magnitudes into descending order
        id0 = np.flipud(np.argsort(catalogue['magnitude'], kind='heapsort'))
        #m = m[id0]
        #catalog_matrix = catalog_matrix[id0, :]
        mag = catalogue['magnitude'][id0]
        longitude = catalogue['longitude'][id0]
        latitude = catalogue['latitude'][id0]
        sw_space = sw_space[id0]
        sw_time = sw_time[id0]
        year_dec = year_dec[id0]
        eqid = eqid[id0]
        flagvector = np.zeros(neq, dtype=int)
        #Begin cluster identification
        clust_index = 0
        for i in range(0, neq - 1):
            if vcl[i] == 0:
                #print vcl[i], mag, neq
                # Find Events inside both fore- and aftershock time windows
                dt = year_dec - year_dec[i]
                vsel = np.logical_and(
                    dt >= (-sw_time[i] * config['fs_time_prop']),
                    dt <= sw_time[i], 
                    flagvector == 0)
                # Of those events inside time window, find those inside distance
                # window
                vsel1 = haversine(longitude, latitude, longitude[i], 
                                  latitude[i]) <= sw_space[i]
                vsel[vsel] = vsel1
                temp_vsel = np.copy(vsel)
                temp_vsel[i] = False
                if any(temp_vsel):
                    # Allocate a cluster number
                    vcl[vsel] = clust_index + 1
                    flagvector[vsel] = 1
                    # For those events in the cluster before the main event,
                    # flagvector is equal to -1
                    temp_vsel[dt >= 0.0] = False
                    flagvector[temp_vsel] = -1
                    flagvector[i] = 0
                    clust_index += 1

        # Re-sort the catalog_matrix into original order
        id1 = np.argsort(eqid, kind='heapsort')
        eqid = eqid[id1]
        #catalog_matrix = catalog_matrix[id1, :]
        vcl = vcl[id1]
        flagvector = flagvector[id1]
        return vcl, flagvector
        #if config['purge']:
        #    # Now to produce a catalogue with aftershocks purged
        #    vmain_shock = purge_catalogue(catalogue, flagvector)
        #else:
        #    vmain_shock = None
        #
        #return vcl, vmain_shock, flagvector


class Afteran(CatalogueDecluster):
    '''AFTERAN declustering algorithm.
    ||(Musson, 1999, "Probabilistic Seismic Hazard Maps for the North Balkan
       region", Annali di Geofisica, 42(6), 1109 - 1124) ||'''

    def decluster(self, catalogue, config):
        '''catalogue_matrix, window_opt=TDW_GARDNERKNOPOFF, time_window=60.):

        :param catalog_matrix: eq catalog in a matrix format with these columns in
                                order: `year`, `month`, `day`, `longitude`,
                                `latitude`, `Mw`
        :type catalog_matrix: numpy.ndarray
        :keyword window_opt: method used in calculating distance and time windows
        :type window_opt: string
        :keyword time_window: Length (in days) of moving time window
        :type time_window: positive float
        :returns: **vcl vector** indicating cluster number, **vmain_shock catalog**
                  containing non-clustered events, **flagvector** indicating
                  which eq events belong to a cluster
        :rtype: numpy.ndarray
        '''

        #Convert time window from days to decimal years
        time_window = config['time_window'] / 365.

        # Pre-processing steps are the same as for Gardner & Knopoff
        # Get relevent parameters
        mag = catalogue['magnitude']
        neq = np.shape(magnitude)[0]  # Number of earthquakes
        # Get decimal year (needed for time windows)
        year_dec = decimal_year(catalogue['year'], catalogue['month'],
                                catalogue['day'])

        # Get space windows corresponding to each event
        sw_space = time_dist_windows[config['window_opt']].calc(mag)[0]

        eqid = np.arange(0, neq, 1)  # Initial Position Identifier

        # Pre-allocate cluster index vectors
        vcl = np.zeros((neq, 1), dtype=int)
        flagvector = np.zeros((neq, 1), dtype=int)
        # Sort magnitudes into descending order
        id0 = np.flipud(np.argsort(mag, kind='heapsort'))
        mag = mag[id0]
        #catalogue_matrix = catalogue_matrix[id0, :]
        longitude = catalogue['longitude'][id0]
        latitude = catalogue['latitude'][id0]
        sw_space = sw_space[id0]
        year_dec = year_dec[id0]
        eqid = eqid[id0]

        i = 0
        clust_index = 0
        while i < neq:
            if vcl[i] == 0:
                # Earthquake not allocated to cluster - perform calculation
                # Perform distance calculation
                mdist = haversine(longitude, latitude,
                                  longitude[i], latitude[i])

                # Select earthquakes inside distance window and not in cluster
                vsel = np.logical_and(mdist <= sw_space[i], vcl == 0).flatten()
                dtime = year_dec[vsel] - year_dec[i]

                nval = np.shape(dtime)[0] #Number of events inside valid window
                # Pre-allocate boolean array (all True)

                vsel1 = self._find_aftershocks(dtime, nval, time_window)
                vsel2 = self._find_foreshocks(dtime, nval, time_window, vsel1)

                temp_vsel = np.copy(vsel)
                temp_vsel[vsel] = np.logical_or(vsel1, vsel2)
                if np.shape(np.nonzero(temp_vsel)[0])[0] > 1:
                    # Contains clustered events - allocate a cluster index
                    vcl[temp_vsel] = clust_index + 1
                    # Remove mainshock from cluster
                    vsel1[0] = False
                    # Assign markers to aftershocks and foreshocks
                    temp_vsel = np.copy(vsel)
                    temp_vsel[vsel] = vsel1
                    flagvector[temp_vsel] = 1
                    vsel[vsel] = vsel2
                    flagvector[vsel] = -1
                    clust_index += 1
            i += 1

        # Now have events - re-sort array back into chronological order
        # Re-sort the data into original order
        id1 = np.argsort(eqid, kind='heapsort')
        eqid = eqid[id1]
        #catalogue_matrix = catalogue_matrix[id1, :]
        vcl = vcl[id1]
        flagvector = flagvector[id1]
        return vcl.flatten(), flagvector.flatten()

        #if config['purge']:
        #    # Now to produce a catalogue with aftershocks purged
        #    vmain_shock = purge_catalogue(catalogue, flagvector)
        #else:
        #    vmain_shock = None
        #return vcl.flatten(), vmain_shock, flagvector.flatten()

    def _find_aftershocks(self, dtime, nval, time_window):
        """
        Searches for aftershocks within the moving
        time window
        :param dtime: time since main event
        :type dtime: numpy.ndarray
        :param nval: number of events in search window
        :type nval: int
        :param time_window: Length (in days) of moving time window
        :type time_window: positive float
        :returns: **vsel** index vector for aftershocks
        :rtype: numpy.ndarray
        """

        vsel = np.array(np.ones(nval), dtype=bool)
        initval = dtime[0]  # Start with the mainshock

        j = 1
        while j < nval:
            ddt = dtime[j] - initval
            # Is event after previous event and within time window?
            vsel[j] = np.logical_and(ddt >= 0.0, ddt <= time_window)
            if vsel[j]:
                # Reset time window to new event time
                initval = dtime[j]
            j += 1
        return vsel


    def _find_foreshocks(self, dtime, nval, time_window, vsel_aftershocks):
        """
        Searches for foreshocks within the moving
        time window
        :param dtime: time since main event
        :type dtime: numpy.ndarray
        :param nval: number of events in search window
        :type nval: int
        :param time_window: Length (in days) of moving time window
        :type time_window: positive float
        :param vsel_aftershocks: index vector for aftershocks
        :type vsel_aftershocks: numpy.ndarray
        :returns: **vsel** index vector for foreshocks
        :rtype: numpy.ndarray
        """

        j = 1
        vsel = np.array(np.zeros(nval), dtype=bool)
        initval = dtime[0]

        while j < nval:
            if vsel_aftershocks[j]:
            # Event already allocated as an aftershock - skip
                j += 1
            else:
                ddt = dtime[j] - initval
                # Is event before previous event and within time window?
                vsel[j] = np.logical_and(ddt <= 0.0,
                                          ddt >= -(time_window))
                if vsel[j]:
                # Yes, reset time window to new event
                    initval = dtime[j]
            j += 1

        return vsel
