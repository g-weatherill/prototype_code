# -*- coding: utf-8 -*-
#
# quakepy/pmc/PMC_NZ.py
# $Id: PMC_NZ.py 336 2012-07-09 14:10:49Z tyrone $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2007-2011 by Fabian Euchner and Danijel Schorlemmer     #
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

__version__  = '$Id: PMC_NZ.py 336 2012-07-09 14:10:49Z tyrone $'
__revision__ = '$Revision: 336 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import copy
import re
import cPickle
import math
import numpy

import pyRXP
from   xml.sax import saxutils

from mx.DateTime import DateTime, TimeDelta
from mx.DateTime.ISO import ParseDateTimeUTC

sys.path.append('..')

from QPCore import *
from PMC import *
from PMCInventory import *
from PMCData import *
from PMCMetadata import *
from QPCatalog import *

NETWORK_CODE = 'NZ'
SUBNETWORKS = {
        'AKN': { 'stations': ('ABAZ', 'ETAZ', 'HBAZ', 'KAAZ', 'KBAZ', 'KUZ', 'MKAZ', 'RVAZ', 'WTAZ'),
                 'trigger_condition': 5,
                 'description': 'Auckland network'
               },
        'CHN': { 'stations': ('BHHZ', 'BKZ', 'CNZ', 'COVZ', 'DRZ', 'FWVZ', 'KATZ', 'KAVZ', 'KRVZ', 
                              'MGZ', 'MOVZ', 'MTVZ', 'NGZ', 'OTVZ', 'PKVZ', 'RITZ', 'TRVZ', 'TUVZ', 
                              'TWVZ', 'WHVZ', 'WNVZ', 'WPVZ', 'WTVZ'),
                 'trigger_condition': 7,
                 'description': 'Tongariro National Park network'
               },
        'CYN': { 'stations': ('CFC', 'CMCZ', 'LRCZ', 'LSCZ', 'MHZ', 'MMCZ', 'MSCZ', 'SBCZ', 'TBC', 'TLC'),
                 'trigger_condition': 5,
                 'description': 'Clyde network'
               },
        'HBO': { 'stations': ('HNH', 'MAHZ', 'MOH', 'MRH', 'PAHZ', 'TAHZ', 'TEHZ', 'TTH', 'WAHZ', 'WHH'),
                 'trigger_condition': 3,
                 'description': 'Hawke\'s Bay network (old)'
               },
        'HBN': { 'stations': ('BHHZ', 'BKZ', 'CKHZ', 'DVHZ', 'KAHZ', 'KNZ', 'KRHZ', 'MCHZ', 'NMHZ', 
                              'PNHZ', 'PRHZ', 'PXZ', 'RAHZ', 'TSZ', 'WHHZ', 'WPHZ'),
                 'trigger_condition': 6,
                 'description': 'Hawke\'s Bay network (new)' 
               },
        'RTN': { 'stations': ('EDRZ', 'HARZ', 'HLRZ', 'HRRZ', 'KARZ', 'LIRZ', 'MARZ', 'MKRZ', 'OMRZ', 
                              'OPRZ', 'RRRZ', 'TARZ', 'TAZ', 'TGRZ', 'URZ', 'UTU'),
                 'trigger_condition': 6,
                 'description': 'Rotorua/Okataina network' 
               },
        'TPN': { 'stations': ('ALRZ', 'ARAZ', 'BKZ', 'ECNZ', 'HATZ', 'HITZ', 'HUTZ', 'KATZ', 'KETZ', 'KAVZ', 
                              'KRVZ', 'MGZ', 'POIZ', 'RATZ', 'RITZ', 'TUTZ', 'WATZ', 'WHTZ'),
                 'trigger_condition': 5,
                 'description': 'Taupo network' 
               },
        'TRN': { 'stations': ('DFE', 'KHEZ', 'LREZ', 'MHEZ', 'NMEZ', 'PKE', 'RAEZ', 'TKEZ', 'VRZ'),
                 'trigger_condition': 4,
                 'description': 'Taranaki network' 
               },
        'WLN': { 'stations': ('AMW', 'ANWZ', 'BBW', 'BFZ', 'BHW', 'BLW', 'BSWZ', 'CAW', 'CCW', 'CMWZ', 'DIW', 
                              'DUWZ', 'GFW', 'HOWZ', 'KIW', 'MOW', 'MRW', 'MRZ', 'MSWZ', 'MTW', 'NNZ', 'OGWZ', 
                              'OTW', 'PAWZ', 'PLWZ', 'POWZ', 'PRWZ', 'QHW', 'SNZO', 'TCW', 'TIWZ', 'TMWZ', 
                              'TRWZ', 'TUWZ', 'WDW', 'WEL', 'WHW'),
                 'trigger_condition': 7,
                 'description': 'Wellington network' 
               },
        'NWN': { 'stations': ('OUZ', 'WCZ', 'KUZ', 'WLZ', 'TOZ', 'MOZ', 'HIZ'),
                 'trigger_condition': 3,
                 'description': 'Northwest North Island network' 
               },
        'NEN': { 'stations': ('KUZ', 'HBZ', 'MXZ', 'WLZ', 'TOZ', 'PUZ', 'URZ', 'MWZ', 'NOZ', 'KNZ', 'BKZ'),
                 'trigger_condition': 3,
                 'description': 'Northeast North Island network' 
               },
        'WNI': { 'stations': ('MOZ', 'HIZ', 'OIZ', 'VRZ', 'NRZ', 'BSZ', 'WAZ', 'TSZ'),
                 'trigger_condition': 3,
                 'description': 'Western North Island network' 
               },
        'CNZ': { 'stations': ('WEL', 'SNZO', 'MNG', 'KHZ', 'PGZ', 'NNZ', 'THZ', 'BSZ', 
                              'WAZ', 'BFZ', 'MRZ', 'PWZ', 'QRZ', 'TSZ', 'DSZ', 'BKZ'),
                 'trigger_condition': 3,
                 'description': 'Central New Zealand network' 
               },
        'CNI': { 'stations': ('THZ', 'KHZ', 'MQZ', 'CRLZ', 'LTZ', 'DSZ', 'WVZ', 'EWZ', 
                              'RPZ', 'LMZ', 'FOZ', 'BWZ', 'LBZ'),
                 'trigger_condition': 3,
                 'description': 'Central South Island network' 
               },
        'SSI': { 'stations': ('EWZ', 'RPZ', 'MSZ', 'ODZ', 'AXZ', 'EAZ', 'MLZ', 'DCZ', 'WHZ', 
                              'TUZ', 'BWZ', 'LBZ', 'WKZ', 'LMZ', 'FOZ', 'JCZ', 'SIZ', 'HHSZ'),
                 'trigger_condition': 3,
                 'description': 'Southern South Island network' 
               } }

# Stations not to consider, according to Kevin Fenaughty's recommendation
IGNORED_STATIONS = ('OH1', 'OH2', 'OH3', 'OH4', 'WK1', 'WK2', 'WK3', 'WK4', 'COVZ')

# Stations to unselect 
# (added to source code by Danijel, origin of recommendation is unclear)
STATIONS_TO_UNSELECT = ( 'API', 'BUN', 'DND', 'HLL', 'HNZ', 'NDF', 'NPR', 
    'NUE', 'SUV', 'TAK', 'VND', 'WELA', 'WELO' )

class PMC_NZ( PMC ):
    """
    QuakePy: PMC_NZ
    PMC method for Geonet network (New Zealand)
    """

    # these attributes are redefined for Geonet network

    # this is the triggering condition until 1987-01-01
    # from 1987-01-01, the triggering condition is set per subnetwork
    stationsRequiredForTriggering = 3

    probThreshold = 0.999

    # minimum number of arrivals an ontime must contain
    # time intervals containing fewer picks will not be counted
    ontimeMinPickCtr = 20         

    useDistStyle = 3

    # time shift in NZ time zone (hours)
    timeZoneShiftStandard = 12.0
    timeZoneShiftDaylight = 13.0
    timeZoneShift         = timeZoneShiftStandard
    
    networkCode = NETWORK_CODE
    stationOnTimeStartLowerLimit = QPDateTime( [1800, 1, 1, 0, 0, 0] )
    stationOnTimeEndUpperLimit   = QPDateTime( [2100, 1, 1, 0, 0, 0] )

    subnetworks = None
    ignoredStations = IGNORED_STATIONS

    noOfHeaderlinesToSkipInStationfile = 5

    def __init__( self, stationlist_filename=None, dist_style=None, 
        useSubnetworks=False, **kwargs ):
        super( PMC_NZ, self ).__init__( stationlist_filename, 
            self.useDistStyle, **kwargs )

        if useSubnetworks is True:
            self.subnetworks = SUBNETWORKS

        if 'subnet' in kwargs and isinstance( kwargs['subnet'], basestring ):
            self.setSubNetwork( kwargs['subnet'] )

    def _calcMagnitudeFromDistance( self, distance ):

        # Robinson formula
        return math.log10(distance/111.2) + (0.0029 * distance/111.2)


    def _stationToBeUsed( self, station ):
        """Exclude stations that are explicitly in a list to be ignored,
        or stations with 4-character names that end with S, B, or C.
        """

        if (    ( station.stationCode in self.ignoredStations ) 
             or (     len(station.stationCode) == 4 
                  and station.stationCode.endswith(('S', 'B', 'C')) ) ):
            print "Station excluded: %s" % station.stationCode
            return False

        # Select only stations from sub-network, if subnetwork is set
        if (     self.subnetwork is not None 
             and station.stationCode not in SUBNETWORKS[subnetwork]['stations'] ):
            print "Station %s excluded, since it's not in sub-network %s" % (
                station.stationCode, self.subnetwork )
            return False

        return True


    def importStations( self, filename, encoding='ascii' ):
        """
        import station list of Geonet

        example data (note: ruler shown below is not part of data)

        BEGIN IMS1.0
        MSG_TYPE DATA
        MSG_ID 02-DEC-2009_15:09:06 NZL_WEL
        DATA_TYPE STATION IMS1.0
        Net       Sta   Type  Latitude  Longitude Coord Sys     Elev   On Date   Off Date
        NZ        001A       -35.72690  174.31920 NZGD49             1970/04/18 1980/06/13
        NZ        002A       -41.73310  171.80000 NZGD49             1968/05/24 1968/07/23
        ...

                  1         2         3         4         5         6         7         8
        0123456789012345678901234567890123456789012345678901234567890123456789012345678901
        ----------------------------------------------------------------------------------
        XX        YUPC       -38.31750  175.55029 NZGD49       0.540 2001/01/15 2001/06/27
        XX        ZIHA       -43.88700  169.04970 NZGD49             1999/08/01 2005/05/21

        STOP
        """

        ftLines = getQPDataSource( filename )

        lines_processed = 0
        nLines          = 0

        for sLine in ftLines: 

            nLines += 1
            
            # skip first 5 header lines
            if (     (nLines > self.noOfHeaderlinesToSkipInStationfile) 
                 and (len(sLine.strip()) != 0) 
                 and (sLine.strip() != 'STOP') 
                 and (not sLine[10:13].strip().isdigit()) 
                 and (sLine[0:2] == 'NZ') ):

                sNetworkCode = 'NZ'
                sStationCode = sLine[10:15].strip()
                fLatitude    = float(sLine[21:30])
                fLongitude   = float(sLine[31:41])
                sCoordRef    = sLine[42:48].strip()

                try:
                    fElevation   = float(sLine[54:60]) * 1000 # in meters
                except:
                    fElevation = 0

                nOnYear      = int(sLine[61:65])
                nOnMonth     = int(sLine[66:68])
                nOnDay       = int(sLine[69:71])
                dtOnTime = QPDateTime([nOnYear, nOnMonth, nOnDay, 0, 0, 0])
                sOffYear     = sLine[72:76].strip()
                if len(sOffYear) == 0:
                    dtOffTime = self.stationOnTimeEndUpperLimit
                else:
                    nOffMonth = int(sLine[77:79])
                    nOffDay   = int(sLine[80:82])
                    dtOffTime = QPDateTime([int(sOffYear), nOffMonth, nOffDay, 0, 0, 0])

                # look in station list if network/station combination is already there
                sta_idx = self._getStationIdx(sNetworkCode, sStationCode)

                self.stations[sta_idx].networkCode  = sNetworkCode
                self.stations[sta_idx].stationCode  = sStationCode

                self.stations[sta_idx].latitude     = fLatitude
                self.stations[sta_idx].longitude    = fLongitude
                self.stations[sta_idx].elevation    = fElevation

                # new station, create new channel
                curr_channel = PMCChannel()
                self.stations[sta_idx].channels.append( curr_channel )

                # set channel data
                curr_channel.waveformID = WaveformStreamID(sNetworkCode, sStationCode)
                curr_channel.onTime.append( [ dtOnTime, dtOffTime] )

                lines_processed += 1

        # set station-wide station.onTime
        for curr_sta in self.stations:

            # loop over channels, add onTime of all channels to station onTime
            # Note: could be improved by checking for overlapping time intervals
            for curr_ch in curr_sta.channels:
                curr_sta.onTime.extend( curr_ch.onTime )

        print " station import complete, read %s station lines, length of station array %s" % (
            lines_processed, len( self.stations ) )

        # check if read numbers match
        if len(self.stations) != lines_processed:
            print " ATTENTION: numbers of processed stations do not match!"


    def getCatalogFiles( self, metadata ):
        """Add list of catalog files for given start/end date to catalog_file
        file list.

        Naming convention for NZ catalog files (gnuzipped), one file per month:
        nz.YYYY-MM.xml.gz
        """

        metadata.getCatalogFilesAll()

        # startDate and endDate have to be mxDateTime.DateTimeType
        if (     ( metadata.startDate is not None )
             and ( metadata.endDate is not None )
             and isinstance( metadata.startDate, DateTimeType )
             and isinstance( metadata.endDate, DateTimeType ) ):

            for time_year in xrange( metadata.startDate.year, 
                metadata.endDate.year + 1 ):

                # get catfiles for time period
                for curr_catfile in metadata.catalogFilesAll:

                    datecode = re.match( r'^.*nz\.(\d{4})-(\d{2})\.xml\.gz$', 
                        curr_catfile )

                    year = int( str( datecode.group(1) ) )
                    month = int( str( datecode.group(2) ) )

                    if year == time_year:
                        if not (    (     year == metadata.startDate.year 
                                      and month < metadata.startDate.month )
                                 or (     year == metadata.endDate.year 
                                      and month > metadata.endDate.month ) ):

                            metadata.catalogFiles.append( curr_catfile )
                            print " PMC_NZ: added catalog file %s to current"\
                                  " catalog" % curr_catfile

    def getCatalog( self, catalog_file ):

        qpc = QPCatalog( idstyle='numeric' )
        qpc.readXML( catalog_file, compression='gz', minimumDataset=True )
        return qpc

# ------------------------------------------------------------------------------

    def preprocess( self, catalog=None, timePeriod=None ):
        # this function is redefined for Geonet network
        
        # - preprocess inventory
        # (1) exclude some networks
        # (2) exclude all channels which have no 'Z' component
        # (3) exclude all stations not operating in given time period
        #     timePeriod is [ startDate, endDate ]
        
        # - preprocess catalog
        # (4) use only Geonet locations with given parameters (time, lat, lon, magnitude)
        #     depth: if not given, use standard value
        # (5) exclude events that have no picks associated
        # (6) exclude all picks that are not 'P' picks: use only those starting with 'P'
        
        self.preprocessInventory( timePeriod )

        # preprocess catalog only if deck files are read into QPCatalog object
        if catalog is not None:
            self.preprocessCatalog( catalog, timePeriod )

            
    def preprocessInventory( self, timePeriod=None ):
        # this function is redefined for Geonet network
        
        # (1) exclude some networks
        # (2) exclude all channels which have no 'Z' component
        # (3) exclude all stations not operating in given time period
        #     timePeriod is [ startDate, endDate ]

        if ( not self.stations ):
            raise ValueError, "PMC_NZ.preprocessInventory - no station list loaded"
        
        unused_networks = ()

        # over stations
        for curr_sta_idx in reversed( xrange( len( self.stations ) ) ):
            
            # in excluded network?
            if self.stations[curr_sta_idx].networkCode in unused_networks:
                del self.stations[curr_sta_idx]
            else:
                # over channels
                for curr_cha_idx in reversed( xrange( len( self.stations[curr_sta_idx].channels ) ) ):
                    
                    # has channel 'Z' component? / does channel operate in given time period?
                    if not self._channelOperatingDuringPeriod( 
                        self.stations[curr_sta_idx].channels[curr_cha_idx],
                        timePeriod ):
                        del self.stations[curr_sta_idx].channels[curr_cha_idx]
                        
                if len( self.stations[curr_sta_idx].channels ) == 0:
                    del self.stations[curr_sta_idx]
        
        # unselect stations 
        for unselect_sta in STATIONS_TO_UNSELECT:
            self.unselectStation( 'NZ', unselect_sta )

        return True

            
    def preprocessCatalog( self, catalog, timePeriod=None ):
        # this function is redefined for Geonet network

        # (1) if timePeriod is set, exclude events outside time period
        # (2) use only Geonet locations with given parameters (time, lat, lon, magnitude)
        #     depth: if not given, use standard value
        # (3) exclude events that have no picks associated
        # (4) exclude all picks that are not 'p' picks: use only those starting with 'p'

        if catalog is None:
            raise ValueError, "PMC_NZ.preprocessCatalog - no valid catalog given"

        if timePeriod is not None:

            if not (      isinstance( timePeriod[0], DateTimeType )
                     and  isinstance( timePeriod[1], DateTimeType ) ):

                error_str = "PMC_NZ.preprocessCatalog - timePeriod has incorrect type"
                raise TypeError, error_str
            else:
                startTime = QPDateTime( timePeriod[0] )
                endTime   = QPDateTime( timePeriod[1] )
                
        event_cnt = len( catalog.eventParameters.event )
        for curr_ev_idx in reversed( xrange( event_cnt ) ):
            
            event = catalog.eventParameters.event[curr_ev_idx]

            #print " analyzing event no. %s w/ ID %s" % ( curr_ev_idx, event.publicID )
            
            # assume that preferred origin is Geonet origin
            origin = event.getPreferredOrigin()

            # focal time in catalog and start/end times are UTC
            if timePeriod is not None:

                if (    ( origin.time.value < startTime )
                     or ( origin.time.value > endTime ) ):

                    del catalog.eventParameters.event[curr_ev_idx]
                    del event
                    continue

            # check whether there are magnitudes
            # if positive, get first magnitude, assume that this 
            # is the magnitude that rresponds to preferred origin
            ori_mag = event.getMagnitudes( origin )
            if len( ori_mag ) > 0:
                magnitude = ori_mag[0]
                if (not hasattr(magnitude, 'mag')) or (magnitude.mag is None):
                    del catalog.eventParameters.event[curr_ev_idx]

                    #print " skipped event w/o complete magnitude %s, event count is now %s" % ( event.publicID, catalog.size() )
                    del event
                    continue
            else:
                # no magnitude, skip event
                del catalog.eventParameters.event[curr_ev_idx]
                
                #print " skipped event w/o magnitude %s, event count is now %s" % ( event.publicID, catalog.size() )
                del event
                continue

            #rint " event %s, magnitude %s, origin %s no. of picks %s" % ( event.publicID, magnitude.publicID, origin.publicID, len(event.pick) )

            # check if origin is complete
            # - lat, lon, depth
            # focal time is always there and we have already checked for magnitude
            # Delete events with depths > 40 (CSEP testing depths)
            if (    ( ( not hasattr( origin, 'latitude' ) ) or ( origin.latitude is None ) )
                 or ( ( not hasattr( origin, 'longitude' ) ) or ( origin.longitude is None ) )
                 or ( ( not hasattr( origin, 'depth' ) ) or ( origin.depth is None ) )
                 or ( origin.depth.value > 40.0 ) ):

                del catalog.eventParameters.event[curr_ev_idx]

                #print " skipped event w/o coordinates or proper depth %s, event count is now %s" % ( event.publicID, catalog.size() )
                del event
                continue

            # check if there are any picks
            if ( ( event.pick is None ) or ( len(event.pick) == 0 ) ):

                del catalog.eventParameters.event[curr_ev_idx]

                #print " skipped event w/o picks %s, event count is now %s" % ( event.publicID, catalog.size() )
                del event
                continue
            
            elif ( ( event.pick is not None ) and ( len(event.pick) > 0 ) ):

                # check for 'p*' picks
                for curr_pick_idx in reversed( xrange( len(event.pick) ) ):

                    curr_pick         = event.pick[curr_pick_idx]
                    curr_arrivals_idx = origin.getArrivalsIdx( curr_pick )

                    # check if pick has an associated arrival
                    if len( curr_arrivals_idx ) == 0:
                        
                        del event.pick[curr_pick_idx]
                        
                        #print " deleted pick no. %s, id %s since no arrival is associated" %  ( curr_pick_idx, curr_pick.publicID )

                        del curr_pick

                    else:

                        curr_arr_idx   = curr_arrivals_idx[0]
                        curr_arr_phase = origin.arrival[curr_arr_idx].phase.code
                        
                        if not (curr_arr_phase.lower().startswith('p')):
                            
                            del origin.arrival[curr_arr_idx]
                            del event.pick[curr_pick_idx]

                            #print " deleted %s phase for pick no. %s, id %s" %  ( curr_arr_phase, curr_pick_idx, curr_pick.publicID )

                            del curr_pick
                    
                # loop over picks has finished
                # if event has no more picks, delete event
                if len(event.pick) == 0:

                    del catalog.eventParameters.event[curr_ev_idx]
                    
                    #print " deleted event with no remaining picks %s, event count is now %s" % ( event.publicID, catalog.size() )
                    del event
                    continue

                # print " event was not deleted, proceed to next"
                
        return True


    def initAnnotation( self ):
        """
        has been redefined for class PMC_NZ derived from PMC
        set annotation for QPGrid output
        """

        self.annotation.setProperty(
            title = 'Completeness of the Geonet network',
            creator = ( 'Danijel Schorlemmer, Fabian Euchner, and Matthew C. Gerstenberger', ),
            publisher = ( 'Danijel Schorlemmer, Fabian Euchner, and Matthew C. Gerstenberger', ),
            rights = ( 'Copyright by Danijel Schorlemmer, Fabian Euchner, and Matthew C. Gerstenberger',
                       'Licensed under Creative Commons Attribution-Noncommercial-Share Alike 3.0',
                       'http://creativecommons.org/licenses/by-nc-sa/3.0/' ),
            subject = ( 'seismology', 'earthquake', 'detection probability', 'PMC', 'seismic network', 'completeness',
                        'magnitude of completeness' ),
            source = ( 'http://completenessweb.org', 'http://quakepy.org' ),
            bibliographicCitation = ( "D. Schorlemmer and J. Woessner: Probability of Detecting an Earthquake, Bull. Seismol. Soc. Am., 98(5), 2103-2117, 2008, doi:10.1785/0120070105.", ),
            comment = ( "Please cite the given reference when using the data.", )
        )

        # not used at the moment
        # identifier  version  acknowledgment

        # for these we need run-time information
        # date  coverageTemporal  coverageSpatial
        
        return True
    
## ---------------------------------------------------------------------------

