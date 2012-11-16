# -*- coding: utf-8 -*-
#
# quakepy/pmc/PMC_JMA.py
# $Id: PMC_JMA.py 321 2011-10-16 06:24:10Z tyrone $
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

__version__  = '$Id: PMC_JMA.py 321 2011-10-16 06:24:10Z tyrone $'
__revision__ = '$Revision: 321 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import copy
import cPickle
import math
import numpy

import pyRXP
from   xml.sax import saxutils

import matplotlib
matplotlib.use('PS')
from pylab import *

from mx.DateTime import DateTime, TimeDelta

sys.path.append('..')

from QPCore                   import *
from PMC                      import *
from PMCInventory             import *
from PMCData                  import *
from PMCMetadata              import *


POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)
PI_DIST, PI_MAG, PI_PICKED, PI_IDX, PI_MAGFROMDIST, PI_TIME = range(6)


class PMC_JMA( PMC ):
    """
    QuakePy: PMC_JMA
    PMC method for JMA network (Japan)
    """

    # these attributes are redefined for JMA network
    stationsRequiredForTriggering = 4            # 4
    probThreshold                 = 0.999

    # time shift in JST time zone (hours)
    timeZoneShift = 9.0
    
    def __init__( self, stationlist_filename=None, dist_style=None, **kwargs ):
        super( PMC_JMA, self ).__init__( stationlist_filename, dist_style, **kwargs )


    def _calcMagnitudeFromDistance( self, distance ):
        # this function is redefined for JMA network
        return ( 1.73 * math.log10( distance ) )

    
    def _stationToBeUsed( self, station ):
        # this function is redefined for JMA network
        
        # exclude some networks
        unused_networks = ()
        if station.networkCode in unused_networks:
            return False

        return True


    def importStations( self, filename, encoding = 'ascii' ):
        """
        import station list of JMA

        example data (note: ruler shown below is not part of data)

            01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    
            Hokkaido District  ËÌ³¤Æ»ÃÏÊý
            JMA  Station
            code number  Latitude Longitude Height From     To       Seismographs
            -HO-----------------------------------------------------------------------------------
            ABASH2  201  43 50.65  143 51.96   180          20000317 EMT            ÌÖÁö£²        
            AKKESH  522  42 59.92  144 41.55    20                   E93 D93        ¸ü´ß          
            ASHIBE  523  43 30.42  142 13.16   187                   E93 D93        °²ÊÌ          
            ASHORO  518  43 17.96  143 45.94   210                   E93 D93        ½½¾¡Â­´ó      
            BIRAT2  667  42 46.85  142 21.59   170                   E93 D93        Ê¿¼è          
            CHURUI  519  42 36.86  143 21.74   120                   E93 D93        ½½¾¡ÃéÎà      
            ENIWA   509  42 50.56  141 26.77   185                   E93 D93 WV     ·ÃÄí          
            ERIMO   517  42 01.11  143 09.21    40                   E93 D93 WV     ¤¨¤ê¤â        
            FURANO  507  43 09.96  142 35.38   360                   E93 D93        ÆîÉÙÎÉÌî      
            HAKOD2  192  41 50.20  140 46.50    80          19980228 EMT            È¡´Û£²        
            HIYAMA  513  41 40.81  140 03.23    30          19980205 E93 D93        É°»³¾å¥Î¹ñ    
            HOKURY  504  43 44.72  141 43.26   195                   E93 D93        ¶õÃÎËÌÎµ      
            KAMIAS  505  44 07.13  142 35.58   220                   E93 D93 WV     ¾åÀîÄ«Æü      
            KAMIK2  676  43 52.51  142 44.65   350 19971202          E93 D93        ¾åÀî          
            KAMIKA  506  43 48.84  142 50.60   430          19971201 E93 D93        ¾åÀî          
            KAYABE  514  41 53.41  141 01.73    10                   E93 D93        ÅÏÅçÆî³ýÉô    
            KIYOSA 4656  43 41.67  144 32.08   182 20040514 20041124 E93            À¶Î¤À¶Àô      

        station blocks have lines with 86 chars

        Note: from/to times for stations are in JST, have to be converted to UTC
        """
        lines = getQPDataSource( filename )

        stationBlock = False

        lines_processed     = 0
        update_lines        = 0
        lastStaIdx          = -1

        for line in lines:

            # skip all until separator/network code '-NN...' line
            if stationBlock is False:
                if line[0] == '-':
                    locationCode = line[1:3]
                    stationBlock = True
                    print " --- now loading stations from district/institution: %s" % locationCode
                continue
            elif len( line.strip() ) == 0:
                stationBlock = False
                continue
            else:

                # we are in station block

                # NOTE: stations can occur in more than one lines, with different time intervals of operation
                #       station file has no multiple channels, so station.onTime is the same as channel.onTime 

                # network code is preliminary, to be fixed in later version
                # of station list

                # print " processing line: %s" % line.strip()
                networkCode  = 'JMA'
                stationCode  = line[0:6].strip()

                curr_lat_deg = line[12:15].strip()
                curr_lat_min = float( line[16:21] )
                curr_lat     = float( curr_lat_deg ) + ( curr_lat_min / 60.0 )

                curr_lon_deg = line[22:26].strip()
                curr_lon_min = float( line[27:32] )
                curr_lon     = float( curr_lon_deg ) + ( curr_lon_min / 60.0 )

                curr_elev    = float( line[33:38].strip() )

                if len( line[39:47].strip() ) == 0:
                    curr_from_date  = 1900
                    curr_from_month = 1
                    curr_from_days  = 1
                else:
                    curr_from_date  = int( line[39:43] )
                    curr_from_month = int( line[43:45] )
                    curr_from_days  = int( line[45:47] )

                if len( line[48:56].strip() ) == 0:
                    curr_to_date  = 2100
                    curr_to_month = 1
                    curr_to_days  = 1
                else:
                    curr_to_date  = int( line[48:52] )
                    curr_to_month = int( line[52:54] )
                    curr_to_days  = int( line[54:56] )
                    
                startTime = DateTime( curr_from_date, curr_from_month, curr_from_days,
                                      0, 0, 0.0 ) - TimeDelta( self.timeZoneShift )
                endTime   = DateTime( curr_to_date, curr_to_month, curr_to_days,
                                      0, 0, 0.0 ) - TimeDelta( self.timeZoneShift )
                
                # look in station list if network/station combination is already there
                sta_idx = self._getStationIdx( networkCode, stationCode )

                self.stations[sta_idx].networkCode  = networkCode
                self.stations[sta_idx].stationCode  = stationCode
                self.stations[sta_idx].locationCode = locationCode
                
                self.stations[sta_idx].latitude    = curr_lat
                self.stations[sta_idx].longitude   = curr_lon
                self.stations[sta_idx].elevation   = curr_elev
                
                # check if we add to the last created station (no new station created)
                if sta_idx == lastStaIdx:
                    print " adding data to last created station %s w/ index %s in input line %s with %s stations already read" % (
                        stationCode, sta_idx, lines_processed+1, len(self.stations) )
                    curr_channel = self.stations[sta_idx].channels[0]
                    update_lines += 1
                    

                # check if we add to an already created station (no new station created)
                elif sta_idx < len(self.stations)-1:
                    print " adding data to existing station w/ index %s in input line %s with %s stations already read" % (
                        sta_idx, lines_processed+1, len(self.stations) )
                    print "  current station: locationCode %s, stationCode %s" % ( locationCode, stationCode )
                    print " modified station: locationCode %s, stationCode %s" % (
                        self.stations[sta_idx].locationCode, self.stations[sta_idx].stationCode )
                    curr_channel = self.stations[sta_idx].channels[0]
                    update_lines += 1

                else:

                    # new station, create new channel
                    curr_channel = PMCChannel()
                    self.stations[sta_idx].channels.append( curr_channel )
                
                # set channel data
                curr_channel.waveformID = WaveformStreamID( networkCode, stationCode )
                curr_channel.waveformID.locationCode = locationCode
                
                curr_channel.onTime.append( [ QPDateTime( startTime ),
                                              QPDateTime( endTime ) ] )

                # print " added data to station index %s with networkCode %s, stationCode %s" % ( sta_idx, self.stations[sta_idx].networkCode, self.stations[sta_idx].stationCode )

                lastStaIdx = sta_idx
                lines_processed += 1

        # set station-wide station.onTime
        for curr_sta in self.stations:

            # loop over channels, add onTime of all channels to station onTime
            # Note: could be improved by checking for overlapping time intervals
            for curr_ch in curr_sta.channels:
                curr_sta.onTime.extend( curr_ch.onTime )
        
        print " station import complete, read %s station lines, %s update lines, length of station array %s" % (
            lines_processed, update_lines, len( self.stations ) )

        # check if read numbers match
        if ( update_lines + len( self.stations ) ) != lines_processed:
            print " ATTENTION: numbers of processed stations do not match!"


    def getCatalogFiles( self, metadata ):

        metadata.getCatalogFilesAll()

        # startDate and endDate have to be mxDateTime.DateTimeType
        if (     ( metadata.startDate is not None )
             and ( metadata.endDate is not None )
             and isinstance( metadata.startDate, DateTimeType )
             and isinstance( metadata.endDate, DateTimeType ) ):

            localStartDate = metadata.startDate + TimeDelta( self.timeZoneShift )
            localEndDate   = metadata.endDate + TimeDelta( self.timeZoneShift )

            # loop over years from startDate to endDate
            # Note: this is done with local time, so that it matches
            # the local year/month designations in catalog file names
            for time_year in xrange( localStartDate.year, localEndDate.year + 1 ):

                # get catfiles for time period
                for curr_catfile in metadata.catalogFilesAll:

                    # datecode      = re.match( r'^.+(\d{6})(_[ABC])?.deck.Z.bz2$', curr_catfile )
                    datecode      = re.match( r'^.+(\d{6})(_[ABC])?.deck.bz2$', curr_catfile )
                    yearmonth_str = str( datecode.group(1) )
                    year          = int( yearmonth_str[0:4] )
                    month         = int( yearmonth_str[4:6] )

                    if ( year == time_year ):

                        if not (    ( year == localStartDate.year and month < localStartDate.month )
                                 or ( year == localEndDate.year and month > localEndDate.month ) ):

                            metadata.catalogFiles.append( curr_catfile )
                            print " PMC_JMA: added catalog file %s to current catalog" % curr_catfile

    # ---------------------------------------------------------------------------------------------------------
    
    def _assignPicks( self, curr_station_idx, metadata, **kwargs ):
        """
        function to assign picks from EQ catalog in ASCII format
        redefines function in base class: PMC._assignPicks()
        avoids creation of QPCatalog for hypocenter/pick representation

        metadata is a dictionary containing metadata of the computation run

        - metadata.catalogFiles is a list of catalog file names in JMA Deck format
          phase data is taken from these files
          Note that the proper set of catalog files for a time period must be provided by the caller

        kwargs: verbose_level  = 0 - no information
                                 1 - short information
                                 2 - detailed information
        """

        # estimated maximum number of seen or missed events/picks
        # (earthquakes in time of operation) for a station and a year
        pickInfoLengthChunk = 300000
        
        if 'verbose_level' in kwargs.keys():
            verbose_level = kwargs['verbose_level']
        else:
            verbose_level = None

        if len( metadata.catalogFiles ) == 0:
            # metadata.getCatalogFiles()
            metadata.getCatalogFilesAll()

        curr_station  = self.stations[curr_station_idx]
        pickinfo_cols = curr_station.distribution.pickInfoCols

        try:
            locCode = curr_station.locationCode
        except:
            locCode = ''

        if self._stationToBeUsed( curr_station ):

            if verbose_level is not None and ( verbose_level > 1 ):

                print " - station %s, %s %s %s: assign picks" % (
                    curr_station_idx, curr_station.networkCode, locCode, curr_station.stationCode )

            # over events
            picks_used = 0

            # create a new pickInfo with default size
            curr_station.distribution.pickInfo = numpy.zeros( ( pickInfoLengthChunk, pickinfo_cols ),
                                                                dtype=float )

            # loop over events in ASCII catalog files from catfiles list
            for catfile in metadata.catalogFiles:

                # get ASCII lines
                istream = getQPDataSource( catfile, **kwargs )

                # use only first magnitude in origin line (first list entry)
                mag_indices = [ { 'from': 52, 'to': 54, 'magtype_from': 54, 'magtype_to': 55 },
                                { 'from': 55, 'to': 57, 'magtype_from': 57, 'magtype_to': 58 } ]

                line_ctr = 0

                phasesMode     = None
                hypocenterMode = None
                newEventMode   = None
                commentMode    = None

                for line in istream:

                    line_ctr += 1

                    # check if it is hypocenter, comment, pick, or separator line
                    line_type = line[0]

                    if line_type in ( 'J', 'U', 'I' ):

                        # previous mode can be: None, Hypocenter
                        if hypocenterMode is None:
                            newEventMode  = True

                            # init pmc_params
                            pmc_params = { 'time': None, 'lat': None, 'lon': None, 'depth': None, 'mag': None }

                        else:
                            newEventMode = False

                        hypocenterMode = True
                        phasesMode     = False

                    elif line_type == '_':

                        # previous mode can be: Hypocenter, Phases, Comment
                        if ( hypocenterMode is True ) or ( commentMode is True ):

                            # we are in first phase line after hypocenter line(s)
                            pickPositive = False

                        hypocenterMode = False
                        phasesMode     = True

                    elif line_type == 'C':

                        # previous mode can be: Hypocenter, Comment
                        hypocenterMode = False
                        commentMode    = True

                    elif line_type == 'E':
                        phasesMode     = None
                        hypocenterMode = None
                        newEventMode   = None
                        commentMode    = None

                        # END: origin and phases is finished, write to PMC object

                        # use only complete JMA origins
                        # complete: time, magnitude, lat, lon, depth must exist
                        if (     not (    ( pmc_params['lat'] is None )
                                       or ( pmc_params['lon'] is None )
                                       or ( pmc_params['depth'] is None )
                                       or ( pmc_params['time'] is None )
                                       or ( pmc_params['mag'] is None ) )
                             and ( self._eventDuringOntime( pmc_params['time'], curr_station ) ) ):

                            # check if pickInfo needs to be resized
                            if picks_used > curr_station.distribution.pickInfo.shape[0] - 1:
                                curr_station.distribution.pickInfo.resize( ( curr_station.distribution.pickInfo.shape[0] +
                                                                                pickInfoLengthChunk,
                                                                                pickinfo_cols ),
                                                                            refcheck = 0 )

                            distanceXYZ, distanceXY = distanceBetweenPoints( ( pmc_params['lat'],
                                                                               pmc_params['lon'],
                                                                               pmc_params['depth'] ),
                                                                             ( curr_station.latitude,
                                                                               curr_station.longitude,
                                                                               -(curr_station.elevation/1000.0) ),
                                                                              2 )

                            curr_station.distribution.pickInfo[picks_used, PI_DIST] = float( distanceXYZ )
                            curr_station.distribution.pickInfo[picks_used, PI_MAG] = float( pmc_params['mag'] )
                            curr_station.distribution.pickInfo[picks_used, PI_PICKED] = float( pickPositive )
                            curr_station.distribution.pickInfo[picks_used, PI_IDX] = float( picks_used )
                            curr_station.distribution.pickInfo[picks_used, PI_MAGFROMDIST] = float(
                                self._calcMagnitudeFromDistance(distanceXYZ) )

                            # put focal time as decimal year in last column of pickInfo
                            curr_station.distribution.pickInfo[picks_used, PI_TIME] = float( pmc_params['time'].toDecimalYear() )

                            picks_used = picks_used + 1

                        # proceed to next input line
                        continue

                    else:
                        error_str = " PMC_JMA._assignPicks() - illegal catalog line %s: %s" % ( line_ctr, line )
                        raise ValueError, error_str

                    # ------------------------------------------------------------------------------------------

                    if hypocenterMode is True:
                        try:
                            # first get required fields
                            # NOTE: coordinates and depth need not be there
                            curr_source = line[0]

                            curr_year   = int( line[1:5] )
                            curr_month  = int( line[5:7] )
                            curr_day    = int( line[7:9] )

                            curr_hour   = int( line[9:11] )
                            curr_minute = int( line[11:13] )

                            # seconds are given w/o decimal point
                            curr_second_int  = int( line[13:15] )
                            curr_second_frac = line[15:17]

                        except:
                            print " PMC_JMA._assignPicks() - error in hypocenter line %s: %s" % ( line_ctr, line )
                            continue

                        # skip origin line if it is not a JMA origin
                        # note: pmc_params dict can only be filled by JMA origins
                        if curr_source != 'J':
                            continue
                        elif (     ( pmc_params['lat'] is not None )
                                and ( pmc_params['lon'] is not None )
                                and ( pmc_params['depth'] is not None )
                                and ( pmc_params['time'] is not None )
                                and ( pmc_params['mag'] is not None ) ):

                            # check if already a complete JMA origin is there
                            # if so, discard another one, continue
                            continue
                        else:
                            # we have no or an incomplete JMA origin, reset pmc_params and read next origin line
                            pmc_params = { 'time': None, 'lat': None, 'lon': None, 'depth': None, 'mag': None }

                        ## check if time components are well-behaved and fix if necessary
                        # seconds and minutes are sometimes set to 60
                        # possibly hours can be set to 24
                        timeCorrection = fixTimeComponents( curr_hour, curr_minute, curr_second_int )

                        try:
                            curr_second = float( '.'.join( ( str(int(timeCorrection['component'][2])), curr_second_frac ) ) )
                        except:
                            print " PMC_JMA._assignPicks() - illegal time format in hypocenter line %s: %s" % ( line_ctr, line )
                            continue

                        focal_time_utc = DateTime( curr_year,
                                                    curr_month,
                                                    curr_day,
                                                    timeCorrection['component'][0],
                                                    timeCorrection['component'][1],
                                                    curr_second ) - TimeDelta( self.timeZoneShift )

                        focal_time_utc = adjustDateTime( timeCorrection['increaseFlag'], focal_time_utc )

                        # set time only if within start/end date
                        if metadata.startDate is not None and metadata.endDate is not None:

                            if (     focal_time_utc >= metadata.startDate
                                 and focal_time_utc <= metadata.endDate ):
                                pmc_params['time'] = QPDateTime( focal_time_utc )

                        else:
                            pmc_params['time'] = QPDateTime( focal_time_utc )

                        ## get optional fields

                        # latitude is given as degrees + decimal minutes (w/o decimal point)
                        try:
                            curr_lat_deg        = line[22:24].strip()
                            curr_lat_min_int    = line[24:26]
                            curr_lat_min_frac   = line[26:28]
                            curr_lat_min        = float( '.'.join( ( curr_lat_min_int, curr_lat_min_frac ) ) )
                            curr_lat            = float( curr_lat_deg ) + ( curr_lat_min / 60.0 )
                            pmc_params['lat']   = curr_lat

                        except:
                            # go to next event
                            continue

                        # longitude is given as degrees + decimal minutes (w/o decimal point)
                        try:
                            curr_lon_deg        = line[33:36].strip()
                            curr_lon_min_int    = line[36:38]
                            curr_lon_min_frac   = line[38:40]
                            curr_lon_min        = float( '.'.join( ( curr_lon_min_int, curr_lon_min_frac ) ) )
                            curr_lon            = float( curr_lon_deg ) + ( curr_lon_min / 60.0 )
                            pmc_params['lon']   = curr_lon

                        except:
                            # go to next event
                            continue

                        # depth
                        try:
                            curr_depth_int  = line[44:47].strip()
                            curr_depth_frac = line[47:49].strip()

                            # if depth was determined using 'depth slice method', no fraction is there
                            if len( curr_depth_frac ) > 0:
                                curr_depth = float( '.'.join( ( curr_depth_int, curr_depth_frac ) ) )
                            else:
                                curr_depth = float( curr_depth_int )

                            pmc_params['depth'] = curr_depth

                        except:
                            # go to next event
                            continue

                        ## magnitudes: take only first magnitude
                        curr_mag_code = line[mag_indices[0]['from']:mag_indices[0]['to']].strip()

                        # something there?
                        if len( curr_mag_code ) > 0:

                            try:
                                # is magnitude code numeric? (mag. >= -0.9)
                                curr_mag_code_int = int( curr_mag_code )

                                if curr_mag_code_int >= 0:

                                    # mag. >= 0 (F2.1 w/o decimal point)
                                    curr_mag = float( '.'.join( ( curr_mag_code[0], curr_mag_code[1] ) ) )

                                else:

                                    # -1, -2, ..., -9 (-0.9 <= mag. <= 0.1)
                                    curr_mag = float( '.'.join( ( '-0', curr_mag_code[1] ) ) )

                            except:
                                # mag. <= -1.0, use letter code

                                letter_code = { 'A': '-1', 'B': '-2', 'C': '-3' }
                                curr_mag = float( '.'.join( ( letter_code[curr_mag_code[0]], curr_mag_code[1] ) ) )

                            pmc_params['mag'] = curr_mag

                        else:
                            # go to next input line
                            continue

                    elif phasesMode is True:

                        # check if there is a valid (JMA) origin for picks
                        if (    ( pmc_params['lat'] is None )
                                or ( pmc_params['lon'] is None )
                                or ( pmc_params['depth'] is None )
                                or ( pmc_params['time'] is None )
                                or ( pmc_params['mag'] is None ) ):

                            # finished origin block without good origin, skip subsequent phase lines

                            # print " PMC_JMA._assignPicks() - no valid origin to associate phase line, skipping line %s: %s" % ( line_ctr, line )
                            continue

                        # first and second phase type
                        try:
                            curr_sta_code     = line[1:7].strip()
                            curr_phase_first  = line[15:19].strip()
                            curr_phase_second = line[27:31].strip()
                        except:
                            print " PMC_JMA._assignPicks() - error in phase line %s: %s" % ( line_ctr, line )
                            continue

                        # check if station code and phase type match for this line
                        if (     pickPositive is False
                                and ( curr_sta_code == curr_station.stationCode )
                                and (    ( curr_phase_first in ( 'P', 'EP', 'IP' ) )
                                    or ( curr_phase_second in ( 'P', 'EP', 'IP' ) ) ) ):

                            pickPositive = True

                        continue

                    elif commentMode is True:
                        continue

            # ------------------------------------------------------------------------------------------

            # resize pickInfo
            curr_station.distribution.pickInfo.resize( ( picks_used, pickinfo_cols ), refcheck = 0 )

        elif verbose_level is not None and ( verbose_level > 1 ):
            # station intentionally not used
            print " --- station no. %s, %s %s %s intentionally NOT USED" % \
                        ( curr_station_idx, curr_station.networkCode, curr_station.locationCode,
                            curr_station.stationCode )

    # ---------------------------------------------------------------------------------------------------------
            
    def preprocess( self, catalog = None, timePeriod = None ):
        # this function is redefined for JMA network
        
        # - preprocess inventory
        # (1) exclude some networks
        # (2) exclude all stations not operating in given time period
        #     timePeriod is [ startDate, endDate ]
        
        # - preprocess catalog
        # (3) use only JMA locations with complete parameters (time, lat, lon, depth, magnitude)
        # (4) exclude events that have no picks associated
        # (5) exclude all picks that are not 'P' picks: use only 'P', 'EP', 'IP'
        
        self.preprocessInventory( timePeriod )

        # preprocess catalog only if deck files are read into QPCatalog object
        if catalog is not None:
            self.preprocessCatalog( catalog, timePeriod )

            
    def preprocessInventory( self, timePeriod = None ):
        # this function is redefined for JMA network
        
        # (1) exclude all stations not operating in given time period
        #     timePeriod is [ startDate, endDate ]

        if ( not self.stations ):
            raise ValueError, "PMC_JMA::preprocessInventory - no station list loaded"
        
        unused_networks = ()
        
        # over stations
        for curr_sta_idx in reversed( xrange( len(self.stations) ) ):

            # over channels
            for curr_cha_idx in reversed( xrange( len(self.stations[curr_sta_idx].channels) ) ):

                # does station operate in given time period?
                if ( not self._channelOperatingDuringPeriod( self.stations[curr_sta_idx].channels[curr_cha_idx],
                                                             timePeriod,
                                                             self.stations[curr_sta_idx] ) ):
                    del self.stations[curr_sta_idx].channels[curr_cha_idx]
                    
            if len( self.stations[curr_sta_idx].channels ) == 0:
                del self.stations[curr_sta_idx]

        return True

            
    def preprocessCatalog( self, catalog, timePeriod = None ):
        # this function is redefined for JMA network

        # (1) if timePeriod is set, exclude events outside time period
        # (2) use only JMA locations with complete parameters (time, lat, lon, depth, magnitude)
        # (3) exclude events that have no picks associated
        # (4) exclude all picks that are not 'P' picks: use only 'P', 'EP', IP'

        if catalog is None:
            raise ValueError, "PMC_JMA::preprocessCatalog - no valid catalog given"

        if timePeriod is not None:

            if not (      isinstance( timePeriod[0], DateTimeType )
                     and  isinstance( timePeriod[1], DateTimeType ) ):

                error_str = "PMC_JMA::preprocessCatalog - timePeriod has incorrect type"
                raise TypeError, error_str
            else:
                startTime = QPDateTime( timePeriod[0] )
                endTime   = QPDateTime( timePeriod[1] )
                
        event_cnt = len( catalog.eventParameters.event )
        for curr_ev_idx in reversed( xrange( event_cnt ) ):
            
            event = catalog.eventParameters.event[curr_ev_idx]

            # print " analyse event no. %s w/ ID %s" % ( curr_ev_idx, event.publicID )
            
            # get JMA origin and corresponding magnitude
            # check if preferred origin is JMA, if not, skip
            origin = event.getPreferredOrigin()

            # focal time is always present in catalog (ensured in reading method)
            # focal time in catalog and start/end times are UTC
            if timePeriod is not None:

                if (    ( origin.time.value < startTime )
                     or ( origin.time.value > endTime ) ):

                    del catalog.eventParameters.event[curr_ev_idx]
                    del event
                    continue
            
            if not ( origin.creationInfo.agencyID == 'JMA' ):
                
                del catalog.eventParameters.event[curr_ev_idx]
                
                # print " skipped non-JMA event %s, event count is now %s" % ( event.publicID, catalog.size() )
                del event
                continue

            # get corresponding magnitude for JMA origin
            ori_mag = event.getMagnitudes( origin )
            if len( ori_mag ) > 0:
                magnitude = ori_mag[0]
            else:
                # no magnitude, skip event
                del catalog.eventParameters.event[curr_ev_idx]
                
                # print " skipped event w/o magnitude %s, event count is now %s" % ( event.publicID, catalog.size() )
                del event
                continue
            
            # print " event %s, magnitude %s, origin %s no. of picks %s" % ( event.publicID, magnitude.publicID, origin.publicID, len(event.pick) )

            # check if origin is complete
            # - lat, lon, depth
            # focal time is always there and we have already checked for magnitude
            if (    ( ( not hasattr( origin, 'latitude' ) ) or ( origin.latitude is None ) )
                 or ( ( not hasattr( origin, 'longitude' ) ) or ( origin.longitude is None ) )
                 or ( ( not hasattr( origin, 'depth' ) ) or ( origin.depth is None ) ) ):

                del catalog.eventParameters.event[curr_ev_idx]

                # print " skipped event w/o coordinates %s, event count is now %s" % ( event.publicID, catalog.size() )
                del event
                continue

            # check if there are any picks
            if ( ( event.pick is None ) or ( len(event.pick) == 0 ) ):

                del catalog.eventParameters.event[curr_ev_idx]

                # print " skipped event w/o picks %s, event count is now %s" % ( event.publicID, catalog.size() )
                del event
                continue
            
            elif ( ( event.pick is not None ) and ( len(event.pick) > 0 ) ):

                # check for 'P', 'EP', 'IP' picks, delete picks which are not 'P', 'EP', 'IP'
                for curr_pick_idx in reversed( xrange( len(event.pick) ) ):

                    curr_pick         = event.pick[curr_pick_idx]
                    curr_arrivals_idx = origin.getArrivalsIdx( curr_pick )

                    # check if pick has an associated arrival
                    if len( curr_arrivals_idx ) == 0:
                        
                        del event.pick[curr_pick_idx]
                        
                        # print " deleted pick no. %s, id %s since no arrival is associated" %  ( curr_pick_idx, curr_pick.publicID )

                        del curr_pick

                    else:

                        curr_arr_idx   = curr_arrivals_idx[0]
                        curr_arr_phase = origin.arrival[curr_arr_idx].phase.code
                        
                        if curr_arr_phase not in ( 'P', 'EP', 'IP' ):
                            
                            del origin.arrival[curr_arr_idx]
                            del event.pick[curr_pick_idx]

                            # print " deleted %s phase for pick no. %s, id %s" %  ( curr_arr_phase, curr_pick_idx, curr_pick.publicID )

                            del curr_pick
                    
                # loop over picks has finished
                # if event has no more picks, delete event
                if len(event.pick) == 0:

                    del catalog.eventParameters.event[curr_ev_idx]
                    
                    # print " deleted event with no remaining picks %s, event count is now %s" % ( event.publicID, catalog.size() )
                    del event
                    continue

                # print " event was not deleted, proceed to next"
                
        return True


    def _preprocessStationScenario( self ):
        """
        has been redefined for class PMC_JMA derived from PMC
        deselect stations from 'all' list for scenario computation
        """

        ## this is an example of selecting only the HiNet network
        
        #selected_networks = ( 'HN', )
        #self.selectStationsByNetworks( selected_networks )

        return True


    def initAnnotation( self ):
        """
        has been redefined for class PMC_JMA derived from PMC
        set annotation for QPGrid output
        """

        self.annotation.setProperty(
            title = 'Completeness of the JMA Seismic Network',
            creator = ( 'Fabian Euchner and Danijel Schorlemmer', ),
            publisher = ( 'Danijel Schorlemmer and Fabian Euchner', ),
            rights = ( 'Copyright by Danijel Schorlemmer and Fabian Euchner',
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
