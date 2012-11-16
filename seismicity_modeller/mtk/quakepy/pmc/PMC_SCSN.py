# -*- coding: utf-8 -*-
#
# quakepy/pmc/PMC_SCSN.py
# $Id: PMC_SCSN.py 258 2010-03-09 13:35:12Z fab $
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

__version__  = '$Id: PMC_SCSN.py 258 2010-03-09 13:35:12Z fab $'
__revision__ = '$Revision: 258 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import copy
import cPickle
import math
import numpy

import pyRXP
from   xml.sax   import saxutils

from mx.DateTime import DateTime, DateTimeType, TimeDelta

sys.path.append('..')

from QPCore                   import *
from QPCatalog                import *
from PMC                      import *
from PMCInventory             import *
from PMCData                  import *


class PMC_SCSN( PMC ):
    """
    QuakePy: PMC_SCSN
    PMC method for SCSN network (Southern California)
    """

    # these attributes are redefined for SCSN network
    stationsRequiredForTriggering = 4
    probThreshold                 = 0.999

    timeZoneShift = 0.0
    
    def __init__( self, stationlist_filename=None, dist_style=None, **kwargs ):
        super( PMC_SCSN, self ).__init__( stationlist_filename, dist_style, **kwargs )


    def _calcMagnitudeFromDistance( self, distance ):
        # this function is redefined for SCSN network
        return ( -math.log10( 0.3173 * math.exp( -0.00505 * distance ) * math.pow( distance, -1.14 ) ) )
    # _calcMagnitudeFromDistance = Callable( _calcMagnitudeFromDistance )

    
    def _stationToBeUsed( self, station ):
        # this function is redefined for SCSN network
        
        # exclude some networks
        unused_networks = ( 'MX', 'ZY', 'NR', 'RB', 'LI', 'BK' )
        if station.networkCode in unused_networks:
            return False
        
        # go over channels: exclude station if no 'Z' components are there
        hasZ = False
        for channel in station.channels:
            if channel.waveformID.channelCode.endswith('Z'):
                hasZ = True
                break
        
        if not hasZ:
            return False
        
        return True


    def importStations( self, filename, encoding='ascii' ):
        """
        this function is redefined for SCSN network

        example of station file:

        1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        
        CI  POB   EHZ  -- Polly Butte                      33.68699 -116.92402   974 1976/05/31 3000/01/01
        CI  POB   ELZ  -- Polly Butte                      33.68699 -116.92402   974 1982/01/01 1994/06/24
        CI  PSP   EHZ  -- Palm Springs                     33.79371 -116.54968   162 1975/03/12 2006/01/24
        """

        lines = open( filename, 'r' )
        for line in lines:

            # skip lines that are completely whitespace
            if len( line.strip() ) == 0:
                continue
            
            networkCode  = line[0:2].strip()
            stationCode  = line[4:9].strip()
            channelCode  = line[10:13].strip()
            locationCode = line[15:17].strip()

            # we assume that on/off times are UTC
            startTime = DateTime( int(line[77:81]), int(line[82:84]), int(line[85:87]),
                                  0, 0, 0.0 ) - TimeDelta( self.timeZoneShift )
            endTime   = DateTime( int(line[88:92]), int(line[93:95]), int(line[96:98]),
                                  0, 0, 0.0 ) - TimeDelta( self.timeZoneShift )
            
            # look in station list if network/station combination is already there
            sta_idx = self._getStationIdx( networkCode, stationCode )
            
            self.stations[sta_idx].networkCode = networkCode
            self.stations[sta_idx].stationCode = stationCode

            # non-unicode description
            #self.stations[sta_idx].description = line[18:50].strip()

            # if encoding is 'xml', unescape entities
            if ( encoding == 'html' ) or ( encoding == 'xml' ):
                self.stations[sta_idx].description =  saxutils.unescape( unicode( line[18:50].strip(), 'ascii' ) )
            else:
                self.stations[sta_idx].description = unicode( line[18:50].strip(), encoding )
            
            self.stations[sta_idx].latitude    = float(line[51:59].strip())
            self.stations[sta_idx].longitude   = float(line[60:70].strip())
            self.stations[sta_idx].elevation   = float(line[72:77].strip())
            
            # set channel data
            # TODO: check if data is added to already existing channel
            curr_channel = PMCChannel()
            
            curr_channel.waveformID = WaveformStreamID( networkCode, stationCode, channelCode )
            if ( locationCode != '--' ):
                curr_channel.waveformID.locationCode = locationCode
            
            curr_channel.onTime.append( [ QPDateTime( startTime ),
                                          QPDateTime( endTime ) ] )
            self.stations[sta_idx].channels.append( curr_channel )

        # set station-wide station.onTime
        for curr_sta in self.stations:

            # loop over channels, add onTime of all channels to station onTime
            # Note: could be improved by checking for overlapping time intervals
            for curr_ch in curr_sta.channels:
                curr_sta.onTime.extend( curr_ch.onTime )


    def getCatalogFiles( self, metadata ):

        metadata.getCatalogFilesAll()

        # startDate and endDate have to be mxDateTime.DateTimeType
        if (     ( metadata.startDate is not None )
             and ( metadata.endDate is not None )
             and isinstance( metadata.startDate, DateTimeType )
             and isinstance( metadata.endDate, DateTimeType ) ):

            # loop over years from startDate to endDate
            # catalog file dates are UTC
            for time_year in xrange( metadata.startDate.year, metadata.endDate.year + 1 ):

                # get catfiles for time period
                for curr_catfile in metadata.catalogFilesAll:

                    datecode = re.match( r'^.*phase(\d{4}).dat.gz$', curr_catfile )
                    year_str = str( datecode.group(1) )
                    year     = int( year_str[0:4] )

                    if ( year == time_year ):
                        
                        metadata.catalogFiles.append( curr_catfile )
                        print " PMC_SCSN: added catalog file %s to current catalog" % curr_catfile


    def getCatalog( self, catalog_file ):

        qpc = QPCatalog( idstyle='numeric' )
        qpc.importSTPPhase( catalog_file, compression='gz', minimumDataset=True )
        return qpc

    
    def preprocess( self, catalog, timePeriod=None ):
        # this function is redefined for SCSN network
        
        # - preprocess inventory
        # (1) exclude some networks
        # (2) exclude all channels which have no 'Z' component
        # (3) exclude all channels not operating in given time period
        #     timePeriod is [ startDate, endDate ]
        
        # - preprocess catalog
        # (4) exclude all events with no picks
        # (5) exclude all events with magnitude <= 0.0
        
        self.preprocessInventory( timePeriod )
        self.preprocessCatalog( catalog, timePeriod )

            
    def preprocessInventory( self, timePeriod=None ):
        # this function is redefined for SCSN network
        
        # (1) exclude some networks
        # (2) exclude all channels which have no 'Z' component
        # (3) exclude all channels not operating in given time period
        #     timePeriod is [ startDate, endDate ]
        # (4) exclude all events with magnitude <= 0.0
        
        if ( not self.stations ):
            raise ValueError, "PMC_SCSN::preprocessInventory - no station list loaded"
        
        unused_networks = ( 'MX', 'ZY', 'NR', 'RB', 'LI', 'BK' )
        
        # over stations
        for curr_sta_idx in reversed( xrange( len(self.stations) ) ):
            
            # in excluded network?
            if self.stations[curr_sta_idx].networkCode in unused_networks:
                del self.stations[curr_sta_idx]
            else:
                # over channels
                for curr_cha_idx in reversed( xrange( len(self.stations[curr_sta_idx].channels) ) ):
                    
                    # has channel 'Z' component? / does channel operate in given time period?
                    if (    not self.stations[curr_sta_idx].channels[curr_cha_idx].waveformID.channelCode.endswith('Z')
                         or not self._channelOperatingDuringPeriod( self.stations[curr_sta_idx].channels[curr_cha_idx],
                                                                    timePeriod ) ):
                        del self.stations[curr_sta_idx].channels[curr_cha_idx]
                        
                if len(self.stations[curr_sta_idx].channels) == 0:
                    del self.stations[curr_sta_idx]
        
        # unselect station SMTC from Anza network
        self.unselectStation( 'AZ', 'SMTC' )

        return True

            
    def preprocessCatalog( self, catalog, timePeriod=None ):
        # this function is redefined for SCSN network

        # (1) if timePeriod is set, cut catalog according to times
        # (2) exclude events that have no picks associated
        # (3) exclude all events with magnitude <= 0.0
        # (4) exclude all picks that are not 'P' picks
        
        if catalog is None:
            raise ValueError, "PMC_SCSN::preprocessCatalog - no valid catalog given"

        if timePeriod is not None:

            if not (      isinstance( timePeriod[0], DateTimeType )
                     and  isinstance( timePeriod[1], DateTimeType ) ):

                error_str = "PMC_SCSN::preprocessCatalog - timePeriod has incorrect type"
                raise TypeError, error_str
            else:
                catalog.cut( mintime = timePeriod[0], maxtime = timePeriod[1] )

        for curr_ev_idx in reversed( xrange( len(catalog.eventParameters.event) ) ):

            event     = catalog.eventParameters.event[curr_ev_idx]
            magnitude = catalog.eventParameters.event[curr_ev_idx].getPreferredMagnitude()
            origin    = catalog.eventParameters.event[curr_ev_idx].getPreferredOrigin()

            # print " event %s, magnitude %s, origin %s no. of picks %s" % ( event.publicID, magnitude.publicID, origin.publicID, len(event.pick) )

            if ( ( magnitude is None ) or ( magnitude.mag.value <= 0.0 ) ):

                del catalog.eventParameters.event[curr_ev_idx]
                
                #print " deleted event with magnitude <= 0, event count is now %s" % ( len(self.catalog.eventParameters.event) )
                continue
            
            if ( ( event.pick is None ) or ( len(event.pick) == 0 ) ):

                del catalog.eventParameters.event[curr_ev_idx]
                
                #print " deleted event with no picks, event count is now %s" % ( len(self.catalog.eventParameters.event) )
                continue
            
            elif ( ( event.pick is not None ) and ( len(event.pick) > 0 ) ):

                # check for 'P' picks, delete picks which are not 'P'
                # check for 'Z' channels, delete picks that are not from 'Z' channels
                for curr_pick_idx in reversed( xrange( len(event.pick) ) ):
                    curr_pick    = event.pick[curr_pick_idx]
                    curr_arr_idx = origin.getArrivalsIdx( curr_pick )[0]
                    curr_arr     = origin.arrival[curr_arr_idx]

                    if (    ( curr_arr.phase.code != 'P' )
                         or ( not curr_pick.waveformID.channelCode.endswith('Z') ) ):
                        #print " deleted %s phase for pick no. %s, id %s" %  ( curr_arr.phase.code, curr_pick_idx, curr_pick.publicID )

                        del origin.arrival[curr_arr_idx]
                        del event.pick[curr_pick_idx]

                # if event has no more picks, delete event
                if  len(event.pick) == 0:

                    del catalog.eventParameters.event[curr_ev_idx]
                    
                    #print " deleted event with no remaining picks, event count is now %s" % ( len(self.catalog.eventParameters.event) )
                    continue
                
        return True


    def importAliases( self, aliasfile ):
        """
        this function is redefined for SCSN network
        ascii list has 3 columns: networkCode | alias | stationCode
        """
        
        lines = open( aliasfile, 'r' ).read()
        for line in lines.splitlines():

            # skip blank lines
            if len( line.strip() ) == 0:
                continue
            
            ad = line.split()
            self.stationAliases.append( [ ad[0].strip(), ad[1].strip(), ad[2].strip() ] )

            print " alias: %s %s %s" % ( ad[0], ad[1], ad[2] )
            
        return True


    def _preprocessStationScenario( self ):
        """
        has been redefined for class PMC_SCSN derived from PMC
        deselect stations from 'all' list for scenario computation
        """

        ## this is an example of selecting only the Anza network (Anza scenario)
        
        #selected_networks = ( 'AZ', )
        #self.selectStationsByNetworks( selected_networks )

        return True


    def initAnnotation( self ):
        """
        has been redefined for class PMC_SCSN derived from PMC
        set annotation for QPGrid output
        """

        self.annotation.setProperty( title = 'Completeness of the Southern California Seismic Network',
                                     creator = ( 'Fabian Euchner and Danijel Schorlemmer', ),
                                     publisher = ( 'Danijel Schorlemmer and Fabian Euchner', ),
                                     rights = ( 'Copyright by Danijel Schorlemmer and Fabian Euchner',
                                                'Licensed under Creative Commons Attribution-Noncommercial-Share Alike 3.0',
                                                'http://creativecommons.org/licenses/by-nc-sa/3.0/' ),
                                     subject = ( 'seismology', 'earthquake', 'detection probability', 'PMC', 'seismic network', 'completeness', 'magnitude of completeness' ),
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

