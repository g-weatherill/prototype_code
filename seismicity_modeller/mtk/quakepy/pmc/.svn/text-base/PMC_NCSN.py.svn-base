# -*- coding: iso-8859-1 -*-
#
# quakepy/pmc/PMC_NCSN.py
# $Id: PMC_NCSN.py 87 2008-07-23 00:32:13Z fab $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2007 by Fabian Euchner and Danijel Schorlemmer          #
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

__version__  = '$Id: PMC_NCSN.py 87 2008-07-23 00:32:13Z fab $'
__revision__ = '$Revision: 87 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import copy
import cPickle
import math
import numpy

import pyRXP
from   xml.sax   import saxutils

# import matplotlib
# matplotlib.use('PS')
# from pylab import *

from mx.DateTime     import DateTime
from mx.DateTime     import utc
from mx.DateTime     import DateTimeDeltaFromSeconds, DateTimeDelta, DateTimeDeltaFrom, TimeDelta
from mx.DateTime     import RelativeDateTime
from mx.DateTime.ISO import ParseDateTimeUTC
from QPDateTime import *
from TimeQuantity   import TimeQuantity

sys.path.append('..')

from QPCore                   import *
from QPCatalog                import *
from PMC                      import *
from PMCInventory             import *
from PMCData                  import *

POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)

class PMC_NCSN( PMC ):
    """
    QuakePy: PMC_NCSN
    PMC method for NCSN network (Northern California)
    """

    # these attributes are redefined for NCSN network
    stationsRequiredForTriggering = 5
    probThreshold                 = 0.999

    timeZoneShift = 0.0
    
    useDistStyle = 3
    
    def __init__( self, stationlist_filename = None, dist_style = None, **kwargs ):
        super( PMC_NCSN, self ).__init__( stationlist_filename, self.useDistStyle, **kwargs )

    # here we need to check how can we deal with this to respect different distances
    # (direct to epicenter vs. surface to surface distance)
    # for now we just ignore the surface2surface distance...
    def _TODO_calcMagnitudeFromDistance( self, distanceXYZ, distanceXY ):
        # this function is redefined for NCSN network
        # TODO formula correct? [Eaton, 1992]
        if distanceXYZ <=185.3:
           fMagnitude = 0.821* math.log10(distanceXYZ) + 0.00405 * distanceXYZ
        else:
           fMagnitude = 2.55 * math.log10(distanceXYZ)

        if distanceXY < 70:
            fMagnitude += (-0.09*math.sin(0.07*(distanceXY-25)))
        
        return fMagnitude

    def _calcMagnitudeFromDistance( self, distance ):
        # this function is redefined for NCSN network
        # TODO formula correct?
        if distance <=185.3:
           fMagnitude = 0.821* math.log10(distance) + 0.00405 * distance
        else:
           fMagnitude = 2.55 * math.log10(distance)
        
        return fMagnitude
    
    def _stationToBeUsed( self, station ):
        # this function is redefined for NCSN network
        
        # exclude some networks
        used_networks = ( 'BP', 'BK', 'NC', 'NN', 'PG', 'WR', 'CI'  )
        if station.networkCode not in used_networks:
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

    def importStations( self, filename, encoding = 'ascii' ):
        # this function is redefined for NCSN network 
        #the channel code is made part of the station name
        
        lines = open( filename, 'r' )
        for line in lines:
            
            stationCode  = line[0:4].strip()
            networkCode  = line[5:7].strip()
            channelCode  = line[8:11].strip()
            locationCode = line[20:24].strip()
            #stationCode = stationCodeo+channelCode

            if  int(line[112:114].strip()) == 0:
             smon='01'
            else:
             smon=line[112:114]

            if ( len (line[114:116].strip()) == 0):
             sday='28'
            else:
             sday=line[114:116] 
            #print int(line[108:112]),' ',int(line[112:114]),' ',int(sday)
            try:
                Date( int(line[108:112]),int(smon), int(sday) )
            except mx.DateTime.RangeError:
                #print int(line[108:112]), ' ',int(smon), ' ', int(sday)
                startTime1 = Date( int(line[108:112]), int(smon), int(28) )
            else:
                startTime1 = Date( int(line[108:112]), int(smon), int(sday) )
            
            if (  ( int(line[121:123].strip()) == 0 )
               or ( int(line[121:123].strip()) > 12  ) ):
              emon='01'
            else:
              emon=line[121:123]

            if ( len (line[123:125].strip()) == 0):
             eday='01'
            else:
             eday=line[123:125]
            #print stationCode , ' ',int(line[117:121]),' ',int(line[121:123]),' ',int(eday)
            try:
               Date( int(line[117:121]), int(emon), int(eday) )
            except mx.DateTime.RangeError:
                #print int(line[117:121]), ' ',int(emon), ' ', int(eday)
                endTime1 = Date( int(line[117:121]), int(emon), int(01) )
            else:
                endTime1   = Date( int(line[117:121]), int(emon), int(eday) )
            
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
            
            #print line[54:58].strip()
            self.stations[sta_idx].latitude    = float(line[30:32].strip())+float(line[33:40].strip())/60
            self.stations[sta_idx].longitude   = (float(line[41:44].strip())+float(line[45:52].strip())/60)*-1
            if ( len (line[54:58].strip()) == 0):
                self.stations[sta_idx].elevation   = 0
            else:
                self.stations[sta_idx].elevation   = float(line[54:58].strip())
            
            # set channel data
            curr_channel = PMCChannel()
            
            curr_channel.waveformID = WaveformStreamID( networkCode, stationCode, channelCode )
            if ( locationCode != '--' ):
                curr_channel.waveformID.locationCode = locationCode

            # TODO time is avalable in station file:
            # 134-137 4   ON      On time (hr-min GMT) of start date of operation.
            # 139-142  4   OFF     Off time (hr-min GMT) of end date of operation.
            # but there is no explanation if the date is GMT or local
	    startTime = Date( int(line[89:93]), int(line[93:95]), int(line[95:97]), 0, 0, 0.0) - TimeDelta( self.timeZoneShift )
	    endTime = Date( int(line[98:102]), int(line[102:104]), int(line[104:106]), 0, 0, 0.0) - TimeDelta( self.timeZoneShift )
            curr_channel.onTime.append( [ QPDateTime( startTime ), QPDateTime( endTime ) ] )
            self.stations[sta_idx].channels.append( curr_channel )
	    # print "Start: %s-%s-%s\tEnd: %s-%s-%s" % (line[89:93], line[93:95], line[95:97], line[98:102], line[102:104], line[104:106])

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

                    datecode = re.match( r'^.*Phase(\d{4}).*$', curr_catfile )
                    if datecode == None:
                        print 'Cannot match catalog file: %s' % curr_catfile
                    year_str = str( datecode.group(1) )
                    year     = int( year_str[0:4] )

                    if ( year == time_year ):
                        
                        metadata.catalogFiles.append( curr_catfile )
                        print " PMC_NCSN: added catalog file %s to current catalog" % curr_catfile


    def getCatalog( self, catalog_file ):

        qpc = QPCatalog( idstyle='numeric' )
        qpc.importHypoInverse( catalog_file )
        return qpc


    def preprocess( self, catalog, timePeriod = None ):
        # this function is redefined for NCSN network
        
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
            
    def preprocessInventory( self, timePeriod = None ):

        # this function is redefined for NCSN network
        
        # (1) include used networks only
        # (2) exclude all channels which have no 'Z' component
        # (3) exclude all channels not operating in given time period
        #     timePeriod is [ startDate, endDate ]
        # (4) exclude all events with magnitude <= 0.0
        
        if ( not self.stations ):
            raise ValueError, "PMC_NCSN::preprocessInventory - no station list loaded"
        
        used_networks = ( 'BP', 'BK', 'NC', 'NN', 'PG', 'WR', 'CI' )
        
        # over stations
        for curr_sta_idx in reversed( range( len(self.stations) ) ):
            
            # in excluded network?
            if self.stations[curr_sta_idx].networkCode in used_networks:
                # over channels
                for curr_cha_idx in reversed( range( len(self.stations[curr_sta_idx].channels) ) ):

                    # has channel 'Z' component? / does channel operate in given time period?
                    if (    not self.stations[curr_sta_idx].channels[curr_cha_idx].waveformID.channelCode.endswith('Z')
                         or not self._channelOperatingDuringPeriod( self.stations[curr_sta_idx].channels[curr_cha_idx], timePeriod ) ):
                        self.stations[curr_sta_idx].channels.pop( curr_cha_idx )
                        
                if len(self.stations[curr_sta_idx].channels) == 0:
                    self.stations.pop( curr_sta_idx )
                    
            else:
                 self.stations.pop( curr_sta_idx )

        # remove odd Edwards airforce baser station
        #self.unselectStation( 'CI', 'EDW' )

        return True

    def preprocessCatalog( self, catalog, timePeriod = None ):
        # this function is redefined for NCSN network

        # (1) exclude events that have no picks associated
        # (2) exclude all events with magnitude <= 0.0
        # (3) exclude all picks that are not 'P' picks
        # (4) exclude events not in the time period
        
        if catalog is None:
            raise ValueError, "PMC_NCSN::preprocessCatalog - no valid catalog given"

        for curr_ev_idx in reversed( range( len(catalog.eventParameters.event) ) ):
            #print len(catalog.eventParameters.event)
            event     = catalog.eventParameters.event[curr_ev_idx]
            magnitude = catalog.eventParameters.event[curr_ev_idx].getPreferredMagnitude()
            origin    = catalog.eventParameters.event[curr_ev_idx].getPreferredOrigin()

            #print " event %s, magnitude %s, origin %s no. of picks %s" % ( event.publicID, magnitude.mag.value, origin.publicID, len(event.pick) )
            
            #remove events outside time period
            if (mx.DateTime.cmp(origin.time.value.datetime,timePeriod[0],origin.time.value.cmpEpsilon) < 0 ) :
                del catalog.eventParameters.event[curr_ev_idx]
                continue
            if (mx.DateTime.cmp(origin.time.value.datetime,timePeriod[1],origin.time.value.cmpEpsilon) > 0 ) :
                del catalog.eventParameters.event[curr_ev_idx]
                continue

            if ( ( magnitude is None ) or ( magnitude.mag.value <= 0.0 ) ):

                #catalog.eventParameters.event.remove( event )
                del catalog.eventParameters.event[curr_ev_idx]

                #print " deleted event with magnitude <= 0, event count is now %s" % ( len(catalog.eventParameters.event) )
                continue
            
            if ( ( event.pick is None ) or ( len(event.pick) == 0 ) ):

                #catalog.eventParameters.event.remove( event )
                del catalog.eventParameters.event[curr_ev_idx]

                #print " deleted event with no picks, event count is now %s" % ( len(self.catalog.eventParameters.event) )
                continue
            
            elif ( ( event.pick is not None ) and ( len(event.pick) > 0 ) ):

                # check for 'P' picks, delete picks which are not 'P'
                # check for 'Z' channels, delete picks that are not from 'Z' channels
                for curr_pick_idx in reversed( range( len(event.pick) ) ):
                    curr_pick = event.pick[curr_pick_idx]
                    curr_arrivals_idx = origin.getArrivalsIdx( curr_pick )
                    #print curr_pick.publicID
                    #print curr_pick_idx
                    #print curr_pick.waveformID.stationCode
                    #print curr_pick.time.value
                    #print "arr",origin.arrival[curr_pick_idx].pickID
                    curr_arr  = origin.getArrivals( curr_pick )[0]
                    #print curr_arr.phase.code
                    if (    ( curr_arr.phase.code != 'P' )
                         or ( not curr_pick.waveformID.channelCode.endswith('Z') ) ):
                        #print " deleted %s phase for pick no. %s, id %s" %  ( curr_arr.phase.code, curr_pick_idx, curr_pick.publicID )
                        #print "rming",curr_pick.waveformID.stationCode
                        #print "rming id", curr_pick.publicID
                        curr_arr_idx   = curr_arrivals_idx[0]
                        #origin.arrival.remove( curr_arr )
                        del origin.arrival[curr_arr_idx]
                        #event.pick.remove( curr_pick )   #this is not working wrong one gone
                        del event.pick[curr_pick_idx]
                        #for curr_pick_idx1 in reversed( range( len(event.pick) ) ):
                            #curr_pick1 = event.pick[curr_pick_idx1]
                            #print event.pick[curr_pick_idx1].waveformID.stationCode
                            #print event.pick[curr_pick_idx1].time.value
                            #print 'pick',event.pick[curr_pick_idx1].publicID
                            #print 'arr',origin.arrival[curr_pick_idx1].pickID

                # if event has no more picks, delete event
                if  len(event.pick) == 0:

                    #catalog.eventParameters.event.remove( event )
                    del catalog.eventParameters.event[curr_ev_idx]
                    #print " deleted event with no remaining picks, event count is now %s" % ( len(self.catalog.eventParameters.event) )
                    continue
                
        return True
                        
    def importAliases( self, aliasfile ):
        # this function is redefined for NCSN network

        # ascii list has 3 columns: networkCode | alias | stationCode

        lines = open( aliasfile, 'r' ).read()
        for line in lines.splitlines():
            ad = line.split()
            if len ( ad[0].strip() ) < 3: 
               self.stationAliases.append( [ ad[0].strip(), ad[1].strip(), ad[2].strip() ] )
            else:
               self.stationAliases.append( [ "NC", ad[0].strip(), ad[1].strip() ] )
        return True

    def _preprocessStationScenario( self ):
        """
        has been redefined for class PMC_NCSN derived from PMC
        deselect stations from 'all' list for scenario computation
        """

        ## this is an example of selecting only the Anza network (Anza scenario)
        
        #selected_networks = ( 'AZ', )
        #self.selectStationsByNetworks( selected_networks )

        return True

    def initAnnotation( self ):
        """
        has been redefined for class PMC_NCSN derived from PMC
        set annotation for QPGrid output
        """

        self.annotation.setProperty( title = 'Completeness for Northern California Seismic Network',
                                     creator = ( 'Danijel Schorlemmer', 'Fabian Euchner', 'Michael Lewis' ),
                                     publisher = 'Danijel Schorlemmer',
                                     rights = 'Copyright 2008 Danijel Schorlemmer and Fabian Euchner.',
                                     subject = ( 'seismology', 'earthquake', 'detection probability', 'PMC', 'seismic network', 'completeness', 'magnitude of completeness' ),
                                     source = ( 'http://completeness.usc.edu', 'http://quakepy.org' ) )

        # not used at the moment
        # identifier  version  acknowledgment

        # for these we need run-time information
        # date  coverageTemporal  coverageSpatial
        
        return True
    
    
## ---------------------------------------------------------------------------

