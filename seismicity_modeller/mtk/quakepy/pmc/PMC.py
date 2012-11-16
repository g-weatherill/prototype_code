# -*- coding: utf-8 -*-
#
# quakepy/pmc/PMC.py
# $Id: PMC.py 336 2012-07-09 14:10:49Z tyrone $
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

__version__  = '$Id: PMC.py 336 2012-07-09 14:10:49Z tyrone $'
__revision__ = '$Revision: 336 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import time
import cPickle
import math
import numpy

import cStringIO
import pyRXP
from   xml.sax   import saxutils

# import matplotlib
# matplotlib.use('PS')
# from pylab import *

from mx.DateTime     import Date, DateTimeType
from mx.DateTime.ISO import ParseDateTimeUTC

sys.path.append('..')
sys.path.append('.')

from QPCore                   import *
from QPAnnotation             import *
from QPColorMap               import *
from PMCInventory             import *
from PMCMetadata              import *
from PMCData                  import *
from PMCFunctions             import *

from Event                    import Event
from WaveformStreamID         import WaveformStreamID


POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)
PI_DIST, PI_MAG, PI_PICKED, PI_IDX, PI_MAGFROMDIST, PI_TIME = range(6)

PICKINFO_FILE_SUFFIX = 'pickInfo.xml.bz2'
DISTRO_FILE_SUFFIX = 'distro.xml.bz2'

class PMC( QPObject ):
    """
    QuakePy: PMC
    """

    # re-define these in any class derived from PMC
    stationsRequiredForTriggering = 4
    probThreshold                 = 0.999

    # parameters for pickInfo refine
    thresholdPickCtr           = 1000
    timeIntervalScaleFewPicks  = 5.0        # Michael: 10.0
    timeIntervalScaleManyPicks = 20.0       # Michael: 500.0
    timeIntervalMin            = 1.0/365.0  # one day, expressed in years
    ontimeMinPickCtr           = 50         # minimum number of arrivals an ontime must contain
                                            # time intervals containing fewer picks will not be counted

    onTimesAreRefinedFlag      = False

    # set time margin after and before gap in picks to 1 hour (expressed in years)
    timeSafetyMargin = 1.0 / ( 365.0 * 24.0 )

    # default setting: no subnetworks
    subnetworks = None
    subnetwork = None

    def __init__( self, stationlist_filename=None, dist_style=None, **kwargs ):
        """
        kwargs: encoding = 'ascii' | 'xml'
                smooth   = True | false
        """

        encoding = kwargs.get( 'encoding', 'ascii' )
        self.distSmoothing = kwargs.get( 'smooth', True )

        self.distStyle = dist_style
        
        self.stations       = []
        self.stationAliases = []
        
        self.annotation = QPAnnotation()
        self.initAnnotation()

        if stationlist_filename is not None:
            self.importStations( stationlist_filename, encoding )


    def merge( self, T ):
        """
        merge contents of another PMC object with self
        """
        self.stations.extend( T.stations )
        self.stationAliases.extend( T.stationAliases )

        
    def save( self, filename ):
        fh = writeQPData( filename, binary = True )
        try:
            cPickle.dump( self, fh, 2 )
        except cPickle.PickleError:
            raise cPickle.PickleError, "PMC::save - error pickling PMC"
        fh.close()

    # -------------------------------------------------------------------------
    
    def readXML( self, stream ):

        # get XML file locally or from web
        # lines = getQPDataSource( filename, **kwargs ).read()

        # get whole content from stream
        lines = stream.read()

        # check if it is XML
        if not lines.startswith('<?xml'):
            raise IOError, 'PMC::readXML - input stream is not XML'
            
        tree = pyRXP.Parser().parse( lines )
        
        if tree[POS_TAGNAME] != 'PMC':
            raise TypeError, 'PMC::readXML - input stream is not PMC XML'

        self.fromXML( tree )

            
    def writeXML( self, output, **kwargs ):
        """
        serialize PMC object to XML

        kwargs:
            prettyPrint - pretty formatting of output XML
        """
        if isinstance( output, basestring ):
            ostream = writeQPData( output, **kwargs )
        else:
            ostream = output

        if 'prettyPrint' in kwargs and kwargs['prettyPrint'] is False:
            prettyPrint = False
        else:
            prettyPrint = True
            
        if prettyPrint is True:

            # serialize to string stream
            try:
                curr_stream = cStringIO.StringIO()
                self.toXML( 'PMC', curr_stream )
                streamSuccess = True
            except:
                streamSuccess = False
                print "PMC::writeXML - error in StringIO self.toXML()"

            if streamSuccess is True:
                try:
                    xmlPrettyPrint( curr_stream, ostream )
                    return
                except:
                    print "PMC::writeXML - error in xmlPrettyPrint()"

        # write to output stream w/o pretty print
        # fallback if prettify has not succeeded
        try:
            self.toXML( 'PMC', ostream )
        except:
            raise IOError, "PMC::writeXML - error in self.toXML()"

    # -------------------------------------------------------------------------
    
    def fromXML( self, tree, additionalElements = None ):
        
        # XML attributes
        if tree[POS_ATTRS]:
            attr_dict = tree[POS_ATTRS]
            if 'publicID' in attr_dict.keys():
                self.publicID = attr_dict['publicID']
                
        for child in tree[POS_CHILDREN]:
            
            # multiple children of same type
            if child[POS_TAGNAME] == 'PMCStation':
                sta = PMCStation( self.distStyle, smooth = self.distSmoothing )
                sta.fromXML( child, additionalElements )
                self.stations.append( sta )

            elif child[POS_TAGNAME] == 'stationAliases':

                # get values from attributes of stationAlias elements
                for newchild in child[POS_CHILDREN]:
                    if newchild[POS_TAGNAME] == 'stationAlias':

                        # stationAlias element has 3 required attributes
                        attr_dict = newchild[POS_ATTRS]
                        self.stationAliases.append( 
                            [ attr_dict['networkCode'], attr_dict['alias'],
                              attr_dict['stationCode'] ] )


    def toXML( self, tagname, stream ):

        stream.write( '<?xml version="1.0" encoding="utf-8"?>' )
        stream.writelines( [ '<', tagname, ' xmlns="http://quakepy.org/xmlns/pmc/1.0">' ] )

        # multiple children of same type
        if hasattr( self, 'stations' ) and self.stations is not None:
            for sta in self.stations:
                sta.toXML( 'PMCStation', stream )
                
        if hasattr( self, 'stationAliases' ) and self.stationAliases is not None:

            if len( self.stationAliases ) > 0:

                stream.write( '<stationAliases>' )

                for sta_alias in self.stationAliases:
                    stream.write( ''.join ( ( '<stationAlias networkCode="', sta_alias[0],
                                              '" alias="', sta_alias[1],
                                              '" stationCode="', sta_alias[2], '"/>' ) ) )
                stream.write( '</stationAliases>' )

        stream.writelines( [ '</', tagname, '>' ] )
        return True

    ## -----------------------------------------------------------------------

    def assignPicks( self, metadata=None, **kwargs ):
        """
        kwargs: reduce_pmc=True        delete station information after picks have been assigned
                                       and written to disk (default: True)
        """
        
        if len( self.stations ) == 0:
            raise ValueError, "PMC.assignPicks - no station list loaded"

        myStartMemory = memory()
        if ( ( metadata is not None ) and not isinstance( 
            metadata, PMCMetadata ) ):
            raise ValueError, "PMC.assignPicks - metadata must be of type PMCMetadata"

        if self.subnetworks is not None:
            all_stations = []
            for subnet_data in self.subnetworks.values():
                all_stations.extend( subnet_data['stations'] )

        if self.subnetwork is not None:
            print "assigning picks for sub-network %s" % self.subnetwork

        # loop over stations, starting from end of list since we delete from list
        for curr_station_idx in reversed( xrange( len( self.stations ) ) ):

            # check if station is in one of the sub-networks
            if (     self.subnetworks is not None 
                 and self.stations[curr_station_idx].stationCode \
                    not in all_stations ):

                print " SKIP station %s %s: NOT FOUND IN ONE OF THE"\
                    " SUBNETWORKS" % ( curr_station_idx,
                    self.stations[curr_station_idx].stationCode )
                continue

            # set start/end date of catalog, from metadata
            if metadata is not None:
                self.stations[curr_station_idx].catalogStartDate = QPDateTime( 
                    metadata.startDate )
                self.stations[curr_station_idx].catalogEndDate   = QPDateTime( 
                    metadata.endDate )
            
            # assign picks
            if metadata.useStationPicks:
                self._assignPicksFromStationPicksFiles( curr_station_idx, metadata=metadata, **kwargs )
            else:
                self._assignPicksFromCatalog( curr_station_idx, metadata=metadata, **kwargs )

            # if pickInfo is empty, skip station
            if (    self.stations[curr_station_idx].distribution.pickInfo is None 
                 or self.stations[curr_station_idx].distribution.pickInfo.shape[0] == 0 ):
                del self.stations[curr_station_idx]
                continue

            # refine picks
            if metadata is not None and metadata.refinePickInfo is True:

                print "  running refinePicks ..."
                self._refinePicks( curr_station_idx, **kwargs )

                if 'write_picktime_diff' in kwargs and kwargs['write_picktime_diff'] is True:
                    self.writePickTimeDiff( curr_station_idx, metadata, **kwargs )

            # write station with pickInfo to XML
            self.writePickToXML( curr_station_idx, metadata, **kwargs )

            # delete station from PMC object
            if not ( 'reduce_pmc' in kwargs and kwargs['reduce_pmc'] == False ):
                del self.stations[curr_station_idx]
        
        print 'Memory usage: %d (Stations: %d)' % (memory(myStartMemory), len(self.stations) )

    # TODO(tb) beautyfying
    def _assignPicksFromStationPicksFiles( self, curr_station_idx, metadata, **kwargs ):
        # metadata
        verbose_level = int(kwargs.get('verbose_level', 0))
        curr_station = self.stations[curr_station_idx]
        try:
            locCode = curr_station.locationCode
        except:
            locCode = 'NO_LOCATION_CODE'
        full_station_name = '%s.%s (%s)' % (curr_station.networkCode, curr_station.stationCode, locCode)
        if self._stationToBeUsed( curr_station ):
            print '  -- using station: %s' % full_station_name
            stationPicksFile = '%s/picks_%s.%s.dat' % (metadata.stationPicksDir, curr_station.networkCode, curr_station.stationCode)
            try:
                station_pick_data = open('%s' % stationPicksFile)
            except:
                # we do not fall back into catalog mode:
                raise IOError, 'Cannot open precomputed station picks file: %s' % stationPicksFile
            else:
                # TODO remove?
                if verbose_level > 1:
                    print '  --- opened station picks file: %s' % stationPicksFile

            start_year = QPDateTime(metadata.startDate).toDecimalYear()
            end_year = QPDateTime(metadata.endDate).toDecimalYear()
            if verbose_level > 1:
                print '  -- requested time span: %s (%.8f) - %s (%.8f)' % (metadata.startDate, start_year, metadata.endDate, end_year)
            used_line_count = 0
            off_time_line_count = 0
            all_line_count = 0
            valid_line_count = 0
            # TODO(tb) this is somewhat bloated, maybe using split() instead?
            # 2005.32703792 2.441 0.938853930489 0 914.351578184 914.351576425 smi:geonet.org.nz/ori/348516/GROPE
            # TODO(tb) look for a goog regex uri pattern
            myPattern = {'floatPattern': '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)', 'boolIntPattern': '([0|1])', 'uriPattern': '(.*)'}
            dataPattern = re.compile('^{floatPattern} {floatPattern} {floatPattern} {boolIntPattern} {floatPattern} {floatPattern} {uriPattern}\n$'.format(**myPattern))
            # create a base numpy array; this line will be removed when we're done
            self.stations[curr_station_idx].distribution.pickInfo = numpy.array([[PI_DIST, PI_MAG, PI_PICKED, PI_IDX, PI_MAGFROMDIST, PI_TIME]])
            myStartTime = utc()
            for line in station_pick_data:
                all_line_count += 1
                # data = line.split()
                data = dataPattern.match(line)
                try:
                    eventTime = float(data.group(1))
                except:
                    raise ValueError, 'ERROR in station pick file - invalid year in station pickfile line: %s' % line
                else:
                    valid_line_count += 1
                if start_year <= eventTime <= end_year:
                    if self._eventDuringOntime(fromDecimalYear(eventTime), curr_station):
                        used_line_count += 1
                        # print '    data: ', data.groups()
                        # '{time} {magnitude} {magnitudeFromDistance} {picked} {distanceXY} {distanceXYZ} {eventId}\n'
                        magnitude = float(data.group(5))
                        magnitudeFromDistance = float(data.group(9))
                        picked = float(data.group(13))
                        distanceXY = float(data.group(14))
                        distanceXYZ = float(data.group(18))
                        eventId = data.group(22)
                        self.stations[curr_station_idx].distribution.pickInfo = numpy.append(self.stations[curr_station_idx].distribution.pickInfo,
                                    [[distanceXYZ, magnitude, picked, used_line_count, magnitudeFromDistance, eventTime]], 0)
                    else:
                        off_time_line_count += 1
            # we remove the first line
            self.stations[curr_station_idx].distribution.pickInfo = self.stations[curr_station_idx].distribution.pickInfo[1:]
            print '  -- read %d lines (%s), lines:%d used:%d skipped:%d ' % (
                        all_line_count, str( utc() - myStartTime ), ( used_line_count + off_time_line_count),
                        used_line_count, off_time_line_count )

        else:
            print '  -- not using station: %s' % full_station_name
        return True

    def _assignPicksFromCatalog( self, curr_station_idx, metadata, **kwargs ):
        """
        kwargs: reduce_catalog = True - delete already assigned picks in catalog (default: False)
                verbose_level  = 0 - no information
                                 1 - short information
                                 2 - detailed information
        """

        verbose_level = int(kwargs.get('verbose_level', 0))

        curr_station  = self.stations[curr_station_idx]
        pickinfo_cols = self.stations[curr_station_idx].distribution.pickInfoCols

        try:
            locCode = curr_station.locationCode
        except:
            locCode = ''

        if self._stationToBeUsed( curr_station ):

            if verbose_level > 1:

                print " - station %s, %s %s %s: assign picks" % (
                    curr_station_idx, curr_station.networkCode, locCode, curr_station.stationCode )
                myStartTime = utc()

            # loop over catalog files in metadata
            for curr_catfile in metadata.catalogFiles:

                catalog = self.getCatalog( curr_catfile )
                self.preprocessCatalog( catalog, [ metadata.startDate, metadata.endDate ] )

                if verbose_level > 0:
                    print " -- assign picks: number of events = %s" % catalog.size()

                # look if we append to an existing pickInfo
                if self.stations[curr_station_idx].distribution.pickInfo is not None:
                    # resize existing pickInfo, add as many rows as there are events

                    # get number of rows in existing pickInfo, add no. of events in current catalog
                    old_rows = self.stations[curr_station_idx].distribution.pickInfo.shape[0]
                    new_rows = old_rows + len(catalog.eventParameters.event)
                    self.stations[curr_station_idx].distribution.pickInfo.resize( 
                        (new_rows, pickinfo_cols), refcheck = 0 )
                else:
                    # create a new pickInfo
                    old_rows = 0
                    new_rows = len(catalog.eventParameters.event)
                    self.stations[curr_station_idx].distribution.pickInfo = numpy.zeros( 
                        (new_rows, pickinfo_cols), dtype=numpy.float )

                picks_used = 0
                
                # only used well-behaved events (ensure this in preprocessCatalog)
                for curr_event_idx, curr_event in enumerate( catalog.eventParameters.event ):

                    # check if event is in one of the stations' "on" time windows
                    if self._eventDuringOntime( curr_event, curr_station ):

                        curr_ori = curr_event.getPreferredOrigin()

                        #distanceXYZ, distanceXY = distanceBetweenPoints( ( curr_ori.latitude.value,
                                                                            #curr_ori.longitude.value,
                                                                            #curr_ori.depth.value ),
                                                                            #( curr_station.latitude,
                                                                            #curr_station.longitude,
                                                                            #-(curr_station.elevation/1000.0) ),
                                                                            #2 )

                        distanceXYZ = distanceBetweenPointsWGS84( ( curr_ori.latitude.value,
                                                                    curr_ori.longitude.value,
                                                                    curr_ori.depth.value ),
                                                                  ( curr_station.latitude,
                                                                    curr_station.longitude,
                                                                    -(curr_station.elevation/1000.0) ) )

                        distanceXY = distanceBetweenPointsWGS84( ( curr_ori.latitude.value,
                                                                    curr_ori.longitude.value,
                                                                    0 ),
                                                                  ( curr_station.latitude,
                                                                    curr_station.longitude,
                                                                    0 ) )

                        magnitude = curr_event.getPreferredMagnitude()
                        if magnitude is None:
                            raise ValueError, "PMC._assignPicks - event has no preferred magnitude"

                        self.stations[curr_station_idx].distribution.pickInfo[old_rows + picks_used, PI_DIST] = float( distanceXYZ )
                        self.stations[curr_station_idx].distribution.pickInfo[old_rows + picks_used, PI_MAG,] = float( magnitude.mag.value )
                        self.stations[curr_station_idx].distribution.pickInfo[old_rows + picks_used, PI_PICKED] = float(
                            self._eventPicked(curr_event, curr_station.networkCode, curr_station.stationCode) )
                        self.stations[curr_station_idx].distribution.pickInfo[old_rows + picks_used, PI_IDX] = float( curr_event_idx )
                        self.stations[curr_station_idx].distribution.pickInfo[old_rows + picks_used, PI_MAGFROMDIST] = float(
                            self._calcMagnitudeFromDistance(distanceXYZ) )

                        # put focal time as decimal year in last column of pickInfo
                        self.stations[curr_station_idx].distribution.pickInfo[old_rows + picks_used, PI_TIME] = float(
                            curr_ori.time.value.toDecimalYear() )

                        picks_used = picks_used + 1

                    # delete all picks that belong to the current station
                    if ( 'reduce_catalog' in kwargs and kwargs['reduce_catalog'] == True ):
                        self._deletePicks( curr_event, curr_station.networkCode, curr_station.stationCode )

                del catalog

                # resize pickInfo
                self.stations[curr_station_idx].distribution.pickInfo.resize( ( old_rows + picks_used,
                                                                                pickinfo_cols ),
                                                                                refcheck = 0 )

                if verbose_level > 1:
                    print " -- time passed: %s" % ( str ( utc() - myStartTime ) )
                    if picks_used > 1:
                        print " ----- catalog file %s: %s picks / total pickInfo size is %s" % (
                        	os.path.basename( curr_catfile ), picks_used, old_rows + picks_used )

        elif verbose_level > 1:
            # station intentionally not used
            print " --- station no. %s, %s %s %s intentionally NOT USED" % \
                        ( curr_station_idx, curr_station.networkCode, locCode, curr_station.stationCode )
        return True


    def _refinePicks( self, station_idx, **kwargs ):

        verbose_level = int(kwargs.get('verbose_level', 0))
            
        self.onTimesAreRefinedFlag = True
        
        curr_station = self.stations[station_idx]

        try:
            locCode = curr_station.locationCode
        except:
            locCode = ''
                
        if self._stationToBeUsed( curr_station ):

            currStart         = []
            currEnd           = []
            missedPickIndices = []

            # sort pickInfo, ascending in time
            indices = self.stations[station_idx].distribution.pickInfo[:, PI_TIME].argsort( axis=0 )
            self.stations[station_idx].distribution.pickInfo = \
                self.stations[station_idx].distribution.pickInfo[indices]

            # create temp array for positive picks
            sel          = ( self.stations[station_idx].distribution.pickInfo[:, PI_PICKED] == 1.0 )
            seenIndices  = numpy.arange( self.stations[station_idx].distribution.pickInfo.shape[0] )
            seenArrivals = numpy.hstack(
                [ seenIndices[sel][:, numpy.newaxis],
                  self.stations[station_idx].distribution.pickInfo[sel, PI_TIME][:, numpy.newaxis] ] )

            # initialize selector for events to be kept
            keepEventSelector = numpy.ones( 
                self.stations[station_idx].distribution.pickInfo.shape[0],
                    dtype=numpy.bool ) * True



            if verbose_level > 1:
                print "  station %s, %s %s %s: %s events, %s picked between %s, %s (%s yr)" % (
                                                station_idx,
                                                curr_station.networkCode,
                                                locCode,
                                                curr_station.stationCode,
                                                self.stations[station_idx].distribution.pickInfo.shape[0],
                                                seenArrivals.shape[0],
                                                mxDateTime2ISO( fromDecimalYear(
                                                    self.stations[station_idx].distribution.pickInfo[0, PI_TIME] ),
                                                    showsecfrac=False ),
                                                mxDateTime2ISO( fromDecimalYear(
                                                    self.stations[station_idx].distribution.pickInfo[-1, PI_TIME] ),
                                                    showsecfrac=False ),
                                                "%.3f" % (
                                                    self.stations[station_idx].distribution.pickInfo[-1, PI_TIME] -\
                                                    self.stations[station_idx].distribution.pickInfo[0, PI_TIME] ) )

            if seenArrivals.shape[0] > 1:

                curr_station.onTimesAreRefinedFlag = True

                # get mean time difference between two events that have been registered
                self.stations[station_idx].meanInterPickTime = 0.0
                self.stations[station_idx].interPickTime     = []
                
                for idx in xrange( seenArrivals.shape[0]-1 ):
                    curr_diff = 365.0 * ( seenArrivals[idx+1, 1] - seenArrivals[idx, 1] )

                    # append to interPickTime in days
                    self.stations[station_idx].interPickTime.append( curr_diff )
                    self.stations[station_idx].meanInterPickTime += curr_diff

                self.stations[station_idx].meanInterPickTime = self.stations[station_idx].meanInterPickTime / (
                    seenArrivals.shape[0]-1 )

                if verbose_level > 1:
                    print "    mean time between arrivals: %.5f (hours)" % (
                        24.0 * self.stations[station_idx].meanInterPickTime )

                # get intervals that are larger than scaled mean time between two arrivals
                # scale factor depends on number of arrivals
                if seenArrivals.shape[0] < self.thresholdPickCtr:
                    scaleFactor = self.timeIntervalScaleFewPicks
                else:
                    scaleFactor = self.timeIntervalScaleManyPicks

                self.stations[station_idx].timeIntervalScale = scaleFactor
                self.stations[station_idx].ontimeMinPickCtr  = self.ontimeMinPickCtr

                # compute gap in pick times above which a station off-time will be assumed
                # this is the mean inter-pick time multiplied with a factor
                # that depends on the number of picks recorded at station
                # few picks recorded: gaps are already large, small factor
                # many picks recorded: gaps are small, large factor
                gap_limit = scaleFactor * self.stations[station_idx].meanInterPickTime / 365.0

                # if mean inter-pick time is so small that critical gap becomes smaller that threshold
                # value (default: one day), use this threshold value
                if gap_limit < self.timeIntervalMin:
                    gap_limit = self.timeIntervalMin
                        
                # set last ontime interval end after last pick
                currEnd.append( [ seenArrivals[-1, 1] + self.timeSafetyMargin,
                                  int( seenArrivals[-1, 0] ) ] )

                # loop over seen arrivals, starting from end of array (most recent events)
                for idx in reversed( xrange( 1, seenArrivals.shape[0] ) ):

                    # get gaps in pick
                    if ( ( seenArrivals[idx, 1] - seenArrivals[idx-1, 1] ) > gap_limit ):

                        # gap detected
                        # 1) set end time of current 'on' interval (+ 1 hr)
                        # 2) set start time of next 'on' interval (-1 hr)
                        # 3) delete all negative picks in-between

                        currStart.append( [ seenArrivals[idx, 1] - self.timeSafetyMargin,
                                            int( seenArrivals[idx, 0] ) ] )
                        currEnd.append( [ seenArrivals[idx-1, 1] + self.timeSafetyMargin,
                                          int( seenArrivals[idx-1, 0] ) ] )

                        # add indices of missed picks in pickInfo that fall in gap to index vector
                        # upper index is one below index of upper limit of gap
                        # lower index is one above index of lower limit of gap
                        cutIndices = range( int( seenArrivals[idx, 0] ) -1,
                                            int( seenArrivals[idx-1, 0] ),
                                            -1 )
                        missedPickIndices.extend( cutIndices )

                        if verbose_level > 1:
                            print "      -> found a gap between %s and %s, contains %s missed events in %.2f days" % (
                                mxDateTime2ISO( fromDecimalYear( seenArrivals[idx-1, 1] ), showsecfrac=False ),
                                mxDateTime2ISO( fromDecimalYear( seenArrivals[idx, 1] ), showsecfrac=False ),
                                len( cutIndices ),
                                365.0 * ( seenArrivals[idx, 1] - seenArrivals[idx-1, 1] ) )

                # set first ontime interval start before first pick
                currStart.append( [ seenArrivals[0, 1] - self.timeSafetyMargin,
                                    int( seenArrivals[0, 0] ) ] )

                # set new onTimes in station
                self.stations[station_idx].onTimeRefined = []
                for idx in reversed( xrange( len( currStart ) ) ):
                    
                    pick_in_interval_ctr = currEnd[idx][1] - currStart[idx][1] + 1
                    enoughPicksFlag = False

                    intervalStartTime = QPDateTime( fromDecimalYear( currStart[idx][0] ) )
                    intervalEndTime   = QPDateTime( fromDecimalYear( currEnd[idx][0] ) )
                    
                    # check if minimum of picks is in time interval
                    if pick_in_interval_ctr >= self.ontimeMinPickCtr:

                        self.stations[station_idx].onTimeRefined.append( [ intervalStartTime,
                                                                           intervalEndTime ] )
                        enoughPicksFlag = True

                    if verbose_level > 1:

                        interval_data_str = "    on time %s: %s, %s (%s picks in %.2f days)" % (
                                                    idx,
                                                    intervalStartTime.toISO( showsecfrac=False ),
                                                    intervalEndTime.toISO( showsecfrac=False ),
                                                    pick_in_interval_ctr,
                                                    365.0 * ( currEnd[idx][0] - currStart[idx][0] ) )
                                                    
                        if enoughPicksFlag is True:
                            print interval_data_str
                        else:
                            print ' '.join( ( interval_data_str, "IGNORED: TOO FEW PICKS" ) )
                            

                missedPickIndices.reverse()

                if len( missedPickIndices ) > 0:

                    # delete missed picks in pickInfo that fall in gap
                    keepEventSelector[missedPickIndices] = False
                    self.stations[station_idx].distribution.pickInfo = \
                        self.stations[station_idx].distribution.pickInfo[keepEventSelector,:]

                    if verbose_level > 1:
                        print "    after removing %s events in gaps: %s events remain" % (
                            len( missedPickIndices ),
                            self.stations[station_idx].distribution.pickInfo.shape[0] )

    ## -----------------------------------------------------------------------

    def writePickTimeDiff( self, station_idx, metadata, **kwargs ):

        # write time difference between picks to ASCII file
        curr_station = self.stations[station_idx]

        try:
            locCode = curr_station.locationCode
        except:
            locCode = ''
            
        if self._stationToBeUsed( curr_station ):

            # get positive picks
            sel          = ( self.stations[station_idx].distribution.pickInfo[:, PI_PICKED] == 1.0 )
            seenArrivals = self.stations[station_idx].distribution.pickInfo[sel, PI_TIME]
                  
            picks_used   = seenArrivals.shape[0]

            if picks_used > 1:

                # open file
                if len( locCode ) > 0:
                    pickdiff_filename = '.'.join( ( curr_station.networkCode, locCode, curr_station.stationCode,
                                                    'pickTimeDiff.dat.gz' ) )
                else:
                    pickdiff_filename = '.'.join( ( curr_station.networkCode, curr_station.stationCode,
                                                    'pickTimeDiff.dat.gz' ) )
                pickdiff_file     = os.path.join( metadata.pickTimeDiffDir, pickdiff_filename )

                fh = writeQPData( pickdiff_file, compression='gz' )
                
                for idx in xrange( picks_used-1 ):

                    fh.write( str( 365.0 * ( seenArrivals[idx+1] - seenArrivals[idx] ) ) + '\n' )

                fh.close()

    
    def writePicksToXML( self, metadata, **kwargs ):

        for curr_station_idx in reversed( xrange( len( self.stations ) ) ):
            self.writePickToXML( curr_station_idx, metadata, **kwargs )
            
        
    def writePickToXML( self, station_idx, metadata, **kwargs ):

        verbose_level = kwargs.get('verbose_level', 0)
            
        # write pickInfo to XML
        curr_station = self.stations[station_idx]

        try:
            locCode = curr_station.locationCode
        except:
            locCode = ''
            
        if self._stationToBeUsed( curr_station ):

            picks_used    = curr_station.distribution.pickInfo.shape[0]
            pickinfo_cols = curr_station.distribution.pickInfoCols

            if picks_used > 0:

                if (     ( metadata.discardPickInfo is True )
                     and ( verbose_level > 1 ) ):

                    print " --- assigned picks for station %s, %s %s %s, with %s picks / NOT writing to XML file" % \
                        ( station_idx, curr_station.networkCode, locCode,
                          curr_station.stationCode, picks_used )

                else:

                    ## write station's pick info to disk
                    if len( locCode ) > 0:
                        sta_filename_xml = '.'.join( ( 
                            curr_station.networkCode, locCode, 
                            curr_station.stationCode,
                            PICKINFO_FILE_SUFFIX ) )
                    else:
                        sta_filename_xml = '.'.join( ( 
                            curr_station.networkCode, 
                            curr_station.stationCode,
                            PICKINFO_FILE_SUFFIX ) )
                    sta_file_xml     = os.path.join( metadata.pickInfoDir, sta_filename_xml )

                    if verbose_level > 1:

                        print " --- assigned picks for station %s, %s %s %s, with %s picks / writing to XML file %s" % \
                            ( station_idx, curr_station.networkCode, locCode,
                              curr_station.stationCode, picks_used, sta_filename_xml )

                    curr_stream = cStringIO.StringIO()

                    curr_stream.write( '<?xml version="1.0" encoding="utf-8"?>' )
                    curr_stream.write( '<PMC xmlns="http://quakepy.org/xmlns/pmc/1.0">' )
                    curr_station.toXML( 'PMCStation', curr_stream )
                    curr_stream.write( '</PMC>' )

                    xmlPrettyPrint( curr_stream, writeQPData( sta_file_xml, compression='bz2' ) )

                # if desired, plot pick info (raw distribution)
                if metadata.plotPickInfo is True:

                    if verbose_level > 1:
                        print " --- creating plot of pickInfo"

                    self.plotRawDistribution( curr_station, metadata )

                # if desired, do fillup (compute probdistro), write to XML file, and plot distro
                if metadata.combinePickInfoProbDistro is True:

                    # fillup
                    curr_station.distribution.restoreDistStyle( metadata.useDistStyle )
                    curr_station.distribution.setSmooth( metadata.smoothDist )

                    if verbose_level > 1:
                        print " --- now calling fillup ..."

                    # compute prob distro
                    curr_station.fillup( self._calcMagnitudeFromDistance )

                    # write prob distro do disk
                    self.writeProbDistro( station_idx, metadata, **kwargs )

            elif verbose_level > 1:
                print " --- assigned picks for station %s, %s %s %s / NO PICKS FOUND" % \
                        ( station_idx, curr_station.networkCode, locCode,
                          curr_station.stationCode )


    def writeProbDistros( self, metadata, **kwargs ):

        for curr_station_idx in reversed( xrange( len( self.stations ) ) ):
            self.writeProbDistro( curr_station_idx, metadata, **kwargs )

            
    def writeProbDistro( self, station_idx, metadata, **kwargs ):

        # TODO(tb) verbose_level handling
        if 'verbose_level' in kwargs:
            verbose_level = kwargs['verbose_level']
        else:
            verbose_level = None

        if not isinstance( metadata, PMCMetadata ):
            raise ValueError, "PMC.writeProbDistro() - metadata must be of type PMCMetadata"

        curr_station = self.stations[station_idx]

        try:
            locCode = curr_station.locationCode
        except:
            locCode = ''
            
        if self._stationToBeUsed( curr_station ):

            # write distro to disk, if not empty
            if (     hasattr( curr_station.distribution, 'distro' )
                 and curr_station.distribution.distro is not None
                 and curr_station.distribution.distro.shape[0] > 0 ):

                # get filename for probdistro XML
                if len( locCode ) > 0:
                    distro_filename_xml = '.'.join( ( 
                        curr_station.networkCode, locCode, 
                        curr_station.stationCode, DISTRO_FILE_SUFFIX ) )
                else:
                    distro_filename_xml = '.'.join( ( 
                        curr_station.networkCode, curr_station.stationCode,
                        DISTRO_FILE_SUFFIX ) )
                distro_file_out     = os.path.join( metadata.distroDir, distro_filename_xml )

                # delete pickInfo, if not already empty
                if (     hasattr( curr_station.distribution, 'pickInfo' )
                     and curr_station.distribution.pickInfo is not None ):
                    del curr_station.distribution.pickInfo

                if verbose_level is not None and ( verbose_level > 1 ):

                    print " --- writing prob distro XML file %s for station %s, %s %s %s" % \
                        ( distro_filename_xml, station_idx, curr_station.networkCode,
                          locCode, curr_station.stationCode )

                curr_stream = cStringIO.StringIO()

                curr_stream.write( '<?xml version="1.0" encoding="utf-8"?>' )
                curr_stream.write( '<PMC xmlns="http://quakepy.org/xmlns/pmc/1.0">' )
                curr_station.toXML( 'PMCStation', curr_stream )
                curr_stream.write( '</PMC>' )

                xmlPrettyPrint( curr_stream, writeQPData( distro_file_out, compression='bz2' ) )

                # if desired, plot prob distro
                if metadata.plotProbDistro is True:

                    if verbose_level is not None and ( verbose_level > 1 ):
                        print " --- creating plot of prob distro"

                    self.plotDistribution( curr_station,
                                           metadata,
                                           colortable=metadata.probDistroColorMapFile )

                # if no grid is computed, delete distro
                if metadata.skipGrid is True:
                    del curr_station.distribution.distro

            elif verbose_level is not None and ( verbose_level > 1 ):
                print " --- station %s, %s %s %s / PROB DISTRO EMPTY" % \
                        ( station_idx, curr_station.networkCode, locCode,
                          curr_station.stationCode )

    ## -----------------------------------------------------------------------
    
    def importStations( self ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return True


    def preprocess( self ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return True

    
    def preprocessInventory( self ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return True

    
    def preprocessCatalog( self ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return True

    
    def postprocess( self ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return True

    
    def importAliases( self ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return True

    
    def computePMC( self ):
        # uses methods that have been redefined in derived classes
        return True

    
    def initAnnotation( self ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return True

    ## -----------------------------------------------------------------------
    
    def plotRawDistribution( self, station, metadata, imgfile=None, **kwargs ):
        """

        default colors for plotting events:
          not picked: grey
          picked:     black
          
        kwargs:
            imgformat: 'png' - convert to png, 'eps' no conversion
            plot_missed_color: Matplotlib color as string specification
            plot_picked_color: Matplotlib color as string specification
        """
        imgformat  = 'png'

        plot_missed_color = '0.6'
        plot_picked_color = 'k'
        
        if 'plot_missed_color' in kwargs and kwargs['plot_missed_color'] is not None:
            plot_missed_color = kwargs['plot_missed_color']

        if 'plot_picked_color' in kwargs and kwargs['plot_picked_color'] is not None:
            plot_picked_color = kwargs['plot_picked_color']

        ax = subplot(111)
        
        abscissa_picked    = []
        abscissa_notpicked = []
        ordinate_picked    = []
        ordinate_notpicked = []
        
        for curr_pick in station.distribution.pickInfo:
            if curr_pick[PI_PICKED]:
                abscissa_picked.append( curr_pick[PI_MAG] )
                ordinate_picked.append( curr_pick[PI_DIST] )
            else:
                abscissa_notpicked.append( curr_pick[PI_MAG] )
                ordinate_notpicked.append( curr_pick[PI_DIST] )

        plot( abscissa_notpicked, ordinate_notpicked, color = plot_missed_color, marker = '.',
              linestyle = 'None', zorder = 1 )
        plot( abscissa_picked, ordinate_picked, color = plot_picked_color, marker = '.',
              linestyle = 'None', zorder = 2 )
        
        xlabel( 'Magnitude' )
        ylabel( 'Distance / km' )
        title( "Network " + station.networkCode + ", Station " + station.stationCode )
        
        # get filename
        if imgfile is None:

            # get filename from station name
            if hasattr( station, 'locationCode' ) and station.locationCode is not None:
                imgfile_name = '.'.join( ( station.networkCode, station.locationCode, station.stationCode, 'pickInfo' ) )
            else:
                imgfile_name = '.'.join( ( station.networkCode, station.stationCode, 'pickInfo' ) )

            
        if imgformat in ( 'eps', 'EPS' ):
            
            if imgfile is None:
                imgfile = os.path.join( metadata.pickInfoPlotDir, '.'.join( ( imgfile_name, 'eps' ) ) )
            
            epsfile = imgfile
            
        else:
            
            epsfile = os.path.join( metadata.pickInfoPlotDir, 'temp.eps' )

            if imgfile is None:
                imgfile = os.path.join( metadata.pickInfoPlotDir, '.'.join( ( imgfile_name, 'png' ) ) )

        savefig( epsfile )
        close()

        # convert to png
        if imgformat in ( 'png', 'PNG' ):
            os.system( 'convert -density 144 %s %s' % ( epsfile, imgfile ) )

            # remove EPS file
            os.remove( epsfile )


    def plotDistribution( self, station, metadata, imgfile=None, **kwargs ):
        """

        kwargs:
            imgformat: 'png' - convert to png, 'eps' no conversion
            colortable: has to be a colortable file in povray format
        """

        imgformat  = 'png'
        colortable = None

        if 'colortable' in kwargs and kwargs['colortable'] is not None:
            colortable = kwargs['colortable']
            
        if 'imgformat' in kwargs and kwargs['imgformat'] is not None:
            imgformat = kwargs['imgformat']
            
        ax = subplot(111)

        # load colortable
        # colordict = gmtColormap( colortable )

        cm = QPColorMap()
        cm.importPovray( colortable, flip=True, extend=True )

        my_cmap = cm.getMatplotlib()
        
        xysteps = station.distribution.distro_binning[self.distStyle-1]
        Xsteps  = xysteps[0]
        Ysteps  = xysteps[1]

        X1 = range( 0, int( (Xsteps[2]-Xsteps[0]) / Xsteps[1]+2 ), 1 )
        Y1 = range( 0, int( (Ysteps[2]-Ysteps[0]) / Ysteps[1]+2 ), 1 )

        X = numpy.zeros( ( len(X1), 1 ), dtype=numpy.float )
        for indx, value in enumerate(X1):
            X[indx] = ( value * Xsteps[1] ) + Xsteps[0] - Xsteps[1]/2

        Y = numpy.zeros( ( len(Y1), 1 ), dtype=numpy.float )
        for indx, value in enumerate(Y1):
            Y[indx]= ( value * Ysteps[1] ) + Ysteps[0] - Ysteps[1]/2

        x, y = meshgrid(X,Y)

        pcolor( x, y, station.distribution.distro, cmap = my_cmap, vmin = 0, vmax = 1 )
        colorbar()
        
        xlabel( 'Magnitude' )
        ylabel( 'Distance / km' )
        title( "Network " + station.networkCode + ", Station " + station.stationCode )
        axis( 'tight' )

        # get filename
        if imgfile is None:

            # get filename from station name
            if hasattr( station, 'locationCode' ) and station.locationCode is not None:
                imgfile_name = '.'.join( ( station.networkCode, station.locationCode, station.stationCode, 'distro' ) )
            else:
                imgfile_name = '.'.join( ( station.networkCode, station.stationCode, 'distro' ) )

            
        if imgformat in ( 'eps', 'EPS' ):
            
            if imgfile is None:
                imgfile = os.path.join( metadata.distroPlotDir, '.'.join( ( imgfile_name, 'eps' ) ) )
            
            epsfile = imgfile
            
        else:
            
            epsfile = os.path.join( metadata.distroPlotDir, 'temp.eps' )

            if imgfile is None:
                imgfile = os.path.join( metadata.distroPlotDir, '.'.join( ( imgfile_name, 'png' ) ) )

        savefig( epsfile )
        close()

        # convert to png
        if imgformat in ( 'png', 'PNG' ):
            os.system( 'convert -density 144 %s %s' % ( epsfile, imgfile ) )

            # remove EPS file
            os.remove( epsfile )

    ## -----------------------------------------------------------------------
    
    def _eventDuringOntime( self, event, station, channel=None ):
        """
        event can either be
         - of type Event (default)
         - event time of type mx.DateTime or QPDateTime
           (if catalog is not read into QPCatalog object, e.g., for Japanese catalog)
        """

        if isinstance( event, Event ):
            event_time = event.getPreferredOrigin().time.value
        elif isinstance( event, DateTimeType ):
            event_time = QPDateTime( event )
        elif isinstance( event, QPDateTime ):
            event_time = event
        else: 
            raise ValueError, "PMC._eventDuringOntime - illegal type for input parameter event"
            
        # check if onTime has been refined
        if (     self.onTimesAreRefinedFlag is True
             and station.onTimesAreRefinedFlag is True ):

            for curr_ontime in station.onTimeRefined:
                if ( event_time >= curr_ontime[0] and event_time <= curr_ontime[1] ):
                    return True

        else:
            if channel is None:
                for curr_ontime in station.onTime:
                    if ( event_time >= curr_ontime[0] and event_time <= curr_ontime[1] ):
                        return True

            else:
                for curr_channel in station.channels:
                    if channel == curr_channel.waveformID.channelCode:
                        for curr_ontime in curr_channel.onTime:
                            if ( event_time >= curr_ontime[0] and event_time <= curr_ontime[1] ):
                                return True
        
        return False

    
    def _channelOperatingDuringPeriod( self, channel, timePeriod, station=None ):
        """
        timePeriod is list of [ startTime, endTime ] or None
        startTime and endTime are of type QPDateTime

        if a PMCStation object is provided for parameter 'station', operation is checked on
        station level (station.onTime or station.onTimeRefined), channel is ignored in that case
        
        if pickInfo has been refined and new onTimes are used, station object has to be provided
        """
        
        # if function is called without valid timePeriod, don't apply (return true)
        if timePeriod is None:
            return True
        else:
            startTime = QPDateTime( timePeriod[0] )
            endTime   = QPDateTime( timePeriod[1] )

        if station is not None:

            if self.onTimesAreRefinedFlag is True:
                useOnTimes = station.onTimeRefined
            else:
                useOnTimes = station.onTime
        else:
            useOnTimes = channel.onTime
            
        # if channel/station has no start/end time, set as not operating
        if not ( ( useOnTimes is not None ) and len( useOnTimes ) > 0 ):
            return False

        # if there is overlap between the two time windows, set channel/station as operating
        for curr_ontime in useOnTimes:
            if (    ( curr_ontime[0] <= startTime and curr_ontime[1] >= endTime )
                 or ( curr_ontime[0] >= startTime and curr_ontime[1] <= endTime )
                 or (     curr_ontime[0] <= startTime and curr_ontime[1] <= endTime
                      and curr_ontime[1] >= startTime )
                 or (     curr_ontime[0] >= startTime and curr_ontime[1] >= endTime
                      and curr_ontime[0] <= endTime ) ):
                return True
            
        # no overlap found
        return False

    
    def _eventPicked( self, event, network_code, station_code ):

        # over picks for given event
        for curr_pick in event.pick:

            # is current pick from station in question?
            if (     curr_pick.waveformID.stationCode == station_code 
                 and curr_pick.waveformID.networkCode == network_code ):
                return True

        return False


    def _deletePicks( self, event, network_code, station_code ):
        """
        delete all pick and arrival information for a given event and station
        from catalog
        """

        # over picks for given event
        for curr_pick_idx in reversed( xrange( len(event.pick) ) ):
            
            # if current pick is from station in question, delete pick and arrival
            if (     event.pick[curr_pick_idx].waveformID.stationCode == station_code 
                 and event.pick[curr_pick_idx].waveformID.networkCode == network_code ):

                # get arrival for current pick, delete
                curr_ori = event.getPreferredOrigin()
                for curr_arr_idx in reversed( xrange( len(curr_ori.arrival) ) ):
                    if curr_ori.arrival[curr_arr_idx].pickID == event.pick[curr_pick_idx]:
                        del curr_ori.arrival[curr_arr_idx]

                # delete pick
                del event.pick[curr_pick_idx]


    def _stationToBeUsed( self, station ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return True


    def _calcMagnitudeFromDistance( self, distanceXYZ ):
        # this function has to be provided by the user
        # derive from PMC and fill in your code
        return 1.0
    # _calcMagnitudeFromDistance = Callable( _calcMagnitudeFromDistance )


    def _renameStation( self, station, newname ):
        station.stationCode = newname
        for ch in station.channels:
            ch.waveformID.stationCode = newname

            
    def _copyChannels( self, fromstation, tostation ):
        for ch in fromstation.channels:
            tostation.channels.append( ch )

    
    def mergeAliasStations( self ):
        
        # loop over alias list: networkCode | alias | stationCode
        for curr_alias in self.stationAliases:
        
            # look for alias in station list, can occur only once
            for curr_station_idx in reversed( xrange( len( self.stations ) ) ):

                curr_station = self.stations[curr_station_idx]
                
                # is alias matching entry in station list?
                # if not, do nothing, try next
                if (     curr_alias[0] == curr_station.networkCode
                     and curr_alias[1] == curr_station.stationCode ):

                    # gotcha: look if we have already an entry for the original station name
                    originalFound = False
                    for check_station in self.stations:
                        if (     curr_alias[0] == check_station.networkCode
                             and curr_alias[2] == check_station.stationCode ):
                            
                            # main station name is there
                            # add alias station's channels to original station
                            self._copyChannels( curr_station, check_station )
                            
                            # delete entry for alias station
                            del self.stations[curr_station_idx]
                            originalFound = True

                            print " mergeAliases: copy data from station %s to %s" % ( 
                                curr_station.stationCode, check_station.stationCode )
                            break

                    # main station name is not there
                    # rename stationCode in station and channels
                    if not originalFound:
                        self._renameStation( curr_station, curr_alias[2] )

                        print " mergeAliases: rename station from %s to %s" % ( 
                            curr_station.stationCode, curr_alias[2] )
                            
                    # found alias, end station loop
                    break

        return True


    def setDistroStyle( self, style ):

        if style != self.distStyle:
            self.distStyle = style

            for curr_sta in self.stations:
                curr_sta.distribution.restoreDistStyle( self.distStyle )
            
        return True

        
    def _getStationIdx( self, networkCode, stationCode ):
        for curr_sta_idx, curr_sta in enumerate(self.stations):
            if ( curr_sta.networkCode == networkCode and curr_sta.stationCode == stationCode ):
                return curr_sta_idx
        
        # if not yet returned, we have a new station
        self.stations.append( PMCStation( self.distStyle, smooth = self.distSmoothing ) )
        return len(self.stations)-1


    def selectStationsByNetworks( self, networks ):
        """
        use only stations that have a network code from those given in the
        networks parameter
        networks is a list of network codes

        TODO(fab): use only stations that have a subnetwork code from those given in
        the subnetworks parameter (if applicable)
        """

        # over stations in PMC
        for curr_sta_idx in reversed( xrange( len( self.stations ) ) ):
            
                # if not in list of networks, delete
                if self.stations[curr_sta_idx].networkCode not in networks:
                    del self.stations[curr_sta_idx]
                    
        return True


    def unselectStation( self, net_code, sta_code ):
        """
        unselect station with net_code, sta_code
        """

        # over stations in PMC
        for curr_sta_idx in reversed( xrange( len(self.stations) ) ):
            
            if (     ( self.stations[curr_sta_idx].networkCode == net_code )
                 and ( self.stations[curr_sta_idx].stationCode == sta_code ) ):
                del self.stations[curr_sta_idx]
                    
        return True

            
    def _preprocessStationScenario( self ):
        """
        has to be redefined in classes derived from PMC
        deselect stations from 'all' list for scenario computation
        """
        pass


    def _preprocessSubnetwork( self, subnetwork ):
        """
        use only stations that are in the given subnetwork
        """

        # over stations in PMC
        if self.subnetworks is not None:
            for curr_sta_idx in reversed( xrange( len( self.stations ) ) ):
                
                    # if not in list of networks, delete
                    sta_code = self.stations[curr_sta_idx].stationCode
                    if sta_code not in self.subnetworks[subnetwork]['stations']:
                        del self.stations[curr_sta_idx]
                    
        return True


    def _selectStationsInOperation( self, time ):
        """
        select stations by time
        time has to be mxDateTime object
        """

        curr_time = QPDateTime( time )
        
        # over stations
        for curr_sta_idx in reversed( xrange( len( self.stations ) ) ):
            
            station_operating = False

            # refined onTimes, do not loop over individual channels
            if self.onTimesAreRefinedFlag is True:

                # over time intervals
                for curr_time_interval in self.stations[curr_sta_idx].onTimeRefined:
                    if ( curr_time >= curr_time_interval[0] and curr_time <= curr_time_interval[1] ):
                        station_operating = True
                        break

            else:
                
                # over channels
                for curr_channel in self.stations[curr_sta_idx].channels:

                    # over time intervals
                    for curr_time_interval in curr_channel.onTime:

                        if ( curr_time >= curr_time_interval[0] and curr_time <= curr_time_interval[1] ):

                            # found an active channel
                            station_operating = True
                            break
                    
                    if station_operating:
                        break

            # if no operating channel found, delete station
            if not station_operating:
                del self.stations[curr_sta_idx]
                
        return True


    def _computeStationDistances( self, location ):
        """
        compute array of station distances for one location
        input: location = ( lat, lon, depth )
        output: distances[] 
        """

        distances = []

        # over stations
        for curr_sta in self.stations:

            #d, dummy = distanceBetweenPoints( location,
                                              #( curr_sta.latitude,
                                                #curr_sta.longitude,
                                                #-(curr_sta.elevation/1000.0) ),
                                              #2 )

            d = distanceBetweenPointsWGS84( location,
                ( curr_sta.latitude, curr_sta.longitude,
                    -(curr_sta.elevation/1000.0) ) )

            # print d
            distances.append( d )

        return distances


    def _computeStationProbabilities( self, magnitude, distances ):
        """
        compute probability for a given magnitude and given distance vector
        """
        probabilities = numpy.zeros( ( len(self.stations), 2 ), dtype=numpy.float )
        prob_ctr = 0

        # over stations
        for curr_sta_idx, curr_sta in enumerate( self.stations ):

            curr_prob = curr_sta.distribution.getProbability( magnitude,
                                                              distances[curr_sta_idx],
                                                              self._calcMagnitudeFromDistance )
            # print curr_sta_idx, prob_ctr, curr_prob
            if curr_prob > 0.0:
                probabilities[prob_ctr, 0] = curr_prob
                probabilities[prob_ctr, 1] = 1.0 - curr_prob
                prob_ctr = prob_ctr+1

        probabilities.resize( ( prob_ctr, 2 ) )
        return probabilities

    
    def getProbability( self, location, magnitude, time, maxStations ):
        """
        get combined probability for (location, magnitude, time)

        1) preprocess station scenario
        2) get distance vector for given location
        3) get probability vector for given magnitude and distance vector
        4) do combinatorics: get combined probability from probability vector
        """

        self._preprocessStationScenario()
        self._selectStationsInOperation( time )

        distances     = self._computeStationDistances( location )
        probabilities = self._computeStationProbabilities( magnitude, distances )

        # sort probabilities (largest first)
        # probabilities is numpy array [ prob, 1-prob ]
        # sort both columns
        indices = probabilities[:,0].argsort(axis=0)
        sortedProbabilities = probabilities[numpy.flipud(indices)]

        # if first (requiredForTriggering) probabilities are 1.0, return 1.0
        # otherwise: compute combined probability
        if floatEqual( sortedProbabilities[(self.stationsRequiredForTriggering-1), 0], 1.0 ):
            return 1.0
        else:
            combinedProbability = computeCombinedProbabilities( probabilities,
                                                                self.stationsRequiredForTriggering)
            return combinedProbability


    def _getProbabilitySet( self, location, time, magnitudes, maxStations ):
        """
        get probability array for (location, time) and a list of magnitudes
        does not compute station scenario beforehand
        
        input:  list of magnitudes
        output: list of [ sorted magnitude | combined probability ]
        """

        distances = self._computeStationDistances( location )

        combinedProbabilities = []
        probThresholdReached  = False
        
        for mag in sorted( magnitudes ):
            
            if probThresholdReached is True:
                combinedProbabilities.append( [ mag, 1.0 ] )
                continue
                
            else:
                probabilities = self._computeStationProbabilities( mag, distances )

                # if numpy array is empty (all probabilities are 0), return 0.0
                if probabilities.shape[0] < self.stationsRequiredForTriggering:
                    combinedProbabilities.append( [ mag, 0.0 ] )
                    continue
                else:
                
                    # sort probabilities (largest first)
                    # probabilities is numpy array [ prob, 1-prob ]
                    # sort both columns
                    indices = probabilities[:,0].argsort(axis=0)
                    sortedProbabilities = probabilities[numpy.flipud(indices)]

                    # if first (requiredForTriggering) probabilities are 1.0, return 1.0
                    # otherwise: compute combined probability
                    if floatEqual( sortedProbabilities[(self.stationsRequiredForTriggering-1), 0], 1.0 ):
                        probThresholdReached = True
                        combinedProbabilities.append( [ mag, 1.0 ] )
                        
                    else:
                        combinedProbability = computeCombinedProbabilities( 
                            probabilities,
                            self.stationsRequiredForTriggering )
                        if combinedProbability >= self.probThreshold:
                            # set all further probabilities to 1.0
                            probThresholdReached = True
                            combinedProbabilities.append( [ mag, 1.0 ] )
                            continue

                        else:
                            combinedProbabilities.append( [ mag, combinedProbability ] )

        return combinedProbabilities

                
    def getProbabilitySet( self, location, time, magnitudes ):
        """
        get probability array for (location, time) and a list of magnitudes
        compute station scenario beforehand
        
        input:  list of magnitudes
        output: list of [ sorted magnitude | combined probability ]
        """

        self._preprocessStationScenario()
        self._selectStationsInOperation( time )

        return self._getProbabilitySet( location, time, magnitudes )


    def getProbabilityMap( self, grid, metadata, **kwargs ):
        """
        compute PMCData on grid nodes for given sets of magnitudes and probabilities

        preprocessing steps:
        1) select stations from sub-network, if applicable
        2) compute station scenario, if applicable

        adds PMCData nodes to input QPGrid

        input:
            grid         - PMCGridN on which PMC data is computed
            metadata     - metadata for computation run

             used from metadata:
               .targetMagArray  - list of magnitudes for which probabilities are computed
               .targetProbArray - list of probabilities for which completeness magnitudes are computed
               .useDate         - mx.DateTime object, point in time for which PMC computation is made
               .currentBucket   - compute only nodes in bucket # n of m
               .bucketCount     - # of buckets
        
        kwargs: dumpfile - dump grid in ASCII format, lon | lat | probability
                           one file for each magnitude
                depth    - use this depth instead of value in depthLayer (float)
                verbose  - self explaining (bool)
        """

        self._preprocessSubnetwork( metadata.subnetwork )
        self._preprocessStationScenario()
        self._selectStationsInOperation( metadata.useDate )
        
        if ( 'dumpfile' in kwargs ) and ( kwargs['dumpfile'] is not None ):
            dumpfilename = kwargs['dumpfile']
            dumpFile = True
        else:
            dumpFile = False 

        depth = kwargs.get( 'depth', None )
        verbose = kwargs.get( 'verbose', False )

        asciimap = []
        for mag_idx in xrange( len( metadata.targetMagArray ) ):
            asciimap.append( [] )

        # over QPGrid depth layers
        for curr_depthlayer_idx in xrange( grid.depthLayer.shape[0] ):

            if depth is not None:
                curr_depth = depth
            else:
                curr_depth = grid.depthLayer[curr_depthlayer_idx, 0]

            if verbose:
                print '  next depth layer, depth: %s' % curr_depth

            # TODO(fab): select cells for current depth layer
            for curr_cell_idx in xrange( grid.cells.shape[0] ):

                curr_lon = grid.cells[curr_cell_idx, 0]
                curr_lat = grid.cells[curr_cell_idx, 1]

                curr_location = ( curr_lat, curr_lon, curr_depth )

                if verbose:
                    print '   compute cell: lat = %s, lon = %s, depth = %s' % curr_location

                # compute array of probabilities for target magnitudes
                curr_mag_prob_list = self._getProbabilitySet( curr_location,
                                                              metadata.useDate,
                                                              metadata.targetMagArray,
                                                              metadata.useMaxStationCnt )

                # loop over computed probabilities and write to pmcDataProb array
                for curr_prob_idx, curr_prob_val in enumerate( curr_mag_prob_list ):

                    grid.pmcDataProb[curr_cell_idx, curr_prob_idx] = curr_prob_val[1]

                    if dumpFile:
                        asciimap[curr_prob_idx].append( '\t'.join( ( "%6.2f" % curr_lon,
                                                                     "%6.2f" % curr_lat,
                                                                     "%6.4f" % curr_prob_val[1] ) ) )

                if metadata.targetProbArray is not None:

                    # loop over target probabilities, compute Mp value, and write to
                    # pmcDataMP array
                    for mpp_idx, mpp in enumerate( metadata.targetProbArray ):

                        grid.pmcDataMP[curr_cell_idx, mpp_idx] = self.getPMC( 
                            mpp, curr_mag_prob_list )

        if dumpFile is True:
            for mag_idx in xrange( len(asciimap) ):

                filename = "%s.%s" % ( metadata.targetMagArray[mag_idx], 
                                       dumpfilename )
                fh = writeQPData( filename, **kwargs )
                try:
                    fh.write( '\n'.join( asciimap[mag_idx] ) )
                except:
                    raise IOError, "PMC::getProbabilityMap - error writing asciimap"
                fh.close()

        return True


    def getPMC( self, targetProbability, magProbList ):
        """
        input: [ mag_array, prob_arr ], target probability
        output: magnitude of completeness
        """
        for curr_mag_prob in magProbList:
            if targetProbability <= curr_mag_prob[1]:
                return curr_mag_prob[0]
        return numpy.nan


    def getMagnitudeOfCompleteness( self ):
        # TODO(fab):
        # calls getProbabilityArray()
        #       getPMC()
        pass

    def computeDataMP( self, grid, metadata ):
        """Re-compute array of magnitudes of completeness from array of
        probabilities for each cell."""

        if metadata.targetProbArray is not None:
            for curr_depthlayer_idx in xrange( grid.depthLayer.shape[0] ):
                for curr_cell_idx in xrange( grid.cells.shape[0] ):
                    
                    # loop over target probabilities, compute Mp value, 
                    # and write to pmcDataMP array
                    for mpp_idx, mpp in enumerate( metadata.targetProbArray ):
                        grid.pmcDataMP[curr_cell_idx, mpp_idx] = self.getPMC( 
                            mpp, zip( sorted( metadata.targetMagArray ),
                                      grid.pmcDataProb[curr_cell_idx, :] ) )

        else:
            # reset array to None
            grid.pmcDataMP = None


    def setSubNetwork( self, subnet ):
        """Set trigger condition for given sub-network."""
        try:
            self.stationsRequiredForTriggering = \
                self.subnetworks[subnet]['trigger_condition']
        except KeyError:
            raise RuntimeError, \
                "PMC_NZ: %s is not a valid sub-network" % subnet
        
        self.subnetwork = subnet

