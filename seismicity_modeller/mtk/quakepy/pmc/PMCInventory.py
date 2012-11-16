# -*- coding: iso-8859-1 -*-
#
# quakepy/pmc/PMCInventory.py
# $Id: PMCInventory.py 248 2009-11-04 12:28:15Z fab $
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

__version__  = '$Id: PMCInventory.py 248 2009-11-04 12:28:15Z fab $'
__revision__ = '$Revision: 248 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import copy
import cPickle
import math
import numpy

import matplotlib
matplotlib.use('PS')
from pylab import *

from mx.DateTime     import Date
from mx.DateTime.ISO import ParseDateTimeUTC

import pyRXP
from   xml.sax   import saxutils

sys.path.append('..')

from QPCore        import *
from PMCFunctions  import *

from WaveformStreamID import WaveformStreamID


POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)


class PMCStation( QPObject ):
    
    def __init__( self, dist_style = None, **kwargs ):

        # holds objects of type 'PMCChannel'
        self.channels = []

        # onTime per station = time intervals in which at least one channel was operating
        self.onTime = []

        # PMCProbabilityDistribution object, holds pickInfo and distro
        self.distribution = PMCProbabilityDistribution( dist_style, **kwargs )

        # start and end date of catalog used for pickInfo
        self.catalogStartDate      = None   # type QPDateTime
        self.catalogEndDate        = None   # type QPDateTime
        
        # inforamtion for refined picks
        self.onTimesAreRefinedFlag = False
        
        self.onTimeRefined         = None
        self.interPickTime         = None   # list, times between seen arrivals at station, in days 
        
        self.meanInterPickTime     = None   # mean time between seen arrivals at station, in days
        self.timeIntervalScale     = None

        self.ontimeMinPickCtr      = None
        
    ## -----------------------------------------------------------------------
    
    def fromXML( self, tree, additionalElements = None ):

        # XML attributes: networkCode, stationCode, locationCode, latitude, longitude, elevation (mandatory)
        attr_dict = tree[POS_ATTRS]
        self.networkCode = str( attr_dict['networkCode'] )
        self.stationCode = str( attr_dict['stationCode'] )

        if 'locationCode' in attr_dict.keys():
            self.locationCode = str( attr_dict['locationCode'] )
            
        self.latitude    = float( attr_dict['latitude'] )
        self.longitude   = float( attr_dict['longitude'] )
        self.elevation   = float( attr_dict['elevation'] )
        
        for child in tree[POS_CHILDREN]:
            
            # multiple children of same type
            if child[POS_TAGNAME] == 'PMCChannel':
                cha = PMCChannel()
                cha.fromXML( child, additionalElements )
                self.channels.append( cha )

            elif child[POS_TAGNAME] == 'onTime':

                # onTime element has 2 required attributes, 'start' and 'end'
                # convert to QPDateTime objects
                attr_dict = child[POS_ATTRS]
                self.onTime.append( [ QPDateTime( attr_dict['start'] ),
                                      QPDateTime( attr_dict['end'] ) ] )
                                             
            elif child[POS_TAGNAME] == 'onTimeRefined':

                if self.onTimeRefined is None:
                    self.onTimeRefined = []
                    self.onTimesAreRefinedFlag = True

                attr_dict = child[POS_ATTRS]
                self.onTimeRefined.append( [ QPDateTime( attr_dict['start'] ),
                                             QPDateTime( attr_dict['end'] ) ] )

            elif child[POS_TAGNAME] == 'interPickTime':

                inter_pick_times_list = map( float, child[POS_CHILDREN].pop().split() )
                if len( inter_pick_times_list ) > 0:
                    self.interPickTime = inter_pick_times_list

            # children of complex type
            elif child[POS_TAGNAME] == 'PMCProbabilityDistribution':
                self.distribution = PMCProbabilityDistribution()
                self.distribution.fromXML( child, additionalElements )

            # children of basic type
            elif child[POS_TAGNAME] == 'description':

                # convert description to Unicode
                description_unicode   = unicode( child[POS_CHILDREN].pop(), 'utf-8' )
                description_sanitized = saxutils.unescape( description_unicode )
                #description_sanitized = description_unicode

                # ----- this works, for non-unicode self.description
                #description_sanitized = saxutils.unescape( child[POS_CHILDREN].pop() )

                self.description = description_sanitized

            elif child[POS_TAGNAME] == 'catalogStartDate':

                self.catalogStartDate = QPDateTime( child[POS_CHILDREN].pop() )

            elif child[POS_TAGNAME] == 'catalogEndDate':

                self.catalogEndDate = QPDateTime( child[POS_CHILDREN].pop() )
                    
            elif child[POS_TAGNAME] == 'meanInterPickTime':

                self.meanInterPickTime = float( child[POS_CHILDREN].pop() )

            elif child[POS_TAGNAME] == 'timeIntervalScale':

                self.timeIntervalScale = float( child[POS_CHILDREN].pop() )

            elif child[POS_TAGNAME] == 'ontimeMinPickCtr':

                self.ontimeMinPickCtr = int( child[POS_CHILDREN].pop() )
                

    def toXML( self, tagname, stream, **kwargs ):
        """
        kwargs: noStationChildElements=True: write only <PMCStation> tag with attributes (default: False) 
        """
        
        stream.writelines( [ '<', tagname ] )

        # XML attributes: networkCode, stationCode, latitude, longitude, elevation (mandatory)
        stream.write( ''.join( (' networkCode="', self.networkCode, '"') ) )
        stream.write( ''.join( (' stationCode="', self.stationCode, '"') ) )

        if hasattr( self, 'locationCode' ) and self.locationCode is not None:
            stream.write( ''.join( (' locationCode="', self.locationCode, '"') ) )

        stream.write( ''.join( (' latitude="', str(self.latitude), '"') ) )
        stream.write( ''.join( (' longitude="', str(self.longitude), '"') ) )
        stream.write( ''.join( (' elevation="', str(self.elevation), '"') ) )

        stream.write( '>' )

        if not ( 'noStationChildElements' in kwargs.keys() and kwargs['noStationChildElements'] is True ):
            
            # children of basic type
            if hasattr( self, 'description' ) and self.description is not None:

                # sanitize description string for XML

                # self.description is Unicode!
                description_str_xml = saxutils.escape( self.description ).encode( 'utf-8', 'xmlcharrefreplace' )

                # ----- this works, for non-unicode self.description
                #description_str_xml = saxutils.escape( self.description )

                stream.write( ''.join( ('<description>', description_str_xml, '</description>') ) )

            if hasattr( self, 'catalogStartDate' ) and self.catalogStartDate is not None:

                stream.write( ''.join( ('<catalogStartDate>', self.catalogStartDate.toISO( showsecfrac=False ), '</catalogStartDate>') ) )
                
            if hasattr( self, 'catalogEndDate' ) and self.catalogEndDate is not None:

                stream.write( ''.join( ('<catalogEndDate>', self.catalogEndDate.toISO( showsecfrac=False ), '</catalogEndDate>') ) )
                
            if hasattr( self, 'meanInterPickTime' ) and self.meanInterPickTime is not None:

                stream.write( ''.join( ('<meanInterPickTime>', str( self.meanInterPickTime ), '</meanInterPickTime>') ) )
                
            if hasattr( self, 'timeIntervalScale' ) and self.timeIntervalScale is not None:

                stream.write( ''.join( ('<timeIntervalScale>', str( self.timeIntervalScale ), '</timeIntervalScale>') ) )

            if hasattr( self, 'ontimeMinPickCtr' ) and self.ontimeMinPickCtr is not None:

                stream.write( ''.join( ('<ontimeMinPickCtr>', str( self.ontimeMinPickCtr ), '</ontimeMinPickCtr>') ) )

            if hasattr( self, 'interPickTime' ) and self.interPickTime is not None and len( self.interPickTime ) > 0:

                stream.write( ''.join( ('<interPickTime>', ' '.join( map( str, self.interPickTime ) ), '</interPickTime>') ) )
                
            ## multiple children of same type

            # channels
            for cha in self.channels:
                cha.toXML( 'PMCChannel', stream )
                    
            # onTime
            for ot in self.onTime:
                stream.writelines( ''.join( ( '<onTime start="', ot[0].toISO( showsecfrac=False ),
                                            '" end="', ot[1].toISO( showsecfrac=False ), '"/>' ) ) )
                                            
            # onTimeRefined (optional)
            if hasattr( self, 'onTimeRefined' ) and self.onTimeRefined is not None:
                for ot in self.onTimeRefined:
                    stream.writelines( ''.join( ( '<onTimeRefined start="', ot[0].toISO( showsecfrac=False ),
                                                '" end="', ot[1].toISO( showsecfrac=False ), '"/>' ) ) )

            # children of complex type
            if hasattr( self, 'distribution' ) and self.distribution is not None:
                self.distribution.toXML( 'PMCProbabilityDistribution', stream )

        stream.writelines( [ '</', tagname, '>' ] )
        return True

    ## -----------------------------------------------------------------------
         
    def fillup( self, calcMagnitudeFromDistance, override=False ):
        self.distribution.fillup( calcMagnitudeFromDistance, override )


class PMCStationList( QPObject ):

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'station', 'PMCStation', 'element', PMCStation, 'multiple' ),
                  ) )
    # <!-- UML2Py end -->

    def __init__( self, **kwargs ):
        super( PMCStationList, self ).__init__( None, **kwargs )
        self.elements.extend( self.addElements )
                        
        self._initMultipleElements()

## ---------------------------------------------------------------------------
    
class PMCChannel( QPObject ):
    
    def __init__( self ):

        # onTime holds a list of [ startTime, endTime ] pairs
        # both times are of type QPDateTime
        self.onTime     = []
        self.waveformID = None


    def fromXML( self, tree, additionalElements = None ):

        for child in tree[POS_CHILDREN]:

            # children of complex type
            if child[POS_TAGNAME] == 'waveformID':
                self.waveformID = WaveformStreamID()
                self.waveformID.fromXML( child, additionalElements )
                
            # multiple children of same type
            elif child[POS_TAGNAME] == 'onTime':

                # onTime element has 2 required attributes
                # convert to QPDateTime objects
                attr_dict = child[POS_ATTRS]
                self.onTime.append( [ QPDateTime( attr_dict['start'] ),
                                      QPDateTime( attr_dict['end'] ) ] )


    def toXML( self, tagname, stream  ):

        stream.writelines( [ '<', tagname, '>' ] )

        # children of complex type
        if hasattr( self, 'waveformID' ) and self.waveformID is not None:
            self.waveformID.toXML( 'waveformID', stream )
            
        # multiple children of same type
        if hasattr( self, 'onTime' ) and self.onTime is not None:
            for ot in self.onTime:
                stream.writelines( ''.join( ( '<onTime start="', ot[0].toISO( showsecfrac=False ),
                                              '" end="', ot[1].toISO( showsecfrac=False ), '"/>' ) ) )
                                               
        stream.writelines( [ '</', tagname, '>' ] )
        return True

## ---------------------------------------------------------------------------

class PMCProbabilityDistribution( QPObject ):
    """
    QuakePy: PMCProbabilityDistribution

    style:  -1 free
             0 don't compute distribution, use raw data
             1 M=0:0.2:4; D=5:5:200;               # 840 (low-res, quick inspection)
             2 M=0:0.05:4; D=1:1:200;              # 16200 ("classical" resolution)
             3 M=0:0.1:5; D=2:2:600;               # 15300 (medium resolution with large range, Italy & SCSN)
             4 M=0:0.2:4; D=10:10:600;             # 1260 (quick inspection, Japan)
             5 M=-1.0:0.1:5.0; D=5:5:1000;               (Japan)
    """

    distro_binning = ( ( ( 0.0,  0.2, 4.0 ), ( 5.0, 5.0, 200.0 ) ),
                       ( ( 0.0, 0.05, 4.0 ), ( 1.0, 1.0, 200.0 ) ),
                       ( ( 0.0,  0.1, 5.0 ), ( 2.0, 2.0, 600.0 ) ),
                       ( ( 0.0,  0.2, 4.0 ), ( 10.0, 10.0, 600.0 ) ),
                       ( ( -1.0,  0.1, 5.0 ), ( 5.0, 5.0, 1000.0 ) ) )

    pickInfoCols = 6         # no. of columns in pickInfo
    
    def __init__( self, dist_style = None, **kwargs ):
        """
        distro_binning: ( magnitude=( start, step, end ), distance=( start, step, end ) )

        kwargs: smooth = True | False
        """

        if 'smooth' in kwargs.keys() and kwargs['smooth'] is False:
            self.smooth = False
        else:
            # default: smooth distribution
            self.smooth = True

        # init pickInfo (NumPy array), filled in PMC::assignPicks
        self.pickInfo  = None
        self.distro    = None
        self.distStyle = None

        self.restoreDistStyle( dist_style )

    ## -----------------------------------------------------------------------
    
    def fromXML( self, tree, additionalElements=None ):

        # XML attributes: distStyle, smooth (mandatory)
        attr_dict      = tree[POS_ATTRS]
        self.distStyle = int( attr_dict['distStyle'] )
        self.smooth    = bool( attr_dict['smooth'] )
        
        for child in tree[POS_CHILDREN]:
            
            # children of 'basic' type
            if child[POS_TAGNAME] == 'magnitudeGrid':
                
                value_arr = map( float, child[POS_CHILDREN].pop().split() )
                self.magnitudeGrid = numpy.array( value_arr )
                
            elif child[POS_TAGNAME] == 'distanceGrid':
                
                value_arr = map( float, child[POS_CHILDREN].pop().split() )
                self.distanceGrid = numpy.array( value_arr )

            elif child[POS_TAGNAME] == 'pickInfo' and len( child[POS_CHILDREN] ) > 0:

                child = child[POS_CHILDREN].pop()
                
                # exclude 'whitespace only' children
                if len( child.strip() ) > 0:
                    
                    value_arr = map( float, child.split() )
                    self.pickInfo = numpy.array( value_arr )

            elif child[POS_TAGNAME] == 'distro' and len( child[POS_CHILDREN] ) > 0:

                child = child[POS_CHILDREN].pop()

                # exclude 'whitespace only' children
                if len( child.strip() ) > 0:

                    value_arr = map( float, child.split() )
                    self.distro = numpy.array( value_arr )

        # reshape pickInfo and distro
        if self.pickInfo is not None:
            ( pickInfoLines, remainder ) = divmod( len(self.pickInfo), self.pickInfoCols )

            if remainder == 0:
                self.pickInfo = numpy.reshape( self.pickInfo, ( pickInfoLines, self.pickInfoCols ) )
            else:
                raise ValueError, "pickInfo dimensions are wrong"

        if self.distro is not None:

            if ( self.magnitudeGrid is not None ) and ( self.distanceGrid is not None ):
                self.distro = numpy.reshape( self.distro, ( len(self.distanceGrid), len(self.magnitudeGrid) ) )
            else:
                raise ValueError, "PMCInventory::fromXML - no distance and magnitude grids"


    def toXML( self, tagname, stream ):
        stream.writelines( [ '<', tagname ] )

        # XML attributes: distStyle, smooth (mandatory)
        stream.write( ''.join( (' distStyle="', str(self.distStyle), '"') ) )
        stream.write( ''.join( (' smooth="', str(self.smooth), '"') ) )
        
        stream.write( '>' )

        # children of 'basic' type
        if hasattr( self, 'magnitudeGrid' ) and self.magnitudeGrid is not None:

            # could possibly also use numpy.array2string()
            value_str = ' '.join( map( str, self.magnitudeGrid.tolist() ) )
            stream.write( ''.join( ( '<magnitudeGrid>', value_str, '</magnitudeGrid>' ) ) )

        if hasattr( self, 'distanceGrid' ) and self.distanceGrid is not None:

            value_str = ' '.join( map( str, self.distanceGrid.tolist() ) )
            stream.write( ''.join( ( '<distanceGrid>', value_str, '</distanceGrid>' ) ) )

        if hasattr( self, 'pickInfo' ) and self.pickInfo is not None:

            value_str = ' '.join( map( str, self.pickInfo.ravel().tolist() ) )

            # exclude 'whitespace only' value_str
            if len( value_str.strip() ) > 0:
                stream.write( ''.join( ( '<pickInfo>', value_str, '</pickInfo>' ) ) )

        if hasattr( self, 'distro' ) and self.distro is not None:

            value_str = ' '.join( map( str, self.distro.ravel().tolist() ) )

            # exclude 'whitespace only' value_str
            if len( value_str.strip() ) > 0:
                stream.write( ''.join( ( '<distro>', value_str, '</distro>' ) ) )

        stream.writelines( [ '</', tagname, '>' ] )
        return True

    ## -----------------------------------------------------------------------

    def getProbability( self, magnitude, distance, calcMagnitudeFromDistance, override=False ):
        """
        get probability for given magnitude and distance
        if not yet available in grid, compute probability value
        """

        # print magnitude, distance
        if (    ( magnitude < self.magnitudeGrid[0] )
             or ( distance > self.distanceGrid[-1] ) ):
            return 0.0

        if magnitude > self.magnitudeGrid[-1]:
            magnitude = self.magnitudeGrid[-1]

        if distance < self.distanceGrid[0]:
            distance = self.distanceGrid[0]
            
        # get grid indices for (magnitude, distance)
        ( mag_idx, dist_idx ) = self._getMagDistIdx( magnitude, distance )

        #print " mag %s, dist %s -> idx %s, %s" % ( magnitude, distance, mag_idx, dist_idx )
        
        # if prob. for  (magnitude, distance) is NaN, compute it
        if ( numpy.isnan(self.distro[dist_idx, mag_idx]) or ( override is True ) ):

            if self.smooth is True:
                # fill up probability grid with smoothed probability value
                self.distro[dist_idx, mag_idx] = self.calcSmoothedProbability( calcMagnitudeFromDistance,
                                                                               magnitude,
                                                                               distance,
                                                                               mag_idx,
                                                                               dist_idx )
            else:
                # fill up probability grid with raw probability value (not smoothed)
                self.distro[dist_idx, mag_idx], dummy = self.calcRawProbability( 
                                                                            magnitude,
                                                                            distance,
                                                                            calcMagnitudeFromDistance )

        # return prob. value from distribution Grid
        return self.distro[dist_idx, mag_idx]

        
    def _getMagDistIdx( self, magnitude, distance ):
        """
        return tuple of located indices ( mag_idx, dist_idx )
        """

        ## locating magnitude and distance in grid: we could
        ## use QPUtils.locateInArray(), but we compute the index directly
        ## since it is faster
        
        # mag_idx = locateInArray( self.magnitudeGrid, magnitude )
        # dist_idx = locateInArray( self.distanceGrid, distance )

        if self.distStyle == 1:
            mag_idx  = int(round(magnitude * 5))
            dist_idx = int(round(distance/5)) - 1
        elif self.distStyle == 2:
            mag_idx  = int(round(magnitude * 20))
            dist_idx = int(round(distance)) - 1
        elif self.distStyle == 3:
            mag_idx  = int(round(magnitude * 10))
            dist_idx = int(round(distance/2)) - 1
        elif self.distStyle == 4:
            mag_idx  = int(round(magnitude * 5))
            dist_idx = int(round(distance/10)) - 1
        elif self.distStyle == 5:
            mag_idx  = int(round(magnitude * 10)) + 10
            dist_idx = int(round(distance/5)) - 1
                
        if dist_idx == -1:
            dist_idx = 0

        return ( mag_idx, dist_idx )

    
    def calcSmoothedProbability( self, calcMagnitudeFromDistance, magnitude, distance,
                                 mag_idx=None, dist_idx=None ):
        """
        calc smoothed probability for given magnitude and distance
        """

        # if indices for magnitude and distance have not been passed, compute again
        if ( not mag_idx ) or ( not dist_idx ):
            ( mag_idx, dist_idx ) = self._getMagDistIdx( magnitude, distance )
            
        # iterate over magnitudes from smallest magnitude up to mag. in question
        for curr_mag_idx in xrange( mag_idx+1 ):
            
            # iterate over distances from largest distance up to dist. in question
            for curr_dist_idx in reversed( xrange( dist_idx, len(self.distanceGrid) ) ):

                if numpy.isnan( self.distro[curr_dist_idx, curr_mag_idx] ):
                    self.distro[curr_dist_idx, curr_mag_idx], dummy = self.calcRawProbability(
                                                      self.magnitudeGrid[curr_mag_idx],
                                                      self.distanceGrid[curr_dist_idx],
                                                      calcMagnitudeFromDistance )

                    if ( curr_dist_idx == len(self.distanceGrid)-1 ) and ( curr_mag_idx == 0 ):

                        # left upper corner, no smoothing required
                        continue
                
                    elif curr_mag_idx == 0:
                        # left edge, compare only with upper neighbor (curr_dist_idx+1)
                        self.distro[curr_dist_idx, curr_mag_idx] = max( self.distro[curr_dist_idx, curr_mag_idx],
                                                                        self.distro[curr_dist_idx+1, curr_mag_idx] )
                    elif curr_dist_idx == len(self.distanceGrid)-1:
                        # upper edge, compare only with left neighbor (curr_mag_idx-1)
                        self.distro[curr_dist_idx, curr_mag_idx] = max( self.distro[curr_dist_idx, curr_mag_idx],
                                                                        self.distro[curr_dist_idx, curr_mag_idx-1] )
                    else:

                        # compare with upper (curr_dist_idx+1) and left (curr_mag_idx-1) neighbor, set cell to maximum
                        self.distro[curr_dist_idx, curr_mag_idx] = max( self.distro[curr_dist_idx, curr_mag_idx],
                                                                        self.distro[curr_dist_idx+1, curr_mag_idx],
                                                                        self.distro[curr_dist_idx, curr_mag_idx-1] )

        return self.distro[dist_idx, mag_idx]

    
    def calcRawProbability( self, magnitude, distance, calcMagnitudeFromDistance ):
        """
        calc raw probability for given magnitude and distance
        """

        # if pickInfo has no rows, return probability 0.0
        if self.pickInfo.shape[0] == 0:
            return ( 0.0, 0 )
        
        # create numpy array with 4 columns and rows as in pickInfo
        values = numpy.zeros( (self.pickInfo.shape[0], 4), dtype=float )

        # 1st column: picked or not picked?
        values[:,0] = self.pickInfo[:,2]
        
        # 2nd column: magnitude differences
        values[:,1] = self.pickInfo[:,1] - magnitude

        # 3rd column: translate difference in distances into a magnitude difference
        logDistance = calcMagnitudeFromDistance( distance )
        values[:,2] = self.pickInfo[:,4] - logDistance

        # 4rd column: calculate total difference in "magnitude units"
        values[:,3] = numpy.sqrt( values[:,1]**2 + values[:,2]**2 )

        # select only picks with a "magnitude difference" of less than 0.1

        sel = ( values[:,3] < 0.1 )             # selector: TODO use magnitude threshold
            
        # check if selector contains at least one 'True'
        if sel.sum() != 0.0:
            sample = values[sel.T,:]
            numSample = sample.shape[0]
        else:
            numSample   = 0

        if ( numSample < 10 ):                                                       # not enough samples available
            nonSample = values[numpy.logical_not(sel),:]                             # invert selector

            sel = numpy.logical_and( (nonSample[:,1] <= 0), (nonSample[:,2] >= 0) )  # select points w/ larger distance, smaller magnitude
            selection = nonSample[sel.T,:]

            indices = selection[:,3].argsort(axis=0)                                 # sort for distance
            if indices.shape[0] < (10 - numSample):                                  # enough points to fill array up to 10?
                
                probability = 0.0                                                    # no, give up (return zero)
            else:
                sortedSelection = selection[indices]

                if numSample > 0:
                    newSample = numpy.concatenate( (sample, sortedSelection[0:(10-numSample),:]), axis=0 )
                else:
                    newSample = sortedSelection[0:10,:]
                    
                probability = newSample[:,0].sum() / 10.0
                numSample = 10
        else:
            probability = sample[:,0].sum() / numSample

        return ( probability, numSample )

        
    def fillup( self, calcMagnitudeFromDistance, override=False ):

        # if self.pickInfo has only unpicked entries, set whole distro to zeros
        if self.pickInfo[:,2].sum() == 0.0:
            self.distro = numpy.zeros( ( len(self.distanceGrid), len(self.magnitudeGrid) ), dtype=float )
        else:
            
            for dist in self.distanceGrid:
                for mag in self.magnitudeGrid:

                    # print "compute prob: dist %s of %s, %s; mag %s of %s, %s" % ( dist, self.distanceGrid[0], self.distanceGrid[-1], mag, self.magnitudeGrid[0], self.magnitudeGrid[-1] )
                    prob = self.getProbability( mag, dist, calcMagnitudeFromDistance, override )

        return True


    def setSmooth( self, smooth ):
        if smooth != self.smooth:
            if self.distStyle != 0:
                self.resetDistro()
            self.smooth = smooth


    def resetDistro( self ):
        if self.distStyle != 0:
            self.distro = numpy.ones( ( len(self.distanceGrid), len(self.magnitudeGrid) ), dtype=float ) * numpy.nan


    def restoreDistStyle( self, dist_style ):

        if dist_style is None:
            self.distStyle = 0

        elif dist_style != self.distStyle:

            self.distStyle = dist_style

            if self.distStyle > 0:
                self.magnitudeGrid = frange( self.distro_binning[self.distStyle-1][0][0],
                                             self.distro_binning[self.distStyle-1][0][2],
                                             self.distro_binning[self.distStyle-1][0][1] )
                self.distanceGrid  = frange( self.distro_binning[self.distStyle-1][1][0],
                                             self.distro_binning[self.distStyle-1][1][2],
                                             self.distro_binning[self.distStyle-1][1][1] )

                self.distro = numpy.ones( ( len(self.distanceGrid), len(self.magnitudeGrid) ), dtype=float ) * numpy.nan
                                              
            elif self.distStyle == -1:
                # get 'free' distribution
                pass
            elif self.distStyle == 0:
                # do nothing
                pass
            else:
                raise ValueError, "PMCProbabilityDistribution::restoreDistStyle - invalid style (%s)" % ( self.distStyle )

    ## -----------------------------------------------------------------------
    
    def dump( self, filename ):
        outstring_arr = []

        for curr_mag_idx in xrange( len(self.magnitudeGrid ) ):
            for curr_dist_idx in xrange( len(self.distanceGrid) ):
            
                line_arr = ( "%6.2f" % self.magnitudeGrid[curr_mag_idx],
                             "%6.1f" % self.distanceGrid[curr_dist_idx],
                             "%6.4f" % self.distro[curr_dist_idx, curr_mag_idx] )
                        
                outstring_arr.extend( ( '\t'.join( line_arr ), '\n' ) )
           
        outstring = ''.join( outstring_arr )

        fh = writeQPData( filename )
        try:
            fh.write( outstring )
        except:
            raise IOError, "dump - write error"
        fh.close()

        return True
        
        