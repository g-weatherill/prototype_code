# -*- coding: utf-8 -*-
#
# quakepy/pmc/PMCGridN.py
# $Id: PMCGridN.py 335 2012-06-08 12:08:36Z tyrone $
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

PMCGridN defines a geospatial grid structure for PMC computations
PMC data is attached to grid cells
data is held in numpy arrays
"""

__version__  = '$Id: PMCGridN.py 335 2012-06-08 12:08:36Z tyrone $'
__revision__ = '$Revision: 335 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import sys
import numpy
import cStringIO

sys.path.append('../..')
sys.path.append('..')

from QPCore       import *
from QPAnnotation import *
from QPPolygon    import *

from PMCInventory import *
from PMCMetadata  import *

class PMCGridN():
    """
    QuakePy: PMCGridN

    PMCGridN provides a 3-d geospatial grid with PMC data structure for each cell
    """

    # these are the default settings
    gridParameter = { 'lonDelta'              : 0.1,
                      'latDelta'              : 0.1,
                      'lonAlign'              : 0.0,
                      'latAlign'              : 0.0,
                      'includePointOnBoundary': True,
                      'shift'                 : False }

    extent = {}
 
    def __init__( self, metadata, **kwargs ):

        self.targetMagnitudes    = metadata.targetMagArray
        self.targetProbabilities = metadata.targetProbArray
        
        self.depthLayerCols  = 1
        self.cellsCols       = 3
        self.pmcDataProbCols = len( self.targetMagnitudes )
        self.pmcDataMPCols   = len( self.targetProbabilities )
        
        # only one 'at' depthLayer at the moment
        self.depthLayer = numpy.ones( 
            (1, self.depthLayerCols), dtype=float ) * numpy.nan

        # resize those later, no. of rows is no. of cells
        # lon | lat | depthLayerRef
        self.cells = numpy.ones( 
            (1, self.cellsCols), dtype=float ) * numpy.nan 

        # prob values for target mags
        self.pmcDataProb = numpy.ones( 
            (1, self.pmcDataProbCols), dtype=float ) * numpy.nan

        # mp values for target probs
        # these are in principle obsolete, can at any time be re-computed 
        # from self.pmcDataProb
        # see method PMC.computeDataMP()
        self.pmcDataMP = numpy.ones( 
            (1, self.pmcDataMPCols), dtype=float ) * numpy.nan

    # ------------------------------------------------------------------------

    def setGridParameter( self, param ):
        """
        update self.gridParameter dictionary
        param has to be a dictionary itself
        """

        if not isinstance( param, dict ):
            raise TypeError, 'PMCGridN:setGridParameter - wrong type for parameter param, must be Python dict'
            
        for key in param.keys():
            self.gridParameter[key] = param[key]

    # ------------------------------------------------------------------------
    
    def writeXML( self, output, **kwargs ):

        if isinstance( output, basestring ):
            ostream = writeQPData( output, **kwargs )
        else:
            ostream = output
            
        if 'prettyPrint' in kwargs.keys() and kwargs['prettyPrint'] is False:
            prettyPrint = False
        else:
            prettyPrint = True
            
        if prettyPrint is True:

            # serialize to string stream
            streamSuccess = False
            try:
                curr_stream = cStringIO.StringIO()
                self.toXML( 'QPGrid', curr_stream, **kwargs )
                streamSuccess = True
            except:
                print "PMCGridN.writeXML - error in StringIO self.toXML()"

            if streamSuccess is True:
                try:
                    xmlPrettyPrint( curr_stream, ostream )
                    return
                except:
                    print "PMCGridN.writeXML - error in xmlPrettyPrint()"

        # write to output stream w/o pretty print
        # fallback if prettify has not succeeded
        try:
            self.toXML( 'QPGrid', ostream, **kwargs )
        except:
            raise IOError, "PMCGridN.writeXML - error in self.toXML()"
    
    # ------------------------------------------------------------------------
    
    def toXML( self, tagname, stream, **kwargs ):

        stream.write( '<?xml version="1.0" encoding="utf-8"?>' )
        stream.writelines( [ '<', tagname, ' xmlns="http://quakepy.org/xmlns/quakepy/1.0">' ] )
        
        if hasattr( self, 'annotation' ) and self.annotation is not None:
            self.annotation.toXML( 'annotation', stream )

        # write nodes w/ PMC data to XML
        self._writeCells2XML( 'grid', stream )
            
        # write station data
        if hasattr( self, 'stations' ) and self.stations is not None:
            stream.writelines( '<PMCStations>' )
            for sta in self.stations:
                sta.toXML( 'PMCStation', stream, **kwargs )
            stream.writelines( '</PMCStations>' )
            
        stream.writelines( [ '</', tagname, '>' ] )
        return True

    # ------------------------------------------------------------------------

    def _writeCells2XML( self, tagname, stream ):

        # open tagname
        stream.writelines( [ '<', tagname, '>' ] )

        # write depth layer
        stream.writelines( [ '<depthLayer at="', str(self.depthLayer[0,0]), '">' ] )

        # loop over cells
        for cell_idx in xrange( self.cells.shape[0] ):

            # open current cell
            stream.writelines( [ '<cell lon="', str(self.cells[cell_idx,0]),
                                 '" lat="', str(self.cells[cell_idx,1]), '">' ] )

            # open PMCData
            stream.write( '<PMCData>' )
            
            # PMC probs
            for pmc_prob_idx in xrange( self.pmcDataProbCols ):

                stream.writelines( [ '<probability magnitude="', str(self.targetMagnitudes[pmc_prob_idx]),
                                     '">', str(self.pmcDataProb[cell_idx,pmc_prob_idx]), '</probability>' ] )
            

            # PMC mp
            for pmc_mp_idx in xrange( self.pmcDataMPCols ):
                stream.writelines( [ '<mp probability="', str(self.targetProbabilities[pmc_mp_idx]),
                                     '">', str(self.pmcDataMP[cell_idx,pmc_mp_idx]), '</mp>' ] )

            # close PMCData
            stream.write( '</PMCData>' )

            # close current cell
            stream.write( '</cell>' )
        
        # close depth
        stream.write( '</depthLayer>' )
        
        # close tagname
        stream.writelines( [ '</', tagname, '>' ] )

    # ------------------------------------------------------------------------
    
    def writeNodes( self, output ):

        if isinstance( output, basestring ):
            ostream = writeQPData( output )
        else:
            ostream = output
            
        for curr_cell_idx in xrange( self.cells.shape[0] ):

            curr_lon = self.cells[curr_cell_idx, 0]
            curr_lat = self.cells[curr_cell_idx, 1]

            ostream.writelines( [ str(curr_lon), "\t",  str(curr_lat), "\n" ] )

    def writeCells( self, output ):

        if isinstance( output, basestring ):
            ostream = writeQPData( output )
        else:
            ostream = output

        for curr_cell_idx in xrange( self.cells.shape[0] ):

            curr_lon = self.cells[curr_cell_idx,0]
            curr_lat = self.cells[curr_cell_idx,1]
            
            ostream.writelines( [ str(curr_lon - 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_lat - 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ str(curr_lon - 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_lat + 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ str(curr_lon + 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_lat + 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ str(curr_lon + 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_lat - 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ str(curr_lon - 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_lat - 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ '>', "\n" ] )
    
    # ------------------------------------------------------------------------

    def setupPolygon( self, polygon, metadata ):
        """
        polygon must be of type QPPolygon
        """

        if not isinstance( polygon, QPPolygon ):
            raise TypeError, 'PMCGridN:setupPolygon - wrong type for parameter polygon, must be QPPolygon'

        if not isinstance( metadata, PMCMetadata ):
            raise TypeError, 'PMCGridN:setupPolygon - wrong type for parameter metadata, must be PMCMetadata'
        
        nodes = self._createNodes( polygon, metadata )

        # resize numpy arrays with number of cells
        self.cells.resize( (len(nodes), self.cellsCols), refcheck = 0 )

        # depth layer
        self.depthLayer[0,0] = metadata.areaDepth

        # loop over nodes of polygon
        for curr_node_idx, curr_node in enumerate(nodes):

            self.cells[curr_node_idx,0] = curr_node[0]      # lon
            self.cells[curr_node_idx,1] = curr_node[1]      # lat
            self.cells[curr_node_idx,2] = float(0)          # (depth layer reference)

        # resize PMC data arrays
        self._setupPMCData( len( nodes ) )
        
        return True

    def _setupPMCData( self, cellCtr ):

        self.pmcDataProb.resize( (cellCtr, self.pmcDataProbCols), refcheck = 0 )
        self.pmcDataMP.resize( (cellCtr, self.pmcDataMPCols), refcheck = 0 )

        return True
    
    # ------------------------------------------------------------------------
    
    def inGridCell( self, lat, lon, depth ):
        
        # depthLayer bin 0...30 km: depth  0 km is inside bin
        #                           depth 30 km is outside bin
        #
        # Note: check for depth currently not implemented
        #
        # returns ( foundLat, foundLon, foundDepthMin, foundDepthMax ) or None

        for cell_idx in xrange( self.cells.shape[0] ):

            curr_lon = self.cells[cell_idx, 0]
            curr_lat = self.cells[cell_idx, 1]
            
            latMin = curr_lat - 0.5 * self.gridParameter['latDelta']
            latMax = curr_lat + 0.5 * self.gridParameter['latDelta']
            lonMin = curr_lon - 0.5 * self.gridParameter['lonDelta']
            lonMax = curr_lon + 0.5 * self.gridParameter['lonDelta']

            if lat >= latMin and lat < latMax and lon >= lonMin and lon < lonMax:
                return ( curr_lat, curr_lon, self.depthLayer[0,0], self.depthLayer[0,0] )
         
        # depthLayer or cell not found
        return None
            
    def inGrid( self, lat, lon, depth ):
        
        if self.inGridCell( lat, lon, depth ) is None:
            return False
        else:
            return True

    # ------------------------------------------------------------------------
    
    def _createNodes( self, polygon, metadata ):
        """
        create grid nodes within polygon
        """

        nodes = []

        # check and fix (if needed) alignments, need to be smaller that deltas
        # use floating-point modulo operator
        self.gridParameter['lonAlign'] = self.gridParameter['lonAlign'] % self.gridParameter['lonDelta']
        self.gridParameter['latAlign'] = self.gridParameter['latAlign'] % self.gridParameter['latDelta']
        
        # get lat/lon range of grid node centers from alignment parameters
        nodeCenterStartLon = numpy.floor(polygon.extent['lonMin']) + self.gridParameter['lonAlign']
        nodeCenterEndLon   = numpy.ceil(polygon.extent['lonMax'])

        nodeCenterStartLat = numpy.floor(polygon.extent['latMin']) + self.gridParameter['latAlign']
        nodeCenterEndLat   = numpy.ceil(polygon.extent['latMax'])

        print " creating nodes from lon: %s, %s and lat: %s, %s" % ( nodeCenterStartLon, nodeCenterEndLon,
                                                                     nodeCenterStartLat, nodeCenterEndLat )
        # loop over grid of node centers
        for fLon in frange( nodeCenterStartLon, nodeCenterEndLon, self.gridParameter['lonDelta'] ):
            for fLat in frange( nodeCenterStartLat, nodeCenterEndLat, self.gridParameter['latDelta'] ):

                if self.gridParameter['includePointOnBoundary'] is True:

                    # check if node center is in polygon or on boundary
                    if polygon.isInsideOrOnBoundary( fLon, fLat ):
                        
                        # Create nodes
                        nodes.append( [ fLon, fLat ] )
                else:
                    
                    # check if node center is truly inside polygon (not on boundary)
                    if polygon.isInside( fLon, fLat ):

                        # Create nodes
                        nodes.append( [ fLon, fLat ] )

        # just create the requested bnucket:
        print '  computing bucket # %d of %d buckets' % ( metadata.currentBucket, metadata.bucketCount )
        allNodeCount = len(nodes)
        minCellIdx = int(math.floor((metadata.currentBucket - 1) * (float(allNodeCount)/metadata.bucketCount)))
        maxCellIdx = int(math.floor((metadata.currentBucket) * (float(allNodeCount)/metadata.bucketCount)))
        print '  computing nodes [%d, %d) of %d nodes' % ( minCellIdx, maxCellIdx, allNodeCount )
        dummy = []
        for idx in xrange( minCellIdx, maxCellIdx ):
            dummy.append(nodes[idx])
        nodes = dummy
        dummy = None
        print '  all-node-count: %d, this-bucket-node-count: %d' % ( allNodeCount, len(nodes))
        
        if self.gridParameter['shift'] is True:

            for nCnt in xrange( len( nodes ) ):
                if nodes[nCnt][0] > 180.0:
                    nodes[nCnt][0] = nodes[nCnt][0] - 360.0

        self._getExtent( nodes )

        return nodes

    def _getExtent( self, nodes ):
        
        # Determine extent of cell grid
        lons = [ v[0] for v in nodes ]
        lats = [ v[1] for v in nodes ]

        self.extent = { 'lonMin': min( lons ) - 0.5 * self.gridParameter['lonDelta'],
                        'lonMax': max( lons ) + 0.5 * self.gridParameter['lonDelta'],
                        'latMin': min( lats ) - 0.5 * self.gridParameter['latDelta'],
                        'latMax': max( lats ) + 0.5 * self.gridParameter['latDelta'] }
        