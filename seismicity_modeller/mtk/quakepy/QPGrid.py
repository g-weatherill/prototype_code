# -*- coding: iso-8859-1 -*-
#
# quakepy/QPGrid.py
# $Id: QPGrid.py 208 2009-06-19 14:53:10Z fab $
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

__version__  = '$Id: QPGrid.py 208 2009-06-19 14:53:10Z fab $'
__revision__ = '$Revision: 208 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import cStringIO
import pyRXP

# internal includes
from QPCore       import *
from QPAnnotation import *
from QPPolygon    import *

POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)

# -------------------------------------------------------------------------------------

class Cell( QPObject ):
    """
    QuakePy: Cell
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'lat', 'lat', 'attribute', float, 'basic' ),
                    QPElement( 'lon', 'lon', 'attribute', float, 'basic' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, **kwargs ):
        super( Cell, self ).__init__( **kwargs )
        
        self.elements.extend( self.addElements )
        self._initMultipleElements()

# -------------------------------------------------------------------------------------

class DefaultCellDimension( QPObject ):
    """
    QuakePy: DefaultCellDimension
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'latRange', 'latRange', 'attribute', float, 'basic' ),
                    QPElement( 'lonRange', 'lonRange', 'attribute', float, 'basic' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, **kwargs ):
        super( DefaultCellDimension, self ).__init__( **kwargs )
        
        self.elements.extend( self.addElements )
        self._initMultipleElements()

# -------------------------------------------------------------------------------------

class DepthLayer( QPObject ):
    """
    QuakePy: DepthLayer
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'min', 'min', 'attribute', float, 'basic' ),
                    QPElement( 'max', 'max', 'attribute', float, 'basic' ),
                    QPElement( 'at', 'at', 'attribute', float, 'basic' ),
                    QPElement( 'cell', 'cell', 'element', Cell, 'multiple' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, **kwargs ):
        super( DepthLayer, self ).__init__( **kwargs )
        
        self.elements.extend( self.addElements )
        self._initMultipleElements()

# -------------------------------------------------------------------------------------

class Grid( QPObject ):
    """
    QuakePy: Grid
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'defaultCellDimension', 'defaultCellDimension', 'element', DefaultCellDimension, 'complex' ),
                    QPElement( 'depthLayer', 'depthLayer', 'element', DepthLayer, 'multiple' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, **kwargs ):
        super( Grid, self ).__init__( **kwargs )
        
        self.elements.extend( self.addElements )
        self._initMultipleElements()
        
# ----------------------------------------------------------------------------

class QPGrid( QPObject ):

    gridParameter = { 'lonDelta'              : 0.1,
                      'latDelta'              : 0.1,
                      'lonAlign'              : 0.0,
                      'latAlign'              : 0.0,
                      'includePointOnBoundary': True,
                      'shift'                 : True  }

    def __init__( self, input = None, **kwargs ):
        super( QPGrid, self ).__init__( **kwargs )

        # set element axis
        self.setElementAxis( '/QPGrid' )

        if input is not None:
            if isinstance( input, basestring ):
                istream = getQPDataSource( input )
            else:
                istream = input
            
            self.readXML( istream )

    # ------------------------------------------------------------------------
    
    def readXML( self, input, **kwargs  ):

        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs  )
        else:
            istream = input
            
        # get whole content of stream at once
        lines = istream.read()
        
        # check if it is XML
        if not lines.startswith('<?xml'):
            raise IOError, 'QPGrid::readXML - input stream is not XML'
            
        tree = pyRXP.Parser().parse( lines )
        
        if tree[POS_TAGNAME] != 'QPGrid':
            raise TypeError, 'QPGrid::readXML - input stream is not of QPGrid type'
        
        # NOTE: annotation element is not read from XML
        for child in tree[POS_CHILDREN]:
            
            # possible child elements: grid
            if child[POS_TAGNAME] == 'grid':
              
                # QPGrid can only contain one single grid element
                # check if grid already existing
                if hasattr( self, 'grid' ):
                    raise TypeError, 'QPGrid::readXML - only single occurrence of grid element allowed'
                
                self.grid = Grid( parentAxis = self.elementAxis,
                                  elementName = 'grid' )
                                                
                self.grid.fromXML( child, self.elements )

        if not hasattr( self, 'grid' ):
            raise TypeError, 'QPGrid::readXML - no grid tag found'

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
            try:
                curr_stream = cStringIO.StringIO()
                self.toXML( 'QPGrid', curr_stream )
                streamSuccess = True
            except:
                streamSuccess = False
                print "QPGrid::writeXML - error in StringIO self.toXML()"

            if streamSuccess is True:
                try:
                    xmlPrettyPrint( curr_stream, ostream )
                    return
                except:
                    print "QPGrid::writeXML - error in xmlPrettyPrint()"

        # write to output stream w/o pretty print
        # fallback if prettify has not succeeded
        try:
            self.toXML( 'QPGrid', ostream )
        except:
            raise IOError, "QPGrid::writeXML - error in self.toXML()"
    
    # ------------------------------------------------------------------------
    
    def toXML( self, tagname, stream ):

        stream.write( '<?xml version="1.0" encoding="utf-8"?>' )
        stream.writelines( [ '<', tagname, ' xmlns="http://quakepy.org/xmlns/quakepy/1.0">' ] )
        
        if hasattr( self, 'annotation' ) and self.annotation is not None:
            self.annotation.toXML( 'annotation', stream )
            
        if hasattr( self, 'grid' ) and self.grid is not None:
            self.grid.toXML( 'grid', stream )
            
        # NOTE: this is not a typical element of QPGrid
        # should later be placed somewhere else
        if hasattr( self, 'stations' ) and self.stations is not None:
            stream.writelines( '<PMCStations>' )
            for sta in self.stations:
                sta.toXML( 'PMCStation', stream )
            stream.writelines( '</PMCStations>' )
            
        stream.writelines( [ '</', tagname, '>' ] )
        return True

    # ------------------------------------------------------------------------

    def writeNodes( self, output ):

        if isinstance( output, basestring ):
            ostream = writeQPData( output )
        else:
            ostream = output
            
        for curr_cell in self.grid.depthLayer[0].cell:
            ostream.writelines( [ str(curr_cell.lon), "\t",  str(curr_cell.lat), "\n" ] )

    def writeCells( self, output ):

        if isinstance( output, basestring ):
            ostream = writeQPData( output )
        else:
            ostream = output

        for curr_cell in self.grid.depthLayer[0].cell:
            
            ostream.writelines( [ str(curr_cell.lon - 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_cell.lat - 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ str(curr_cell.lon - 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_cell.lat + 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ str(curr_cell.lon + 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_cell.lat + 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ str(curr_cell.lon + 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_cell.lat - 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ str(curr_cell.lon - 0.5 * self.gridParameter['lonDelta']), "\t",
                                  str(curr_cell.lat - 0.5 * self.gridParameter['latDelta']), "\n" ] )
            ostream.writelines( [ '>', "\n" ] )
    
    # ------------------------------------------------------------------------

    def setGridParameter( self, param ):
    
        # update gridParameter dictionary
        for key in param.keys():
            self.gridParameter[key] = param[key]

    # ------------------------------------------------------------------------
    
    def setupBox( self, lonmin, lonmax, latmin, latmax, depthmin, depthmax ):

        self.grid = Grid( parentAxis = self.elementAxis,
                          elementName = 'grid' )
        
        # depth layer
        dl = DepthLayer()
        dl.add( self.grid, 'depthLayer' )
        
        if depthmax == depthmin:
            dl.at = depthmax
        else:
            dl.min = depthmin
            dl.max = depthmax

        # over lon
        for curr_lon in frange( lonmin, lonmax, self.gridParameter['lonDelta'] ):

            # over lat
            for curr_lat in frange( latmin, latmax, self.gridParameter['latDelta'] ):

                curr_cell = Cell()
                curr_cell.add( dl, 'cell' )
                
                curr_cell.lon = curr_lon
                curr_cell.lat = curr_lat

        return True

    def setupPolygon( self, polygon, depthmin, depthmax ):
        """
        polygon must be of type QPPolygon, createNodes() should already be called on polygon
        if node list is empty, call polygon.createNodes() from constructor
        """

        if not isinstance( polygon, QPPolygon ):
            raise TypeError, 'QPGrid::setupPolygon - wrong type for parameter polygon, must be QPPolygon'

        nodes = self._createNodes( polygon )

        self.grid = Grid( parentAxis = self.elementAxis,
                          elementName = 'grid' )
        
        # depth layer
        dl = DepthLayer()
        dl.add( self.grid, 'depthLayer' )

        if depthmax == depthmin:
            dl.at = depthmax
        else:
            dl.min = depthmin
            dl.max = depthmax

        # loop over nodes of polygon
        for curr_node in nodes:
            
            curr_cell     = Cell()
            curr_cell.add( dl, 'cell' )
            
            curr_cell.lon = curr_node[0]
            curr_cell.lat = curr_node[1]

        return True

    # ------------------------------------------------------------------------
    
    def inGridCell( self, lat, lon, depth ):
        
        # depthLayer bin 0...30 km: depth  0 km is inside bin
        #                           depth 30 km is outside bin
        #
        # returns ( foundLat, foundLon, foundDepthMin, foundDepthMax ) or None
    
        for curr_depth_layer in self.grid.depthLayer:
            if depth >= curr_depth_layer.min and depth < curr_depth_layer.max:

                # look for cell
                for curr_cell in curr_depth_layer.cell:
                    
                    latMin = curr_cell.lat - 0.5 * self.grid.defaultCellDimension.latRange
                    latMax = curr_cell.lat + 0.5 * self.grid.defaultCellDimension.latRange
                    lonMin = curr_cell.lon - 0.5 * self.grid.defaultCellDimension.lonRange
                    lonMax = curr_cell.lon + 0.5 * self.grid.defaultCellDimension.lonRange
                    
                    if lat >= latMin and lat < latMax and lon >= lonMin and lon < lonMax:
                        return ( curr_cell.lat, curr_cell.lon, curr_depth_layer.min, curr_depth_layer.max )
         
        # depthLayer or cell not found
        return None
            
    def inGrid( self, lat, lon, depth ):
        
        if self.inGridCell( lat, lon, depth ) is None:
            return False
        else:
            return True

    # ------------------------------------------------------------------------
    
    def _createNodes( self, polygon ):

        polygon.getExtent( self.gridParameter['lonDelta'],
                           self.gridParameter['latDelta'],
                           self.gridParameter['lonAlign'],
                           self.gridParameter['latAlign'] )

        nodes = []
        
        # loop over the grid
        for fLon in frange( polygon.lonMin, polygon.lonMax, self.gridParameter['lonDelta'] ):
            for fLat in frange( polygon.latMin, polygon.latMax, self.gridParameter['latDelta'] ):

                if self.gridParameter['includePointOnBoundary'] is True:

                    # check if point is in polygon or on boundary
                    if polygon.isInsideOrOnBoundary( fLon, fLat ):
                        
                        # Create nodes
                        nodes.append( [ fLon, fLat ] )
                else:
                    
                    # check if point in polygon
                    if polygon.isInside( fLon, fLat ):

                        # Create nodes
                        nodes.append( [ fLon, fLat ] )

        if self.gridParameter['shift'] is True:

            for nCnt in xrange( len( nodes ) ):
                if nodes[nCnt][0] > 180.0:
                    nodes[nCnt][0] = nodes[nCnt][0] - 360.0

        return nodes





