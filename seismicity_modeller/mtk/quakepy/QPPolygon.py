# -*- coding: iso-8859-1 -*-
#
# quakepy/QPPolygon.py
# $Id: QPPolygon.py 255 2009-12-01 17:55:22Z fab $
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

__version__  = '$Id: QPPolygon.py 255 2009-12-01 17:55:22Z fab $'
__revision__ = '$Revision: 255 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import numpy

from shapely.geometry import Polygon, Point

from QPUtils import *

class QPPolygon( object ):
    """
    QuakePy: QPPolygon

    QPPolygon implements a 2-dim polygon, using the Shapely package
    """

    def __init__( self, vertices, **kwargs ):
        """
        vertices can be either filename or stream from which data is read (2-column ASCII, lon | lat),
        or iterable of two components, e.g.
        [ [,], [,], ... ] or ( (,), (,), ... )

        kwargs: swapCoordinates=True: input has lat|lon pairs (default: False)
                verbose_level       : 0 - no information, 1 - verbose
        """

        # Shapely polygon object
        self.polygon = None

        # polygon vertices
        self.vertices = []

        # polygon extent
        self.extent = {}
    
        if isinstance( vertices, ( list, tuple ) ):
            self.vertices = vertices
        else:
            try:
                input_stream = getQPDataSource( vertices, **kwargs )
            except:
                raise IOError, "QPPolygon - cannot open data stream from %s" % ( vertices )

            line_ctr = 0
            vertices_added = 0
            v_old    = None
            
            for line in input_stream:

                line_ctr += 1

                if len( line.strip() ) > 0:

                    v = line.strip().split()

                    # check if line is not duplicated
                    if len( v ) >= 2:

                        try:
                            v0 = float( v[0] )
                            v1 = float( v[1] )
                        except:
                            error_str = "QPPolygon - non-float value in input line no. %s: %s" % ( line_ctr, line )
                            raise ValueError, error_str

                        # add polygon vertex if it differs from last one
                        if ( v_old is None ) or not ( floatEqual( v0, v_old[0] ) and floatEqual( v1, v_old[1] ) ):
                    
                            if 'swapCoordinates' in kwargs.keys() and kwargs['swapCoordinates'] is True:
                                self.vertices.append( [ v1, v0 ] )
                            else:
                                self.vertices.append( [ v0, v1 ] )

                            vertices_added += 1
                            v_old = [ v0, v1 ]

                        else:
                            # go to next line
                            continue

                    else:
                        error_str = "QPPolygon - error in input line no. %s: %s" % ( line_ctr, line )
                        raise ValueError, error_str
                else:
                    # blank line, go to next
                    continue

            if 'verbose_level' in kwargs.keys() and kwargs['verbose_level'] > 0:
                print "QPPolygon - processed %s input lines, added %s vertices" % ( line_ctr, vertices_added )

        if len( self.vertices ) > 0:
            self.polygon = Polygon( self.vertices )
            self.getExtent()

    def getExtent( self ):
        
        # Determine extent of polygon bounding box
        self.extent = { 'lonMin': self.polygon.bounds[0],
                        'lonMax': self.polygon.bounds[2],
                        'latMin': self.polygon.bounds[1],
                        'latMax': self.polygon.bounds[3] }
        
    def isInside( self, fX, fY ):
        """
        check if point( fX, fY ) is truly inside polygon (not on boundary)
        'true'  if inside
        'false' if outside or on boundary
        """

        if self.polygon.contains( Point( fX, fY ) ):
            return True
        else:
            return False

    def isInsideOrOnBoundary( self, fX, fY ):
        """
        check if point( fX, fY ) is inside polygon or on boundary
        'true'  if inside or on boundary
        'false' if outside
        """

        if self.isInside( fX, fY ) or self.isOnBoundary( fX, fY ):
            return True
        else:
            return False

    def isOnBoundary( self, fX, fY ):

        if self.polygon.boundary.contains( Point( fX, fY ) ):
            return True
        else:
            return False
        