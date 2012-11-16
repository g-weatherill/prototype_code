# -*- coding: iso-8859-1 -*-
#
# quakepy/Origin.py
# $Id: Origin.py 208 2009-06-19 14:53:10Z fab $
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

__version__  = '$Id: Origin.py 208 2009-06-19 14:53:10Z fab $'
__revision__ = '$Revision: 208 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

from QPCore import *

from Comment                    import Comment
from CreationInfo               import CreationInfo
from OriginUncertainty          import OriginUncertainty
from OriginQuality              import OriginQuality 
from Arrival                    import Arrival
from StationMagnitude           import StationMagnitude

from RealQuantity               import RealQuantity
from IntegerQuantity            import IntegerQuantity
from TimeQuantity               import TimeQuantity
from CompositeTime              import CompositeTime


class Origin( QPPublicObject ):
    """
    QuakePy: Origin
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'publicID', 'publicID', 'attribute', unicode, 'basic' ),
        
                    QPElement( 'referenceSystemID', 'referenceSystemID', 'element', unicode, 'basic' ),
                    QPElement( 'methodID', 'methodID', 'element', unicode, 'basic' ),
                    QPElement( 'earthModelID', 'earthModelID', 'element', unicode, 'basic' ),

                    QPElement( 'depthType', 'depthType', 'element', unicode, 'enum' ),
                    QPElement( 'type', 'type', 'element', unicode, 'enum' ),
                    QPElement( 'evaluationMode', 'evaluationMode', 'element', unicode, 'enum' ),
                    QPElement( 'evaluationStatus', 'evaluationStatus', 'element', unicode, 'enum' ),

                    QPElement( 'time', 'time', 'element', TimeQuantity, 'complex' ),
                    QPElement( 'longitude', 'longitude', 'element', RealQuantity, 'complex' ),
                    QPElement( 'latitude', 'latitude', 'element', RealQuantity, 'complex' ),
                    QPElement( 'depth', 'depth', 'element', RealQuantity, 'complex' ),
                    QPElement( 'quality', 'quality', 'element', OriginQuality, 'complex' ),
                    QPElement( 'creationInfo', 'creationInfo', 'element', CreationInfo, 'complex' ),
                          
                    QPElement( 'originUncertainty', 'originUncertainty', 'element', OriginUncertainty, 'multiple' ),
                    QPElement( 'arrival', 'arrival', 'element', Arrival, 'multiple' ),
                    QPElement( 'stationMagnitude', 'stationMagnitude', 'element', StationMagnitude, 'multiple' ),
                    QPElement( 'compositeTime', 'compositeTime', 'element', CompositeTime, 'multiple' ),
                    QPElement( 'comment', 'comment', 'element', Comment, 'multiple' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, publicID = None, **kwargs  ):
        super( Origin, self ).__init__( publicID, **kwargs  )

        self.elements.extend( self.addElements )
                        
        # publicID has not been set in parent class
        if self.publicID is None:
            self.publicID = self.createPublicID( self.__class__.__name__, **kwargs )

        self._initMultipleElements()

    def getArrivalsIdx( self, pick ):
        arr_idx = []
        
        for curr_arr_idx, curr_arr in enumerate( self.arrival ):
            if curr_arr.pickID == pick.publicID:
                arr_idx.append( curr_arr_idx )
        return arr_idx
    
    def getArrivals( self, pick ):
        arrivals = []
        
        for arr in self.arrival:
            if arr.pickID == pick.publicID:
                arrivals.append( arr )
        return arrivals