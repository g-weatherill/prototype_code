# -*- coding: iso-8859-1 -*-
#
# quakepy/Magnitude.py
# $Id: Magnitude.py 208 2009-06-19 14:53:10Z fab $
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

__version__  = '$Id: Magnitude.py 208 2009-06-19 14:53:10Z fab $'
__revision__ = '$Revision: 208 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

from QPCore import *

from StationMagnitudeReference  import StationMagnitudeReference
from CreationInfo               import CreationInfo
from Comment                    import Comment

from RealQuantity    import RealQuantity

class Magnitude( QPPublicObject ):
    """
    QuakePy: Magnitude
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'publicID', 'publicID', 'attribute', unicode, 'basic' ),
        
                    QPElement( 'type', 'type', 'element', unicode, 'basic' ),
                    QPElement( 'originID', 'originID', 'element', unicode, 'basic' ),
                    QPElement( 'methodID', 'methodID', 'element', unicode, 'basic' ),
                    QPElement( 'stationCount', 'stationCount', 'element', int, 'basic' ),
                    QPElement( 'azimuthalGap', 'azimuthalGap', 'element', float, 'basic' ),

                    QPElement( 'mag', 'mag', 'element', RealQuantity, 'complex' ),
                    QPElement( 'creationInfo', 'creationInfo', 'element', CreationInfo, 'complex' ),
                          
                    QPElement( 'stationMagnitudeReference', 'stationMagnitudeReference', 'element', StationMagnitudeReference, 'multiple' ),
                    QPElement( 'comment', 'comment', 'element', Comment, 'multiple' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, publicID = None, originID = None, **kwargs  ):
        super( Magnitude, self ).__init__( publicID, **kwargs  )

        self.elements.extend( self.addElements )
                        
        # publicID has not been set in parent class
        if self.publicID is None:
            self.publicID = self.createPublicID( self.__class__.__name__, **kwargs )
        self.setOriginAssociation( originID )

        self._initMultipleElements()
        
    def setOriginAssociation( self, originID ):
        self.originID = originID
        