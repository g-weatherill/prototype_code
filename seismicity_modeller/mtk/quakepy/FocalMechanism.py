# -*- coding: utf-8 -*-
#
# quakepy/FocalMechanism.py
# $Id: FocalMechanism.py 266 2010-09-21 12:18:02Z fab $
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

__version__  = '$Id: FocalMechanism.py 266 2010-09-21 12:18:02Z fab $'
__revision__ = '$Revision: 266 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import pyRXP

from QPCore import *

from NodalPlanes     import NodalPlanes
from PrincipalAxes   import PrincipalAxes
from CreationInfo    import CreationInfo
from Comment         import Comment

from MomentTensor    import MomentTensor
from RealQuantity    import RealQuantity
from IntegerQuantity import IntegerQuantity
from TimeQuantity    import TimeQuantity

POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)

class FocalMechanism( QPPublicObject ):
    """
    QuakePy: FocalMechanism
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'publicID', 'publicID', 'attribute', unicode, 'basic' ),
        
                    QPElement( 'triggeringOriginID', 'triggeringOriginID', 'attribute', unicode, 'basic' ),
                    QPElement( 'azimuthalGap', 'azimuthalGap', 'element', float, 'basic' ),
                    QPElement( 'stationPolarityCount', 'stationPolarityCount', 'element', int, 'basic' ),
                    QPElement( 'misfit', 'misfit', 'element', float, 'basic' ),
                    QPElement( 'stationDistributionRatio', 'stationDistributionRatio', 'element', float, 'basic' ),
                    QPElement( 'methodID', 'methodID', 'attribute', unicode, 'basic' ),

                    QPElement( 'nodalPlanes', 'nodalPlanes', 'element', NodalPlanes, 'complex' ),
                    QPElement( 'principalAxes', 'principalAxes', 'element', PrincipalAxes, 'complex' ),
                    QPElement( 'creationInfo', 'creationInfo', 'element', CreationInfo, 'complex' ),
                          
                    QPElement( 'momentTensor', 'momentTensor', 'element', MomentTensor, 'multiple' ),
                    QPElement( 'comment', 'comment', 'element', Comment, 'multiple' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, publicID = None, **kwargs  ):
        super( FocalMechanism, self ).__init__( publicID, **kwargs  )

        self.elements.extend( self.addElements )
                        
        # publicID has not been set in parent class
        if self.publicID is None:
            self.publicID = self.createPublicID( self.__class__.__name__, **kwargs )

        self._initMultipleElements()

        
