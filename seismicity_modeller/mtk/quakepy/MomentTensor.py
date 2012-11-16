# -*- coding: utf-8 -*-
#
# quakepy/MomentTensor.py
# $Id: MomentTensor.py 266 2010-09-21 12:18:02Z fab $
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

__version__  = '$Id: MomentTensor.py 266 2010-09-21 12:18:02Z fab $'
__revision__ = '$Revision: 266 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

from QPCore import *

from Tensor import Tensor
from DataUsed import DataUsed
from SourceTimeFunction import SourceTimeFunction
from CreationInfo       import CreationInfo
from Comment            import Comment

from RealQuantity    import RealQuantity


class MomentTensor( QPObject ):
    """
    QuakePy: MomentTensor
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'derivedOriginID', 'derivedOriginID', 'attribute', unicode, 'basic' ),
                    QPElement( 'triggeringOriginID', 'triggeringOriginID', 'attribute', unicode, 'basic' ),
                    QPElement( 'momentMagnitudeID', 'momentMagnitudeID', 'attribute', unicode, 'basic' ),
                    QPElement( 'variance', 'variance', 'element', float, 'basic' ),
                    QPElement( 'varianceReduction', 'varianceReduction', 'element', float, 'basic' ),
                    QPElement( 'doubleCouple', 'doubleCouple', 'element', float, 'basic' ),
                    QPElement( 'clvd', 'clvd', 'element', float, 'basic' ),
                    QPElement( 'iso', 'iso', 'element', float, 'basic' ),
                    QPElement( 'greensFunctionID', 'greensFunctionID', 'element', unicode, 'basic' ),
                    QPElement( 'filterID', 'filterID', 'element', unicode, 'basic' ),
                    QPElement( 'methodID', 'methodID', 'attribute', unicode, 'basic' ),
                    QPElement( 'cmtName', 'cmtName', 'element', unicode, 'basic' ),
                    QPElement( 'cmtVersion', 'cmtVersion', 'element', unicode, 'basic' ),
                          
                    QPElement( 'method', 'method', 'element', unicode, 'enum' ),
                    QPElement( 'status', 'status', 'element', unicode, 'enum' ),

                    QPElement( 'scalarMoment', 'scalarMoment', 'element', RealQuantity, 'complex' ),
                    QPElement( 'tensor', 'tensor', 'element', Tensor, 'complex' ),
                    QPElement( 'sourceTimeFunction', 'sourceTimeFunction', 'element', SourceTimeFunction, 'complex' ),
                    QPElement( 'creationInfo', 'creationInfo', 'element', CreationInfo, 'complex' ),
                          
                    QPElement( 'dataUsed', 'dataUsed', 'element', DataUsed, 'multiple' ),
                    QPElement( 'comment', 'comment', 'element', Comment, 'multiple' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, **kwargs ):
        super( MomentTensor, self ).__init__( **kwargs )
        
        self.elements.extend( self.addElements )

        self._initMultipleElements()

        
