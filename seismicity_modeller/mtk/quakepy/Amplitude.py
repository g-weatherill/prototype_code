# -*- coding: iso-8859-1 -*-
#
# quakepy/Amplitude.py
# $Id: Amplitude.py 157 2009-02-16 11:15:52Z fab $
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

__version__  = '$Id: Amplitude.py 157 2009-02-16 11:15:52Z fab $'
__revision__ = '$Revision: 157 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

from QPCore import *

from RealQuantity     import RealQuantity
from IntegerQuantity  import IntegerQuantity
from TimeQuantity     import TimeQuantity
from TimeWindow       import TimeWindow
from WaveformStreamID import WaveformStreamID
from CreationInfo     import CreationInfo
from Comment          import Comment

class Amplitude( QPPublicObject ):
    """
    QuakePy: Amplitude
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'publicID', 'publicID', 'attribute', unicode, 'basic' ),
                    
                    QPElement( 'type', 'type', 'element', unicode, 'basic' ),
                    QPElement( 'pickID', 'pickID', 'element', unicode, 'basic' ),
                    QPElement( 'filterID', 'filterID', 'element', unicode, 'basic' ),
                    QPElement( 'methodID', 'methodID', 'element', unicode, 'basic' ),
                    QPElement( 'magnitudeHint', 'magnitudeHint', 'element', unicode, 'basic' ),

                    QPElement( 'evaluationMode', 'evaluationMode', 'element', unicode, 'enum' ),

                    QPElement( 'amp', 'amp', 'element', RealQuantity, 'complex' ),
                    QPElement( 'period', 'period', 'element', RealQuantity, 'complex' ),
                    QPElement( 'waveformID', 'waveformID', 'element', WaveformStreamID, 'complex' ),
                    QPElement( 'timeWindow', 'timeWindow', 'element', TimeWindow, 'complex' ),
                    QPElement( 'scalingTime', 'scalingTime', 'element', TimeQuantity, 'complex' ),
                    QPElement( 'creationInfo', 'creationInfo', 'element', CreationInfo, 'complex' ),
                          
                    QPElement( 'comment', 'comment', 'element', Comment, 'multiple' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, publicID = None, **kwargs ):
        super( Amplitude, self ).__init__( publicID, **kwargs )

        self.elements.extend( self.addElements )
                        
        # publicID has not been set in parent class
        if self.publicID is None:
            self.publicID = self.createPublicID( self.__class__.__name__, **kwargs )

        self._initMultipleElements()

        
