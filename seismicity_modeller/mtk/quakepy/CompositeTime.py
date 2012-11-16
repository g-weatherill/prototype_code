# -*- coding: iso-8859-1 -*-
#
# quakepy/CompositeTime.py
# $Id: CompositeTime.py 157 2009-02-16 11:15:52Z fab $
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

__version__  = '$Id: CompositeTime.py 157 2009-02-16 11:15:52Z fab $'
__revision__ = '$Revision: 157 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

from QPCore import *

from RealQuantity    import RealQuantity
from IntegerQuantity import IntegerQuantity

class CompositeTime( QPObject ):
    """
    QuakePy: CompositeTime
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'differenceToUTC', 'differenceToUTC', 'element', float, 'basic' ),

                    QPElement( 'year', 'year', 'element', IntegerQuantity, 'complex' ),
                    QPElement( 'month', 'month', 'element', IntegerQuantity, 'complex' ),
                    QPElement( 'day', 'day', 'element', IntegerQuantity, 'complex' ),
                    QPElement( 'hour', 'hour', 'element', IntegerQuantity, 'complex' ),
                    QPElement( 'minute', 'minute', 'element', IntegerQuantity, 'complex' ),
                    QPElement( 'second', 'second', 'element', RealQuantity, 'complex' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, **kwargs ):
        super( CompositeTime, self ).__init__( **kwargs )

        self.elements.extend( self.addElements )
                        
        self._initMultipleElements()
        
