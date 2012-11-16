# -*- coding: iso-8859-1 -*-
#
# quakepy/IntegerQuantity.py
# $Id: IntegerQuantity.py 157 2009-02-16 11:15:52Z fab $
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

__version__  = '$Id: IntegerQuantity.py 157 2009-02-16 11:15:52Z fab $'
__revision__ = '$Revision: 157 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import pyRXP

from QPCore import *

POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)

class IntegerQuantity( QPObject ):
    """
    QuakePy: IntegerQuantity
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'value', 'value', 'element', int, 'basic' ),
                    QPElement( 'uncertainty', 'uncertainty', 'element', int, 'basic' ),
                    QPElement( 'lowerUncertainty', 'lowerUncertainty', 'element', int, 'basic' ),
                    QPElement( 'upperUncertainty', 'upperUncertainty', 'element', int, 'basic' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, val = None, unc = None, unc_lower = None, unc_upper = None, **kwargs ):
        super( IntegerQuantity, self ).__init__( **kwargs )

        self.elements.extend( self.addElements )
                        
        self._initMultipleElements()
        
        self.value            = ( val != None and [int(val)] or [None] )[0]
        self.uncertainty      = ( unc != None and [int(unc)] or [None] )[0]
        self.lowerUncertainty = ( unc_lower != None and [int(unc_lower)] or [None] )[0]
        self.upperUncertainty = ( unc_upper != None and [int(unc_upper)] or [None] )[0]
       