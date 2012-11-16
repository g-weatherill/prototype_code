# -*- coding: iso-8859-1 -*-
#
# quakepy/CreationInfo.py
# $Id: CreationInfo.py 157 2009-02-16 11:15:52Z fab $
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

__version__  = '$Id: CreationInfo.py 157 2009-02-16 11:15:52Z fab $'
__revision__ = '$Revision: 157 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

from QPCore import *
from QPDateTime import *

class CreationInfo( QPObject ):
    """
    QuakePy: CreationInfo
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'agencyID', 'agencyID', 'element', unicode, 'basic' ),
                    QPElement( 'agencyURI', 'agencyURI', 'element', unicode, 'basic' ),
                    QPElement( 'author', 'author', 'element', unicode, 'basic' ),
                    QPElement( 'authorURI', 'authorURI', 'element', unicode, 'basic' ),
                    QPElement( 'version', 'version', 'element', unicode, 'basic' ),
                    QPElement( 'creationTime', 'creationTime', 'element', QPDateTime, 'basic' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self,
                  agencyID = None,
                  agencyURI = None,
                  author = None,
                  authorURI = None,
                  creationTime = None,
                  version = None,
                  **kwargs ):
        super( CreationInfo, self ).__init__( **kwargs )

        self.elements.extend( self.addElements )

        self._initMultipleElements()
        
        self.agencyID     = agencyID
        self.agencyURI    = agencyURI
        self.author       = author
        self.authorURI    = authorURI
        self.version      = version

        if creationTime is not None:
            self.creationTime = QPDateTime( creationTime )
        
