# -*- coding: iso-8859-1 -*-
#
# quakepy/Comment.py
# $Id: Comment.py 157 2009-02-16 11:15:52Z fab $
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

__version__  = '$Id: Comment.py 157 2009-02-16 11:15:52Z fab $'
__revision__ = '$Revision: 157 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

from QPCore import *

from CreationInfo import CreationInfo


class Comment( QPObject ):
    """
    QuakePy: Comment
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'id', 'id', 'attribute', unicode, 'basic' ),
                    QPElement( 'text', 'text', 'element', unicode, 'basic' ),
                    QPElement( 'creationInfo', 'creationInfo', 'element', CreationInfo, 'complex' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, text = None, id = None, creationInfo = None, **kwargs ):
        super( Comment, self ).__init__( **kwargs )

        self.elements.extend( self.addElements )
                        
        self.text         = text
        self.id           = id               # id is an XML attribute
        self.creationInfo = creationInfo
        
