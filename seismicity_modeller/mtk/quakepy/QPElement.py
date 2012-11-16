# -*- coding: iso-8859-1 -*-
#
# quakepy/QPElement.py
# $Id: QPElement.py 159 2009-02-17 13:34:05Z fab $
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

__version__  = '$Id: QPElement.py 159 2009-02-17 13:34:05Z fab $'
__revision__ = '$Revision: 159 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import pyRXP

from QPCore import *

POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)

class QPElementList( list ):
    pass

class QPElement( object ):
    """
    QuakePy: QPElement

    defines characteristics of attributes a class derived from QPObject has
    """
    
    def __init__( self,
                  varname,
                  xmlname,
                  xmltype,
                  pytype,
                  vartype,
                  parentaxis = None,
                  parenttype = None ):

        if varname is not None:
            self.varname = varname
        else:
            error_str = "QPElement constructor: varname must not be None"
            raise ValueError, error_str

        if xmlname is not None:
            self.xmlname = xmlname
        else:
            error_str = "QPElement constructor: xmlname must not be None"
            raise ValueError, error_str

        if xmltype is not None:
            self.xmltype = xmltype
        else:
            error_str = "QPElement constructor: xmltype must not be None"
            raise ValueError, error_str

        if pytype is not None:
            self.pytype = pytype
        else:
            error_str = "QPElement constructor: pytype must not be None"
            raise ValueError, error_str

        if vartype is not None:
            self.vartype = vartype
        else:
            error_str = "QPElement constructor: vartype must not be None"
            raise ValueError, error_str

        self.parentaxis = parentaxis
        self.parenttype = parenttype

        
       
