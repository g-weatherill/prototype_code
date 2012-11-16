#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
#
# quakepy/test/QPUnitTest.py
# $Id: QPUnitTest.py 228 2009-09-29 22:11:03Z fab $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2007 by Fabian Euchner and Danijel Schorlemmer          #
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

__version__  = '$Id: QPUnitTest.py 228 2009-09-29 22:11:03Z fab $'
__revision__ = '$Revision: 228 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import sys
import unittest

sys.path.append( '../../..' )
sys.path.append( '../..' )
sys.path.append( '..' )
sys.path.append( 'test/unitTest' )
sys.path.append( 'unitTest' )

# invoke unit tests
from QPUtilsTest import QPUtilsTest
from QPCatalogTest import QPCatalogTest
#from QPDateTimeTest import QPDateTimeTest
#from PMCTest import PMCTest

# Invoke the module
if __name__ == '__main__':
   
   # Invoke all tests
   unittest.main()
        
# end of main
