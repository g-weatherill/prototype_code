# -*- coding: iso-8859-1 -*-
#
# quakepy/test/QPTestCase.py
# $Id: QPTestCase.py 70 2008-05-28 16:42:59Z fab $
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

__version__  = '$Id: QPTestCase.py 70 2008-05-28 16:42:59Z fab $'
__revision__ = '$Revision: 70 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import datetime
import os
import unittest

class QPTestCase( unittest.TestCase ):

    ## static data
    
    # test date, initialize with current date/time (UTC)
    Date = datetime.datetime.utcnow()
    
    # top-level directory for the tests
    # is directory of this module
    __Dir = os.path.dirname( __file__ )
    
    # directory that stores reference data (e.g., EQ catalogs)
    ReferenceDataDir = os.path.join( __Dir, 'data' )

    # directory for execution of particular test and storing result data
    TestDirPath = os.path.join( __Dir, 'results' )

    # should test result directory be kept?
    KeepTestDir = True
    
    # name of the test currently being invoked
    __TestName = ''
    
    def setUp( self ):
        """
        setup test environment
        create test directory if not there
        """
        
        if os.path.exists( self.TestDirPath ) is False:
            os.mkdir( self.TestDirPath )

    def tearDown( self ):
        """
        clean up test environment
        remove test directory *or* move to a unique directory name
        """

        if os.path.exists( self.TestDirPath ) is True:
                       
            # move test directory to unique location
            if self.KeepTestDir is True:
                new_dir = "%s-%s" % ( self.TestDirPath, self.__TestName )
                new_location = os.path.join( self.__Dir, new_dir )
                os.rename( self.TestDirPath, new_location )
            else:
                os.remove( self.TestDirPath )

    def setTestName( self, name ):
        """
        set name for test
        will be called from specific test module
        """
        self.__TestName = name
        