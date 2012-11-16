#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
#
# quakepy/test/unitTest/QPUtilsTest.py
# $Id: QPUtilsTest.py 228 2009-09-29 22:11:03Z fab $
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

__version__  = '$Id: QPUtilsTest.py 228 2009-09-29 22:11:03Z fab $'
__revision__ = '$Revision: 228 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import sys
import shutil
import os
import unittest
import datetime

import mx.DateTime

sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')

from QPTestCase import QPTestCase

from QPUtils import *


class QPUtilsTest( QPTestCase ):

    ## static data of the class

    # unit tests use sub-directory of global reference data directory
    __referenceDataDir = os.path.join( QPTestCase.ReferenceDataDir,
                                       'unitTest', 'qputils' )
                                       

    def testMxDateTime2ISO( self ):
        """
        test mxDateTime2ISO function
        """
        print
        print " ----- testMxDateTime2ISO -----"

        # setup test name
        QPTestCase.setTestName( self, "QPUtils-MxDateTime2ISO" )

        # cd to the test directory, remember current directory
        cwd = os.getcwd()
        os.chdir( QPTestCase.TestDirPath )

        try:
    
            # test date
            testDate = mx.DateTime.DateTime( 2007, 1, 1, 12, 0, 0.12345678 )

            self.failIf( '2007-01-01T12:00:00.12345678' != mxDateTime2ISO( testDate ),
                         "error: testing without options" )

            self.failIf( '2007-01-01T12:00:00.1234' != mxDateTime2ISO( testDate, secondsdigits=4 ),
                         "error: using option secondsdigits " )

            self.failIf( '2007-01-01T12:00:00.1235' != mxDateTime2ISO( testDate, secondsdigits=4, round=True ),
                         "error: using options secondsdigits and round" )

            self.failIf( '2007-01-01T12:00:00' != mxDateTime2ISO( testDate, showsecfrac=False ),
                         "error: using option showsecfrac" )

            self.failIf( '2007-01-01' != mxDateTime2ISO( testDate, showtime=False ),
                         "error: using option showtime" )

            self.failIf( '2007-01-01$12:00:00' != mxDateTime2ISO( testDate, showsecfrac=False, partsepreplacechar='$' ),
                         "error: using options showsecfrac and partsepreplacechar" )

            self.failIf( '2007-01-01T12!00!00' != mxDateTime2ISO( testDate, showsecfrac=False, timesepreplacechar='!' ),
                         "error: using options showsecfrac and timesepreplacechar" )

            self.failIf( '2007=01=01T12:00:00' != mxDateTime2ISO( testDate, showsecfrac=False, datesepreplacechar='=' ),
                         "error: using options showsecfrac and datesepreplacechar" )
                         
        finally:
            # return to the original directory
            os.chdir( cwd )
        

    def testDecimalYear( self ):
        """
        test decimalYear function
        """
        print
        print " ----- testDecimalYear -----"

        # setup test name
        QPTestCase.setTestName( self, "QPUtils-DecimalYear" )

        # cd to the test directory, remember current directory
        cwd = os.getcwd()
        os.chdir( QPTestCase.TestDirPath )

        try:
    
            self.failIf( 2000.0 != decimalYear( mx.DateTime.DateTime( 2000, 1, 1 ) ),
                         "error: testing date 2000-01-01" )

            self.failIf( 2000.5 != decimalYear( mx.DateTime.DateTime( 2000, 7, 2 ) ),
                         "error: testing date 2000-07-02" )
                         
            self.failIf( 2001.5 != decimalYear( mx.DateTime.DateTime( 2001, 7, 2, 12, 0, 0 ) ),
                         "error: testing date 2001-07-02T12:00:00" )
                         
        finally:
            # return to the original directory
            os.chdir( cwd )


    def testFromDecimalYear( self ):
        """
        test fromDecimalYear function
        """
        print
        print " ----- testFromDecimalYear -----"

        # setup test name
        QPTestCase.setTestName( self, "QPUtils-FromDecimalYear" )

        # cd to the test directory, remember current directory
        cwd = os.getcwd()
        os.chdir( QPTestCase.TestDirPath )

        try:
    
            self.failIf( fromDecimalYear( 2000.0 ) != mx.DateTime.DateTime( 2000, 1, 1 ),
                         "error: testing date 2000-01-01" )

            self.failIf( fromDecimalYear( 2000.5 ) != mx.DateTime.DateTime( 2000, 7, 2 ),
                         "error: testing date 2000-07-02" )
                         
            self.failIf( fromDecimalYear( 2001.5 ) != mx.DateTime.DateTime( 2001, 7, 2, 12, 0, 0 ),
                         "error: testing date 2001-07-02T12:00:00" )
                         
        finally:
            # return to the original directory
            os.chdir( cwd )


    def testFixTimeComponents( self ):
        """
        test fixTimeComponents function
        """
        print
        print " ----- testFixTimeComponents -----"

        # setup test name
        QPTestCase.setTestName( self, "QPUtils-FixTimeComponents" )

        # cd to the test directory, remember current directory
        cwd = os.getcwd()
        os.chdir( QPTestCase.TestDirPath )

        try:
    
            self.failIf( { 'component': ( 0, 0, 0.0 ), 'increaseFlag': ( True, True, True ) } != fixTimeComponents( 24, 60, 60.0 ),
                         "error: checking hour 24, min 60, sec 60.0" )
                         
        finally:
            # return to the original directory
            os.chdir( cwd )


    def testAdjustDateTime( self ):
        """
        test adjustDateTime function
        """
        print
        print " ----- testAdjustDateTime -----"

        # setup test name
        QPTestCase.setTestName( self, "QPUtils-AdjustDateTime" )

        # cd to the test directory, remember current directory
        cwd = os.getcwd()
        os.chdir( QPTestCase.TestDirPath )

        try:

            self.failIf( DateTime( 1999, 3, 2, 1, 0, 59.5  ) != adjustDateTime( ( True, True, True ),
                                                                                DateTime( 1999, 2, 28, 23, 59, 59.5 ) ),
                         "error: checking date 1999-02-28T23:59:59.9" )
                         
        finally:
            # return to the original directory
            os.chdir( cwd )


    def testFrange( self ):
        """
        test frange function
        """
        print
        print " ----- testFrange -----"

        # setup test name
        QPTestCase.setTestName( self, "QPUtils-Frange" )

        # cd to the test directory, remember current directory
        cwd = os.getcwd()
        os.chdir( QPTestCase.TestDirPath )

        try:

            refArray   = [ -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1 ]
            checkArray = frange( -0.1, 1.1, 0.1 ).tolist()
            
            self.failIf( len( refArray ) != len( checkArray ),
                         "error: length is not the same for both lists" )

            for idx in xrange( len( refArray ) ):
                self.failIf( floatEqual( refArray[idx], checkArray[idx] ) is False,
                             "mismatch in list element %s: %s, %s" %  ( idx, refArray[idx], checkArray[idx] ) )
                         
        finally:
            # return to the original directory
            os.chdir( cwd )


    def testLocateInArray( self ):
        """
        test locateInArray function
        """
        print
        print " ----- testLocateInArray -----"

        # setup test name
        QPTestCase.setTestName( self, "QPUtils-LocateInArray" )

        # cd to the test directory, remember current directory
        cwd = os.getcwd()
        os.chdir( QPTestCase.TestDirPath )

        try:

            # test array
            testArray = [ -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1 ]
            N         = len( testArray )

            # check lower out of bounds (-1.5): idx = -1
            self.failIf( -1 != locateInArray( testArray, -1.5 ),
                            "error locating -1.5, lower out of bounds" )

            # check lower bound (-0.1): idx = 0
            self.failIf( 0 != locateInArray( testArray, -0.1 ),
                            "error locating -0.1, value on lower bound" )

            # check regular value (0.25): idx = 3
            self.failIf( 3 != locateInArray( testArray, 0.25 ),
                            "error locating 0.25, regular value" )

            # check value on grid (0.4): idx = 5
            self.failIf( 5 != locateInArray( testArray, 0.4 ),
                            "error locating 0.4, value on grid" )

            # check value in last bin (1.05): idx = 11
            self.failIf( 11 != locateInArray( testArray, 1.05 ),
                            "error locating 1.05, value in last bin" )

            # check upper bound (1.1): idx = 11
            # note the exception in assignment of returned index
            self.failIf( 11 != locateInArray( testArray, 1.1 ),
                            "error locating 1.1, value on upper bound (exception)" )

            # check upper out of bounds (1.35): idx = N-1
            self.failIf( N-1 != locateInArray( testArray, 1.35 ),
                            "error locating 1.35, upper out of bounds" )

        finally:
            # return to the original directory
            os.chdir( cwd )


# Invoke the module
if __name__ == '__main__':
   
   # Invoke all tests
   unittest.main()
        
# end of main