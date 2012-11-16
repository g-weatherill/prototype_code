# -*- coding: iso-8859-1 -*-
#
# quakepy/QPDateTime.py
# $Id: QPDateTime.py 174 2009-03-26 23:32:57Z fab $
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

__version__  = '$Id: QPDateTime.py 174 2009-03-26 23:32:57Z fab $'
__revision__ = '$Revision: 174 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import sys

sys.path.append('..')
sys.path.append('../..')
sys.path.append('../../..')

import mx.DateTime

import QPUtils

class QPDateTime( object ):
    """
    QuakePy date/time class, wrapped around mx.DateTime object
    added methods for comparison and string representation
    """
    
    # standard accuracy (fraction of seconds) for datetime comparison in QuakePy
    cmpEpsilon = 1.0e-09
    
    # standard seconds format string: 6 digits = microseconds
    secondsDigits = 6
    
    def __init__( self, datetime_in = None, **kwargs ):
        """
        initialize with time stamp, which can be
        
        (1) mx.DateTime object
        (2) date/time string in ISO format (understandable by mx.DateTime ISO method) 
        (3) list or tuple with 6 elements: ( year, month, day, hour, minute, second )
        """
        
        if datetime_in is not None:
          
            if isinstance( datetime_in, basestring ):
                try:
                    self.datetime = mx.DateTime.ISO.ParseDateTimeUTC( datetime_in )
                except:
                    error_msg = "QPDateTime constructor: input value not in valid ISO format - %s" % datetime_in
                    raise ValueError, error_msg
            
            elif isinstance( datetime_in, mx.DateTime.DateTimeType ):
                self.datetime = datetime_in
                
            elif (     ( isinstance( datetime_in, tuple ) or isinstance( datetime_in, list ) ) 
                   and ( len( datetime_in ) >= 6 ) ):
                try:
                    self.datetime = mx.DateTime.DateTime( int( datetime_in[0] ), 
                                                          int( datetime_in[1] ), 
                                                          int( datetime_in[2] ), 
                                                          int( datetime_in[3] ), 
                                                          int( datetime_in[4] ), 
                                                          float( datetime_in[5] ) )
                except:
                    error_msg = "QPDateTime constructor: input datetime array not valid - %s" % datetime_in
                    raise ValueError, error_msg
                
            else:
                error_msg = "QPDateTime constructor: no valid input format - %s %s" % ( datetime_in, type( datetime_in ) )
                raise ValueError, error_msg
            
        if ( 'epsilon' in kwargs.keys() ) and isinstance( kwargs['epsilon'], float ):
            self.cmpEpsilon = kwargs['epsilon']
            
        if ( 'secondsdigits' in kwargs.keys() ) and isinstance( kwargs['secondsdigits'], int ):
            self.secondsDigits = kwargs['secondsdigits']
        
    def setEpsilon( self, epsilon ):
        if isinstance( epsilon, float ):
            self.cmpEpsilon = epsilon
            
    def setSecondsDigits( self, digits ):
        if isinstance( digits, int ):
            self.secondsDigits = digits
            
    def diff( self, T ):
        """
        return difference between self and other QPDateTime object
        return value is of type mx.DateTime.DateTimeDelta
        """
        return self.datetime - T.datetime
    
    def toISO( self, **kwargs ):
        """
        return datetime as customisable ISO string
        uses function mxDateTime2ISO() from module QPUtils
        default: print decimal places for seconds as given in self.secondsDigits, no rounding
        """
        if not 'secondsdigits' in kwargs.keys():
            kwargs['secondsdigits'] = self.secondsDigits
          
        return QPUtils.mxDateTime2ISO( self.datetime, **kwargs )

    def toDecimalYear( self ):
        """
        return decimal year / floating point representation of self.datetime
        """
        return QPUtils.decimalYear( self.datetime )
        
    def __str__( self ):
        """
        return datetime as default ISO string representation
        """
        return self.toISO()
    
    def __unicode__( self ):
        """
        return datetime as unicode representation of default ISO string
        """
        return unicode( self.toISO() )

    ## comparison methods

    def __cmp__( self, T ):
        return mx.DateTime.cmp( self.datetime, T.datetime, self.cmpEpsilon )
    
    def __eq__( self, T ):
        if self.__cmp__( T ) == 0:
            return True
        else:
            return False
    
    def __ne__( self, T ):
        if self.__cmp__( T ) != 0:
            return True
        else:
            return False
    
    def __le__( self, T ):
        if self.__cmp__( T ) <= 0:
            return True
        else:
            return False
    
    def __ge__( self, T ):
        if self.__cmp__( T ) >= 0:
            return True
        else:
            return False
    
    def __lt__( self, T ):
        if self.__cmp__( T ) < 0:
            return True
        else:
            return False
    
    def __gt__( self, T ):
        if self.__cmp__( T ) > 0:
            return True
        else:
            return False

def cmpQPDateTime( dt1, dt2, epsilon = 0.0 ):
    return mx.DateTime.cmp( dt1.datetime, dt2.datetime, epsilon )

def diffQPDateTime( dt1, dt2 ):
    """
    return difference between two QPDateTime objects
    return value is of type mx.DateTime.DateTimeDelta
    """
    return dt1.datetime - dt2.datetime