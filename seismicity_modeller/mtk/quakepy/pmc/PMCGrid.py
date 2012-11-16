# -*- coding: iso-8859-1 -*-
#
# quakepy/pmc/PMCGrid.py
# $Id: PMCGrid.py 166 2009-03-24 17:49:56Z fab $
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

__version__  = '$Id: PMCGrid.py 166 2009-03-24 17:49:56Z fab $'
__revision__ = '$Revision: 166 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import copy
import cPickle
import math
import numpy

import pyRXP
from   xml.sax import saxutils

import matplotlib
matplotlib.use('PS')
from pylab import *

from mx.DateTime import Date

sys.path.append('../..')
sys.path.append('..')

from PMCData import *
from PMCInventory import *

from QPCore import *
from QPGrid import *

class PMCGrid( QPGrid ):
    """
    QuakePy: PMCGrid

    PMCGrid adds a PMCData object to each grid cell of a QPGrid
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                    QPElement( 'PMCData', 'PMCData', 'element', PMCData, 'complex',
                               '/QPGrid/grid/depthLayer/cell', Cell ),
                    QPElement( 'stations', 'PMCStations', 'element', PMCStationList, 'complex',
                               '/QPGrid', QPGrid ),
                    QPElement( 'annotation', 'annotation', 'element', QPAnnotation, 'complex',
                               '/QPGrid', QPGrid )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, **kwargs ):
        super( PMCGrid, self ).__init__( None, **kwargs )
        self.elements.extend( self.addElements )
                        
        self._initMultipleElements()

        
