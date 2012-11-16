# -*- coding: iso-8859-1 -*-
#
# quakepy/pmc/PMCData.py
# $Id: PMCData.py 159 2009-02-17 13:34:05Z fab $
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

__version__  = '$Id: PMCData.py 159 2009-02-17 13:34:05Z fab $'
__revision__ = '$Revision: 159 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys

sys.path.append('../..')
sys.path.append('..')

from QPCore import *


class PMCDataMP( QPObject ):
    """
    QuakePy: PMCDataMP
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                      QPElement( 'probability', 'probability', 'attribute', float, 'basic' ),
                      QPElement( 'mag', 'mag', 'cdata', float, 'basic' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, mag = None, prob = None, **kwargs ):
        """
        <PMCData>
            <probability magnitude="1.3">0.4456</probability>
            <probability magnitude="1.4">0.4456</probability>
            <probability magnitude="1.5">0.4456</probability>
            <mp probability="0.99">2.4</mp>
            <mp probability="0.9">2.0</mp>
        </PMCData>
        """
        super( PMCDataMP, self ).__init__( **kwargs )
        self.elements.extend( self.addElements )
        self._initMultipleElements()

        self.probability = prob
        self.mag         = mag

## ---------------------------------------------------------------------------

class PMCDataProbability( QPObject ):
    """
    QuakePy: PMCDataProbability
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                      QPElement( 'magnitude', 'magnitude', 'attribute', float, 'basic' ),
                      QPElement( 'prob', 'prob', 'cdata', float, 'basic' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, prob = None, mag = None, **kwargs ):
        """
        <PMCData>
            <probability magnitude="1.3">0.4456</probability>
            <probability magnitude="1.4">0.4456</probability>
            <probability magnitude="1.5">0.4456</probability>
            <mp probability="0.99">2.4</mp>
            <mp probability="0.9">2.0</mp>
        </PMCData>
        """
        super( PMCDataProbability, self ).__init__( **kwargs )
        self.elements.extend( self.addElements )
        self._initMultipleElements()

        self.magnitude = mag
        self.prob      = prob

## ---------------------------------------------------------------------------

class PMCData( QPObject ):
    """
    QuakePy: PMCData
    """

    # <!-- UML2Py start -->
    addElements = QPElementList( (
                      QPElement( 'probability', 'probability', 'element', PMCDataProbability, 'multiple' ),
                      QPElement( 'mp', 'mp', 'element', PMCDataMP, 'multiple' )
                  ) )
    # <!-- UML2Py end -->
    
    def __init__( self, probability = None, mp = None, **kwargs ):
        """
        <PMCData>
            <probability magnitude="1.3">0.4456</probability>
            <probability magnitude="1.4">0.4456</probability>
            <probability magnitude="1.5">0.4456</probability>
            <mp probability="0.99">2.4</mp>
            <mp probability="0.9">2.0</mp>
        </PMCData>
        """
        super( PMCData, self ).__init__( **kwargs )
        self.elements.extend( self.addElements )
        self._initMultipleElements()

        if probability is not None:
            self.probability = probability

        if mp is not None:
            self.mp = mp

## ---------------------------------------------------------------------------








