# -*- coding: utf-8 -*-
#
# quakepy/pmc/PMCFactory.py
# $Id$
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

__version__  = '$Id$'
__revision__ = '$Revision$'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys

sys.path.append('..')

import PMC_SCSN
import PMC_JMA
import PMC_NZ
import PMC_NCSN


class PMCFactory( object ):
    """
    QuakePy: PMCFactory
    factory to create PMC instance, depending on network
    """
    
    def createPMC( self, network, stationlist_filename=None, dist_style=None, 
        useSubnetworks=False, **kwargs ):

        pmc = None
        
        if network == 'SC':
            pmc = PMC_SCSN.PMC_SCSN( stationlist_filename, dist_style, 
                useSubnetworks=useSubnetworks, **kwargs )
        elif network == 'JP':
            pmc = PMC_JMA.PMC_JMA( stationlist_filename, dist_style, 
                useSubnetworks=useSubnetworks, **kwargs )
        elif network == 'NZ':
            pmc = PMC_NZ.PMC_NZ( stationlist_filename, dist_style, 
                useSubnetworks=useSubnetworks, **kwargs )
        elif network == 'NC':
            pmc = PMC_NCSN.PMC_NCSN( stationlist_filename, dist_style, 
                useSubnetworks=useSubnetworks, **kwargs )
        else:
            error_str = "PMCFactory.object() - unknown network identifier: %s" % ( network )
            raise ValueError, error_str

        return pmc
