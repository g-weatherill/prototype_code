# -*- coding: utf-8 -*-
#
# quakepy/cumuldist.py
# $Id: cumuldist.py 277 2011-03-02 11:20:15Z fab $
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


__version__  = '$Id: cumuldist.py 277 2011-03-02 11:20:15Z fab $'
__revision__ = '$Revision: 277 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import qpplot

class CumulativeDistribution( object ):
    def __init__( self, evpar ):
        
        # over events and preferred origins
        curr_cd = []
        
        for curr_ev in evpar.event:
            
            curr_ori_time = curr_ev.getPreferredOrigin().time.value.datetime
            
            # append timestamp to list
            curr_cd.append( curr_ori_time )
        
        # sort list and add consecutive number
        self.cd = [ [ curr_cd_val.strftime('%Y-%m-%dT%H:%M:%S'), 
            str(curr_cd_idx+1) ] for curr_cd_idx, curr_cd_val in enumerate( 
            sorted( curr_cd ) ) ]
        del curr_cd
        
        
    def plot( self, imgfile=None, **kwargs ):
        return qpplot.QPPlot().plot_vs_date( imgfile, 
            [ curr_data[0] for curr_data in self.cd ], 
            [ curr_data[1] for curr_data in self.cd ],
            **kwargs )
    
def main():
    pass
    
if __name__ == '__main__':
    main()
