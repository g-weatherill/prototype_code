#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# quakepy/pmc/scsn/grid_job_template.py
# $Id: grid_job_template.py 78 2008-06-17 16:20:56Z fab $
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

__version__  = '$Id: grid_job_template.py 78 2008-06-17 16:20:56Z fab $'
__revision__ = '$Revision: 78 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import os
import re
import random
import datetime

import pickle
import cPickle

#import matplotlib
#from matplotlib import rcParams
#from matplotlib.dates import MO
#from matplotlib.dates import YearLocator, MonthLocator, DayLocator, WeekdayLocator, DateFormatter
#matplotlib.use('PS')

#from pylab import *

#import Polygon

from   mx.DateTime     import Date, DateTime
from   mx.DateTime.ISO import ParseDateTimeUTC

sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')


from QPCore  import *
from QPUtils  import *
from QPCatalog import *

from pmc.PMC      import *
from pmc.PMC_SCSN import *
from pmc.PMCData  import *

from QPGrid           import *

"""
template for computing probability grids on cluster
"""

useDate       = DATE_USED
useDistStyle  = 3

maggrid       = frange( 0.0, 4.0, 0.1 )
mp_prob_list  = [ 0.9, 0.95, 0.99, 0.999 ]

# set geographical region for completeness evaluation
geobox = { 'lonmin': -122.0, 'lonmax': -113.5, 
           'latmin': 31.5, 'latmax': 38.0,
           'londelta': 0.1, 'latdelta': 0.1,
           'depthmin': 7.5, 'depthmax': 7.5 }
               
basepath = '/home/geovault-00/fabian/prog/pyprog/quakepy/pmc/'

distro_xml_path   = basepath + 'data/distro/'
grid_xml_path     = basepath + 'data/grid/'
combi_path        = basepath + 'data/'

pmc_xml_in        = distro_xml_path + 'scsn.distro.xml'
combi_pickle      = combi_path + 'combinations-3-375.numpy.pickle'

print " ----- get PMC from XML -----  "

print " load Combinations from pickle"
combi = unpickleObj( combi_pickle )
    
pmc = PMC_SCSN( combinations = combi )

print " importing PMC from XML file: ", pmc_xml_in
pmc.readXML( getQPDataSource( pmc_xml_in ) )

print " loaded PMC from XML: %s stations" % ( len( pmc.stations ) )

print " ----- compute QPGrid -----  "

# compute grid

g = QPGrid()
g.annotation = pmc.annotation

g.setupBox( geobox['lonmin'], geobox['lonmax'], 
            geobox['latmin'], geobox['latmax'], 
            geobox['depthmin'], geobox['depthmax'], 
            geobox['londelta'], geobox['latdelta'] )

mapdate = Date( useDate[0], useDate[1], useDate[2] )
mapdate_str = mxDateTime2ISO( [ mapdate ] )

# set grid annotation
g.annotation.setProperty( date = utc(),
                          starttime = mapdate,
                          endtime = mapdate,
                          latmin = geobox['latmin'],
                          latmax = geobox['latmax'],
                          lonmin = geobox['lonmin'],
                          lonmax = geobox['lonmax']  )
                              
print " compute probability map for Date: %s" % ( mapdate_str )
pmc.getProbabilityMap( g, maggrid, mapdate, mp_prob_list, verbose = False )

# add station data to grid - delete not required fields
g.stations = []
for curr_sta in pmc.stations:
    del curr_sta.channels
    del curr_sta.distribution
    g.stations.append( curr_sta )
        
grid_file_out = grid_xml_path + 'scsn.' + mapdate_str + '.grid.xml'
print " write computed grid file ", grid_file_out
g.writeXML( writeQPData( grid_file_out ) )
