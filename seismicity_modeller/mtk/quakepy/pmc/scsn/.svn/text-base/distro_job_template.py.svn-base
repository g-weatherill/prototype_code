#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

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

import sys
import os
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

from   mx.DateTime     import DateTime
from   mx.DateTime.ISO import ParseDateTimeUTC

sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')

from recipes import *

from QPCore  import *
from QPUtils  import *
from QPCatalog import *

from pmc.PMC      import *
from pmc.PMC_SCSN import *
from pmc.PMCData  import *

from QPGrid           import *

"""
template for computing distros on cluster
1) load file list of pickInfo XML files
2) load PMC for one station with pickInfo from XML file / get filename from file list (idx.NC.STA.pickInfo.xml)
3) reset distStyle to 3
4) call fillup on station
5) write PMC (with only one station) as XML: idx.NC.STA.pickInfo.distro.xml
6) delete pickInfo
7) write PMC (with only one station) as XML: idx.NC.STA.distro.xml
"""

useStationIdx = STATION_INDEX
useDistStyle  = 3

basepath = '/home/geovault-00/fabian/prog/pyprog/quakepy/pmc/'

pickinfo_path     = basepath + 'data/pickInfo/'
pickinfo_filelist = pickinfo_path + 'pmc.sc.pickInfo.filelist.dat'

pickinfo_distro_xml_path = basepath + 'data/distro/'
distro_xml_path          = basepath + 'data/distro/'

# read pickinfo filelist
print " opening %s ..." % ( pickinfo_filelist )
pickinfo_file_arr = getQPDataSource( pickinfo_filelist ).read().split()
print " i could find %s pickinfo files" % ( len(pickinfo_file_arr) )

pmcfile_xml = pickinfo_path + pickinfo_file_arr[useStationIdx]

pmc = PMC_SCSN()
print " importing PMC from XML file: ", pmcfile_xml
pmc.readXML( getQPDataSource( pmcfile_xml ) )
print " read PMC inventory with ", len(pmc.stations), " stations"

if len(pmc.stations) > 1:
    print " TROUBLE: i have too many stations!"

# select station
use_sta = pmc.stations[0]

use_sta.distribution.restoreDistStyle( useDistStyle )
use_sta.distribution.setSmooth( True )

print " now calling fillup ..."
use_sta.fillup( pmc._calcMagnitudeFromDistance )

# get filename
pmc_xml_out = pickinfo_distro_xml_path + str(useStationIdx) + '.' + use_sta.networkCode + '.' + use_sta.stationCode + '.pickInfo.distro.xml'

# write PMC as XML
print " writing XML file: ", pmc_xml_out
fh = writeQPData( pmc_xml_out )
pmc.writeXML( fh )
fh.close()

# delete pick info
del use_sta.distribution.pickInfo

# get filename
pmc_xml_out = distro_xml_path + str(useStationIdx) + '.' + use_sta.networkCode + '.' + use_sta.stationCode + '.distro.xml'

# write PMC as XML
print " writing XML file: ", pmc_xml_out
fh = writeQPData( pmc_xml_out )
pmc.writeXML( fh )
fh.close()
