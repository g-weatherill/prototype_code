#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

#
# quakepy/pmc/jma/ot-compYearlyDistros-all-01.py
# $Id: ot-compYearlyDistros-all-01.py 98 2008-07-30 07:12:44Z fab $
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

__version__  = '$Id: ot-compYearlyDistros-all-01.py 98 2008-07-30 07:12:44Z fab $'
__revision__ = '$Revision: 98 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import sys
import os
import glob
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

from mx.DateTime     import Date
from mx.DateTime     import DateTime
from mx.DateTime.ISO import ParseDateTimeUTC

sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')

from QPCore    import *
from QPUtils   import *
from QPCatalog import *

from pmc.PMC      import *
from pmc.PMC_JMA  import *
from pmc.PMCData  import *

from QPGrid           import *

#from qpplot           import *
#from qpseismicityplot import *

def main():
    """
    compute pickInfo and distros for all JMA stations on dynamic.usc.edu
    in chunks of ~50 stations, starting at station index 00
    
    create PMC_JMA for a chunk of stations, distStyle 0 so that we have no <distro> element in XML
    loop over catalog chunks for a sequence of time periods
    read in catalog chunk from JMA Deck file

    - read complete station file into PMC object, with added region/network information
    - delete stations which are not in desired range from PMC object
    - loop over remaining stations in PMC object
      - assign picks (= fill pickInfo for station)
      - write station as XML file NO.LOC.STA.pickInfo.xml
      - write XML station file name to stationlist file
    """
    
    station_chunk = ( 0, 199 ) # ( 200, 399 ), ..., ( 1400, 1659 )

    time_intervals = ( ( '200001', '200012' ),
                       ( '200005', '200104' ),
                       ( '200009', '200108' ),

                       ( '200101', '200112' ),
                       ( '200105', '200204' ),
                       ( '200109', '200208' ),

                       ( '200201', '200212' ),
                       ( '200205', '200304' ),
                       ( '200209', '200308' ),

                       ( '200301', '200312' ),
                       ( '200305', '200404' ),
                       ( '200309', '200408' ),

                       ( '200401', '200412' ),
                       ( '200405', '200504' ),
                       ( '200409', '200508' ),
                       
                       ( '200501', '200512' ),
                       ( '200505', '200604' ),
                       ( '200509', '200608' ),

                       ( '200601', '200612' ),
                       ( '200605', '200704' ),
                       ( '200609', '200708' ),

                       ( '200701', '200712' ) )

    files_per_interval = 12
    
    useDistStyle  = 5
    
    basepath = '/home/fab/prog/pyprog/quakepy/pmc/jma/'

    pickInfoDir = basepath + 'data/pickInfo-monthly/'
    distroDir   = basepath + 'data/distro-yearly/'

    stationfile = basepath + 'data/jma_station_list_with_network_code.dat'
    
    print " importing PMC inventory file: ", stationfile
    pmc = PMC_JMA( stationfile )
    print " full inventory has ", len(pmc.stations), " stations"
    
    # select only stations from desired chunk
    print " selecting desired stations from inventory: %s to %s" % ( station_chunk[0], station_chunk[1] )

    for sta_idx in reversed( xrange( station_chunk[1]+1, len(pmc.stations) ) ):
        del pmc.stations[sta_idx]

    for sta_idx in reversed( xrange( 0, station_chunk[0] ) ):
        del pmc.stations[sta_idx]    
        
    print " after selecting stations: inventory has ", len(pmc.stations), " stations"
    
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"
    pmc.preprocessInventory( [ DateTime(2000, 1, 1, 0, 0, 0), DateTime(2008, 4, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    for curr_sta_idx in reversed( range( len( pmc.stations ) ) ):

        sta_number = station_chunk[0] + curr_sta_idx
        curr_sta   = pmc.stations[curr_sta_idx]

        print " processing station %s: %s %s" % ( sta_number, curr_sta.locationCode, curr_sta.stationCode )
        
        sta_dirname = curr_sta.locationCode + '.' + curr_sta.stationCode + '/'

        # get list of all filenames
        pickfiles = sort( glob.glob( pickInfoDir + sta_dirname + '*.pickInfo.xml.bz2' ) )
          
        # over time intervals
        for curr_ti in time_intervals:

            print " computing distro for time interval %s, %s" % ( curr_ti[0], curr_ti[1] )

            curr_pmc = PMC_JMA()
            
            # get monthly distro files that we need
            chunk_ctr = 0
            for pickfile in pickfiles:

                # print " checking file %s for chunk %s" % ( pickfile, curr_chunk )

                datecode = re.match( r'^.+(\d{6}).pickInfo.xml.bz2$', pickfile )
                yearmonth_str = str( datecode.group(1) )

                if ( not (    int( yearmonth_str ) >= int( curr_ti[0] )
                          and int( yearmonth_str ) <= int( curr_ti[1] ) ) ):

                    # if file is not in desired interval, try next one
                    continue
                else:

                    # gotcha, file found
                    chunk_ctr += 1
            
                    # if it's first chunk, create PMC object for curr station
                    if chunk_ctr == 1:
                        curr_pmc.readXML( getQPDataSource( pickfile, compression='bz2' ) )
                    else:
                        
                        # read temporary pickInfo
                        pmc_temp = PMC_JMA()
                        pmc_temp.readXML( getQPDataSource( pickfile, compression='bz2' ) )

                        # copy over entries to master PMC, if not empty
                        if ( pmc_temp.stations[0].distribution.pickInfo is not None ):

                            # concatenate with existing distro if it is not None
                            if ( curr_pmc.stations[0].distribution.pickInfo is not None ):

                            
                                pickInfo_temp = numpy.concatenate( ( curr_pmc.stations[0].distribution.pickInfo,
                                                                     pmc_temp.stations[0].distribution.pickInfo ) )
                            else:

                                # existing distro is None, replace with current one
                                pickInfo_temp = pmc_temp.stations[0].distribution.pickInfo

                            curr_pmc.stations[0].distribution.pickInfo = pickInfo_temp
                            

                        del pmc_temp

            # no files found, no distro for that time interval
            if chunk_ctr == 0:
                print " NO FILES FOUND, CONTINUE ..."
                continue
                
            ## compute distro
            if curr_pmc.stations[0].distribution.pickInfo  == None:
                print " NO PICK INFO, CONTINUE ..."
                continue

            if chunk_ctr != files_per_interval:
                print " NOTE: MISSING FILES, only %s pickInfo files for this interval" % chunk_ctr
                
            curr_pmc.stations[0].distribution.restoreDistStyle( useDistStyle )
            curr_pmc.stations[0].distribution.setSmooth( True )

            print " now calling fillup ..."
            curr_pmc.stations[0].fillup( pmc._calcMagnitudeFromDistance )
            
            # delete pick info
            del curr_pmc.stations[0].distribution.pickInfo
            
            # get dir
            distro_dir = distroDir + curr_ti[0] + '-' + curr_ti[1] + '/'

            if not os.path.isdir( distro_dir ):
                os.mkdir( distro_dir )
            
            # get filename for probdistro XML
            distro_filename_xml = str(sta_number) + '.' + curr_sta.locationCode + '.' + curr_sta.stationCode + '.' + curr_ti[0] + '-' + curr_ti[1] + '.distro.xml.bz2'
            
            distro_file_out = distro_dir + distro_filename_xml

            # write PMC as XML
            print " writing distro XML file: ", distro_file_out
            fh = writeQPData( distro_file_out, compression='bz2' )
            curr_pmc.writeXML( fh )
            fh.close()

            # delete station just processed 
            del curr_pmc


main()