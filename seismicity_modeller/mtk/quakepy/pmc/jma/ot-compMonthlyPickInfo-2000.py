#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

#
# quakepy/pmc/jma/ot-compMonthlyPickInfo-2000.py
# $Id: ot-compMonthlyPickInfo-2000.py 98 2008-07-30 07:12:44Z fab $
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

__version__  = '$Id: ot-compMonthlyPickInfo-2000.py 98 2008-07-30 07:12:44Z fab $'
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
    
    time_year   = 2000 # 2001, ..., 2008

    basepath    = '/home/fab/prog/pyprog/quakepy/pmc/jma/'
    catdir      = basepath + 'data/catalog-bz2/'
    pickInfoDir = basepath + 'data/pickInfo-monthly/'

    stationfile      = basepath + 'data/jma_station_list_with_network_code.dat'

    # catfiles       =  ( catdir + '200007_A.deck.Z.bz2', )
    catfiles_full = sort( glob.glob( catdir + '*.Z.bz2' ) )

    # get catalog decks for monthly chunk
    print " ----- processing year %s" % time_year

    catfiles = []
    for idx in xrange( 12 ):
        catfiles.append( { 'name': '', 'files': [] } )

    # get catfiles for time period
    for curr_catfile in catfiles_full:

        # get year and month
        # print " search in %s" % curr_catfile

        datecode = re.match( r'^.+(\d{6})_[ABC].deck.Z.bz2$', curr_catfile )
        yearmonth_str = str( datecode.group(1) )
        year     = int( yearmonth_str[0:4] )
        month    = int( yearmonth_str[4:6] )

        if ( year == time_year ):
            catfiles[month-1]['name'] = yearmonth_str
            catfiles[month-1]['files'].append( curr_catfile )

            print " added catalog file %s to current catalog" % curr_catfile

    # delete missing months from catfiles vector
    for idx in reversed( xrange( len( catfiles ) ) ):
        if catfiles[idx]['name'] == '':
            del catfiles[idx]

    # loop over catalog files in months
    print " ----- looping over %s monthly chunks" % len( catfiles )
    
    for curr_month in xrange( len( catfiles ) ):

        # new inventory object
        print " importing PMC inventory file: ", stationfile
        pmc = PMC_JMA( stationfile, 0, smooth = True, encoding = 'html' )
        print " full inventory has ", len(pmc.stations), " stations"
        
        print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"
        pmc.preprocessInventory( [ DateTime(2000, 1, 1, 0, 0, 0), DateTime(2008, 4, 1, 0, 0, 0) ] )
        print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"
        
        for catfile in catfiles[curr_month]['files']:
            print "   ----- importing catalog chunk from file: ", catfile

            qpc = QPCatalog( idstyle='numeric' )
            qpc.importJMADeck( catfile, compression='bz2', minimumDataset=True )
            print "   catalog chunk has %s events" % qpc.size()

            pmc.preprocessCatalog( qpc, [ DateTime(2000, 1, 1, 0, 0, 0), DateTime(2008, 4, 1, 0, 0, 0) ] )
            print "   after catalog preprocess: catalog chunk has %s events" % qpc.size()

            print "   assigning picks ... "
            pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

            del qpc

        print " now writing station data"
        for curr_sta_idx in reversed( range( len( pmc.stations ) ) ):

            curr_sta = pmc.stations[curr_sta_idx]

            print "   ----- processing station %s: %s %s" % ( curr_sta_idx, curr_sta.locationCode, curr_sta.stationCode )

            # see if directory is already there
            sta_file_dir = pickInfoDir + curr_sta.locationCode + '.' + curr_sta.stationCode + '/'
            if not os.path.isdir( sta_file_dir ):
                print "   creating path %s" % sta_file_dir
                os.mkdir( sta_file_dir )

            sta_filename_xml = str(curr_sta_idx) + '.' + curr_sta.locationCode + '.' + curr_sta.stationCode + '.' + catfiles[curr_month]['name'] + '.pickInfo.xml.bz2'
            sta_file_xml     = sta_file_dir + sta_filename_xml

            print "   ----- writing to XML"
            fh = writeQPData( sta_file_xml, compression='bz2' )
            fh.write( '<?xml version="1.0" encoding="utf-8"?>' )
            fh.write( '<PMC>' )
            curr_sta.toXML( 'PMCStation', fh )
            fh.write( '</PMC>' )
            fh.close()

            # delete station just processed
            del pmc.stations[curr_sta_idx]

main()