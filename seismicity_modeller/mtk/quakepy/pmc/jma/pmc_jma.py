#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# quakepy/pmc/jma/pmc_jma.py
# $Id: pmc_jma.py 174 2009-03-26 23:32:57Z fab $
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

__version__  = '$Id: pmc_jma.py 174 2009-03-26 23:32:57Z fab $'
__revision__ = '$Revision: 174 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
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

def task_1():
    """
    compute pickInfo for all JMA stations
    
    create PMC_JMA with full station list, distStyle 0 so that we have no <distro> element in XML
    loop over catalog chunks, read in catalog chunk from JMA Deck file

    loop over stations:
     - assign picks (= fill pickInfo for station)
     - write station as XML snippet NO.NC.STA.pickInfo.xml
     - write XML station file name to stationlist file

    write PMC object with pickInfo as XML

    NOTE: it became clear that 16 GB memory is not sufficient to hold
    the PMC object with pickInfo data for all ~1600 stations and full catalog

    for a computation of pickInfo for all stations and the time span
    of 2007-200803 see task_2
    
    we compute the pickInfos separately for chunks of ~200 stations
    see task_5
    """

    catdir         = './data/catalog-bz2/'

    # catfiles       =  ( catdir + '200007_A.deck.Z.bz2', )
    catfiles = sort( glob.glob( catdir + '*.Z.bz2' ) )

    stationfile   = './data/jma_station_list.dat'
    
    pmc_xml_file     = './data/pickInfo/pmc.jma.pickInfo.xml.bz2'
    stationlist_file = './data/pickInfo/pmc.jma.pickInfo.filelist.dat'

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_JMA( stationfile, 0, smooth = True, encoding = 'html' )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"

    pmc.preprocessInventory( [ DateTime(2000, 1, 1, 0, 0, 0), DateTime(2008, 4, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    # loop over catalog chunk pickles
    for catfile in catfiles:
        print " ----- importing catalog chunk from file: ", catfile

        qpc = QPCatalog()
        qpc.importJMADeck( catfile, compression='bz2' )
        print " catalog has %s events" % qpc.size()

        pmc.preprocessCatalog( qpc, [ DateTime(2000, 1, 1, 0, 0, 0), DateTime(2008, 4, 1, 0, 0, 0) ] )
        print " after catalog preprocess: catalog has %s events" %  qpc.size()

        print " assigning picks ... "
        pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

        del qpc

    print " after assign picks: inventory has ", len(pmc.stations), " stations"

    # over stations - write station with pickInfo separately as XML
    print " +++++ NOW WRITING STATION DATA AS XML +++++ "

    # open stationlist file
    flist = writeQPData( stationlist_file )
    
    for curr_sta_idx, curr_sta in enumerate( pmc.stations ):

        sta_filename_xml = str(curr_sta_idx) + '.' + curr_sta.networkCode + '.' + curr_sta.stationCode + '.pickInfo.xml.bz2'
        sta_file_xml     = './data/pickInfo/' + sta_filename_xml
    
        print " processing station %s: %s %s" % ( curr_sta_idx, curr_sta.networkCode, curr_sta.stationCode )

        print "  writing to XML"
        fh = writeQPData( sta_file_xml, compression='bz2' )
        fh.write( '<?xml version="1.0" encoding="utf-8"?>' )
        fh.write( '<PMC>' )
        curr_sta.toXML( 'PMCStation', fh )
        fh.write( '</PMC>' )
        fh.close()

        # write station filename to list file
        flist.write( ''.join( ( sta_filename_xml, '\n' ) ) )

    flist.close()
        
    # write PMC as XML
    print "  writing PMC to XML"
    pmc.writeXML( writeQPData( pmc_xml_file, compression='bz2' ) )

def task_2():
    """
    compute pickInfo for all JMA stations, for the time span of 2007-200803
    
    create PMC_JMA with full station list, distStyle 0 so that we have no <distro> element in XML
    loop over subset of catalog chunks, read in catalog chunk from JMA Deck file

    loop over stations:
     - assign picks (= fill pickInfo for station)
     - write station as XML snippet NO.NC.STA.pickInfo.xml
     - write XML station file name to stationlist file

    write PMC object with pickInfo as XML
    """

    catdir   = './data/catalog-bz2/'

    catfiles =  ( catdir + '200701_A.deck.Z.bz2',
                  catdir + '200701_B.deck.Z.bz2',
                  catdir + '200701_C.deck.Z.bz2',
                  catdir + '200702_A.deck.Z.bz2',
                  catdir + '200702_B.deck.Z.bz2',
                  catdir + '200702_C.deck.Z.bz2',
                  catdir + '200703_A.deck.Z.bz2',
                  catdir + '200703_B.deck.Z.bz2',
                  catdir + '200703_C.deck.Z.bz2',
                  catdir + '200704_A.deck.Z.bz2',
                  catdir + '200704_B.deck.Z.bz2',
                  catdir + '200704_C.deck.Z.bz2',
                  catdir + '200705_A.deck.Z.bz2',
                  catdir + '200705_B.deck.Z.bz2',
                  catdir + '200705_C.deck.Z.bz2',
                  catdir + '200706_A.deck.Z.bz2',
                  catdir + '200706_B.deck.Z.bz2',
                  catdir + '200706_C.deck.Z.bz2',
                  catdir + '200707_A.deck.Z.bz2',
                  catdir + '200707_B.deck.Z.bz2',
                  catdir + '200707_C.deck.Z.bz2',
                  catdir + '200708_A.deck.Z.bz2',
                  catdir + '200708_B.deck.Z.bz2',
                  catdir + '200708_C.deck.Z.bz2',
                  catdir + '200709_A.deck.Z.bz2',
                  catdir + '200709_B.deck.Z.bz2',
                  catdir + '200709_C.deck.Z.bz2',
                  catdir + '200710_A.deck.Z.bz2',
                  catdir + '200710_B.deck.Z.bz2',
                  catdir + '200710_C.deck.Z.bz2',
                  catdir + '200711_A.deck.Z.bz2',
                  catdir + '200711_B.deck.Z.bz2',
                  catdir + '200711_C.deck.Z.bz2',
                  catdir + '200712_A.deck.Z.bz2',
                  catdir + '200712_B.deck.Z.bz2',
                  catdir + '200712_C.deck.Z.bz2',
                  catdir + '200801_A.deck.Z.bz2',
                  catdir + '200801_B.deck.Z.bz2',
                  catdir + '200801_C.deck.Z.bz2',
                  catdir + '200802_A.deck.Z.bz2',
                  catdir + '200802_B.deck.Z.bz2',
                  catdir + '200802_C.deck.Z.bz2',
                  catdir + '200803_A.deck.Z.bz2',
                  catdir + '200803_B.deck.Z.bz2',
                  catdir + '200803_C.deck.Z.bz2' )

    #catfiles = sort( glob.glob( catdir + '*.Z.bz2' ) )

    stationfile   = './data/jma_station_list.dat'
    
    pmc_xml_file     = './data/pickInfo/pmc.jma.2007-200803.pickInfo.xml.bz2'
    stationlist_file = './data/pickInfo/pmc.jma.2007-200803.pickInfo.filelist.dat'

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_JMA( stationfile, 0, smooth = True, encoding = 'html' )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"

    pmc.preprocessInventory( [ DateTime(2000, 1, 1, 0, 0, 0), DateTime(2008, 4, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    # loop over catalog chunk pickles
    for catfile in catfiles:
        print " ----- importing catalog chunk from file: ", catfile

        qpc = QPCatalog()
        qpc.importJMADeck( catfile, compression='bz2' )
        print " catalog has %s events" % qpc.size()

        pmc.preprocessCatalog( qpc, [ DateTime(2000, 1, 1, 0, 0, 0), DateTime(2008, 4, 1, 0, 0, 0) ] )
        print " after catalog preprocess: catalog has %s events" %  qpc.size()

        print " assigning picks ... "
        pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

        del qpc

    print " after assign picks: inventory has ", len(pmc.stations), " stations"

    # over stations - write station with pickInfo separately as XML
    print " +++++ NOW WRITING STATION DATA AS XML +++++ "

    # open stationlist file
    flist = writeQPData( stationlist_file )
    
    for curr_sta_idx, curr_sta in enumerate( pmc.stations ):

        sta_filename_xml = str(curr_sta_idx) + '.' + curr_sta.networkCode + '.' + curr_sta.stationCode + '.2007-200803.pickInfo.xml.bz2'
        sta_file_xml     = './data/pickInfo/' + sta_filename_xml
    
        print " processing station %s: %s %s" % ( curr_sta_idx, curr_sta.networkCode, curr_sta.stationCode )

        print "  writing to XML"
        fh = writeQPData( sta_file_xml, compression='bz2' )
        fh.write( '<?xml version="1.0" encoding="utf-8"?>' )
        fh.write( '<PMC>' )
        curr_sta.toXML( 'PMCStation', fh )
        fh.write( '</PMC>' )
        fh.close()

        # write station filename to list file
        flist.write( ''.join( ( sta_filename_xml, '\n' ) ) )

    flist.close()
        
    # write PMC as XML
    print "  writing PMC to XML"
    pmc.writeXML( writeQPData( pmc_xml_file, compression='bz2' ) )

def task_3():
    """
    compute distros for all stations
    catalog from 2007 to 2008-03
    
    1) load file list of pickInfo XML files (pickInfos for 2007-200803 have been computed previously, task_2)
    2) loop over pickInfo XML files (= over stations)
    3) load PMC for one station with pickInfo from XML file / get filename from file list (idx.NC.STA.pickInfo.xml.bz2)
    4) reset distStyle to 4
    5) call fillup on station (= compute distro)
    6) delete pickInfo
    7) write PMC as XML per station: idx.NC.STA.distro.xml.bz2

    NOTE: we split up the list of pickInfo files into 8 parts and started
    8 processes for computing the distros in parallel

    -> compdistro-2007-200803-NNN.py (8 chunks) on othello.ethz.ch
    -> 2007-200803.pickInfo.files.001.dat, ..., 2007-200803.pickInfo.files.008.dat
    """

    useDistStyle  = 4

    pickinfo_filelist = './data/pickInfo/pmc.jma.2007-200803.pickInfo.filelist.dat'
    pickinfo_path     = './data/pickInfo/'

    pickinfo_distro_xml_path = './data/distro/'
    distro_xml_path          = './data/distro/'
    distro_eps_path          = './data/distro-eps/'

    colortable_file = '../scsn/data/colortables/rastafari.cpt'

    # read pickinfo filelist
    pickinfo_file_arr = getQPDataSource( pickinfo_filelist ).read().split()
    print " i could find %s pickinfo files" % len(pickinfo_file_arr)

    staCtr = 0
    for curr_pickinfo_file in pickinfo_file_arr:

        staCtr += 1
        
        pmcfile_xml = pickinfo_path + curr_pickinfo_file

        pmc = PMC_JMA()
        
        print " importing PMC from XML file: ", pmcfile_xml
        pmc.readXML( getQPDataSource( pmcfile_xml, compression='bz2' ) )
        print " read PMC inventory with ", len(pmc.stations), " stations"

        if len(pmc.stations) > 1:
            print " TROUBLE: i have too many stations!"

        # select station
        use_sta = pmc.stations[0]

        if use_sta.distribution.pickInfo  == None:
            continue
        
        use_sta.distribution.restoreDistStyle( useDistStyle )
        use_sta.distribution.setSmooth( True )

        print " now calling fillup ..."
        use_sta.fillup( pmc._calcMagnitudeFromDistance )

        # get filename
        pmc_xml_out = pickinfo_distro_xml_path + str(staCtr) + '.' + use_sta.networkCode + '.' + use_sta.stationCode + '.2007-200803.distro.xml.bz2'

        # delete pick info
        del use_sta.distribution.pickInfo

        # write PMC as XML
        print " writing XML file: ", pmc_xml_out
        fh = writeQPData( pmc_xml_out, compression='bz2' )
        pmc.writeXML( fh )
        fh.close()

        # write probdistro
        # get filename
        pmc_probdist_out = distro_eps_path + str(staCtr) + '.' + use_sta.networkCode + '.' + use_sta.stationCode + '.2007-200803.distro.eps'
        
        pmc.plotDistribution( pmc.stations[0], pmc_probdist_out, colortable_file )
    
def task_4():
    """
    compute combination list class
    use NumPy version
    class has changed from revision 57, old pickles won't work any more
    
    Southern California: M=3, N=375
    Northern California: M=3, N=630
    JMA (Japan)        : M=3, N=100
    """
    
    max_N = 200   # 100, 200, 375, 630
    max_M = 3

    pickle_file_py   = 'combinations-' + str(max_M) + '-' + str(max_N) + '.py.pickle'
    pickle_file_numpy = 'combinations-' + str(max_M) + '-' + str(max_N) + '.numpy.pickle'

    # use numpy - MEMORY ERROR ON 32 bit 4GB machine (4 / 3.25 GB) FOR 3 / 630 when creating numpy.ones for m=3 !!
    print " NumPy: compute combinations for M=%s, N=%s" % ( max_M, max_N )
    c = CombinationsNumPy( max_M, max_N )

    print " pickling combinations object"
    pickleObj( c, pickle_file_numpy )
    
    del c
    
    # use native python
    print " Native Python: compute combinations for M=%s, N=%s" % ( max_M, max_N )
    c = Combinations( max_M, max_N )

    print " pickling combinations object"
    pickleObj( c, pickle_file_py )

    del c

def task_5():
    """
    compute pickInfo for all JMA stations, in chunks of ~200 stations
    
    create PMC_JMA for a chunk of stations, distStyle 0 so that we have no <distro> element in XML
    loop over catalog chunks, read in catalog chunk from JMA Deck file

    - read complete station file into PMC object, with added region/network information
    - delete stations which are not in desired range from PMC object
    - loop over remaining stations in PMC object
      - assign picks (= fill pickInfo for station)
      - write station as XML file NO.LOC.STA.pickInfo.xml
      - write XML station file name to stationlist file

    run this task in parallel on many nodes
    -> comppickInfo-NNN.py (200 station chunks) on othello.ethz.ch, needs ~8 GB RAM
    -> dy-comppickInfo-NNN.py (50 station chunks) on dynamic.usc.edu, needs ~2 GB RAM
    """
    
    station_chunk = ( 10, 19 )
    
    catdir         = './data/catalog-bz2/'

    catfiles       =  ( catdir + '200007_A.deck.Z.bz2', )
    # catfiles = sort( glob.glob( catdir + '*.Z.bz2' ) )

    stationfile      = './data/jma_station_list_with_network_code.dat'
    stationlist_file = './data/pickInfo/pmc.jma.%s-%s.pickInfo.filelist.dat' % ( station_chunk[0], station_chunk[1] )

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_JMA( stationfile, 0, smooth = True, encoding = 'html' )
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

    # loop over catalog chunk pickles
    for catfile in catfiles:
        print " ----- importing catalog chunk from file: ", catfile

        qpc = QPCatalog( idstyle='numeric' )
        qpc.importJMADeck( catfile, compression='bz2', minimumDataset=True )
        print " catalog chunk has %s events" % qpc.size()

        pmc.preprocessCatalog( qpc, [ DateTime(2000, 1, 1, 0, 0, 0), DateTime(2008, 4, 1, 0, 0, 0) ] )
        print " after catalog preprocess: catalog chunk has %s events" %  qpc.size()

        print " assigning picks ... "
        pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

        del qpc

    print " after assign picks: inventory has ", len(pmc.stations), " stations"

    # over stations - write station with pickInfo separately as XML
    print " +++++ NOW WRITING STATION DATA AS XML +++++ "

    # open stationlist file
    flist = writeQPData( stationlist_file )
    
    for curr_sta_idx in reversed( range( len( pmc.stations ) ) ):

        sta_number = station_chunk[0] + curr_sta_idx
        curr_sta = pmc.stations[curr_sta_idx]

        sta_filename_xml = str(sta_number) + '.' + curr_sta.locationCode + '.' + curr_sta.stationCode + '.pickInfo.xml.bz2'
        sta_file_xml     = './data/pickInfo/' + sta_filename_xml
    
        print " processing station %s: %s %s" % ( sta_number, curr_sta.locationCode, curr_sta.stationCode )

        print "  writing to XML"
        fh = writeQPData( sta_file_xml, compression='bz2' )
        fh.write( '<?xml version="1.0" encoding="utf-8"?>' )
        fh.write( '<PMC>' )
        curr_sta.toXML( 'PMCStation', fh )
        fh.write( '</PMC>' )
        fh.close()

        # write station filename to list file
        flist.write( ''.join( ( sta_filename_xml, '\n' ) ) )

        # delete station just processed
        del pmc.stations[curr_sta_idx]

    flist.close()
    
def task_6():
    """
    import JMA catalog from 200001-200803 from JMA deck format
    do not use phase data
    export to 10-column ZMAP format, bz2-compressed
    """
    
    catdir         = './data/catalog-bz2/'
    catfile_zmap   = './data/jma.200001-200803.zmap.bz2'
    
    catfiles = sort( glob.glob( catdir + '*.Z.bz2' ) )
    
    qpc = QPCatalog( idstyle='numeric' )

    curr_event_ct = 0
    
    # loop over catalog chunk pickles
    for catfile in catfiles:
        print " ----- importing catalog chunk from file: ", catfile
        qpc.importJMADeck( catfile, compression='bz2', minimumDataset=True, nopicks=True )

        print " catalog chunk has %s events, %s overall" % ( qpc.size() - curr_event_ct, qpc.size() )
        curr_event_ct = qpc.size()

    # export catalog to ZMAP
    print " +++++ NOW WRITING CATALOG DATA TO ZMAP FORMAT, FILE %s +++++ " % catfile_zmap
    qpc.exportZMAP( catfile_zmap, compression='bz2' )
    
def task_7():
    """
    convert distros of all stations from xml to ASCII
    catalog from 2007-04 to 2008-03
    """

    useDistStyle  = 4

    basepath = '/home/fab/prog/pyprog/quakepy/pmc/jma/'
    
    distro_xml_path = basepath + 'data/distro-yearly/200704-200803/'
    distro_dat_path = basepath + 'data/distro-ascii/200704-200803/'
    
    distrofiles = sort( glob.glob( distro_xml_path + '*.xml.bz2' ) )
    
    staCtr = 0
    for curr_distrofile in distrofiles:

        staCtr += 1
        
        pmc = PMC_JMA()
        
        print " importing PMC from XML file: ", curr_distrofile
        pmc.readXML( getQPDataSource( curr_distrofile, compression='bz2' ) )
        print " read PMC inventory with ", len(pmc.stations), " stations"

        if len(pmc.stations) > 1:
            print " TROUBLE: i have too many stations!"

        # select station
        use_sta = pmc.stations[0]

        if use_sta.distribution.distro == None:
            print " no distro for station no. %s: %s %s" % ( staCtr, use_sta.locationCode, use_sta.stationCode )
            continue

        # get filename
        distro_data_file = distro_dat_path + str(staCtr) + '.' + use_sta.locationCode + '.' + use_sta.stationCode + '.200704-200803.distro.dat.bz2'

        # write distro as ASCII
        print " writing station no. %s: %s %s to data file %s" % ( staCtr, use_sta.locationCode, use_sta.stationCode, distro_data_file )
        
        fh = writeQPData( distro_data_file, compression='bz2' )

        for curr_mag_idx in xrange( len(use_sta.distribution.magnitudeGrid ) ):

            line_arr = []
            for curr_dist_idx in xrange( len(use_sta.distribution.distanceGrid) ):

                line_arr.append( "%7.5e" % use_sta.distribution.distro[curr_dist_idx, curr_mag_idx]  )
                        

            fh.write( ( '\t'.join( line_arr ) + '\n' ) )
        
        fh.close()

        del pmc

def task_8():
    """ 
    read PMC distros for all stations from PMC run
    plot distros (create EPS, convert to PNG, remove EPS)
    """

    metadata = PMCMetadata()
    
    # this gets the absolute base path  == path where this script resides
    # subsequent directory hierarchy will be below this directory
    scriptname = os.path.basename( os.path.abspath( sys.argv[0] ) )
    basepath   = os.path.dirname( os.path.abspath( sys.argv[0] ) ) # no trailing slash

    useDate = '2007-04-02'
    
    # directory where data from this run will be written to
    rundir_base = os.path.join( basepath, 'data/runs/old' )
    runPath     = os.path.join( rundir_base, 't002' )

    # directory where yearly probability distributions will be written to
    distroDir = os.path.join( runPath, useDate, 'distro/' )

    # directory for eps figures of distros
    distroEPSDir = os.path.join( runPath, useDate, 'distro-plot/' )

    # create eps path
    if not os.path.isdir( distroEPSDir ):
        print "   creating path %s" % distroEPSDir
        os.makedirs( distroEPSDir )

    # distStyle for prob distros
    useDistStyle = 5

    # color table
    # colortable_file = '../scsn/data/colortables/rastafari.cpt'
    metadata.probDistroColorMapFile = '../../data/colormaps/pov/rastafari.pov'
    
    # ---------------------------------------------------------------------------------

    # import PMC distro data: loop over distro files for stations
    distrofiles = sorted( glob.glob( os.path.join( distroDir, '*.distro.xml.bz2' ) ) )

    for sta_ctr, curr_sta_file in enumerate( distrofiles ):

        print " importing PMC from XML file: %s" % curr_sta_file

        curr_pmc = PMC_JMA()
        curr_pmc.readXML( getQPDataSource( curr_sta_file, compression='bz2' ) )

        print " --> imported %s station(s), %s of %s input files" % ( len( curr_pmc.stations ),
                                                                      sta_ctr+1,
                                                                      len( distrofiles ) )

        curr_pmc.setDistroStyle( useDistStyle )

        # get name for eps plot
        pmc_png_filename = '.'.join( ( curr_pmc.stations[0].locationCode, curr_pmc.stations[0].stationCode,
                                       'distro.png' ) )

        pmc_png_out = os.path.join( distroEPSDir, pmc_png_filename )
        
        curr_pmc.plotDistribution( curr_pmc.stations[0], metadata, pmc_png_out,
                                   colortable=metadata.probDistroColorMapFile )

        del curr_pmc


def main():
    
    #task_1()
    #task_2()
    #task_3()
    #task_4()
    #task_5()
    #task_6()
    #task_7()
    task_8()
    #task_9()
    #task_10()
    #task_11()
    #task_12()
    #task_13()
    #task_14()
        
main()