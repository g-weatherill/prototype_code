#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# quakepy/pmc/scsn/pmc_scsn.py
# $Id: pmc_scsn.py 78 2008-06-17 16:20:56Z fab $
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

__version__  = '$Id: pmc_scsn.py 78 2008-06-17 16:20:56Z fab $'
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

from mx.DateTime     import Date
from mx.DateTime     import DateTime
from mx.DateTime.ISO import ParseDateTimeUTC

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

#from qpplot           import *
#from qpseismicityplot import *

def task_1():
    """
    create PMC with catalog for all years (in chunks)
    import full station list, alias list
    preprocess inventory & catalog, assign picks, merge alias stations, postprocess
    pickle inventory
    ---
    OBSOLETE -- THIS WAS WRITTEN FOR AN OLDER VERSION OF QUAKEPY
    DOES NOT RUN ANY LONGER
    """
    exit()
    catfiles       =  ( '../data/phase2001.dat',
                        '../data/phase2002.dat',
                        '../data/phase2003.dat',
                        '../data/phase2004.dat',
                        '../data/phase2005.dat',
                        '../data/phase2006.dat',
                        '../data/phase200701-06.dat' )
    
    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'

    raw_inv_pickle = 'inventory.raw.pickle'
    inv_pickle     = 'inventory.nopickinfo.pickle'

    m0 = memory() 
    
    print " importing inventory file: ", stationfile
    pmci   = quakepy.pmc.pmc.PMCInventory( stationfile )
    print " before preprocess: inventory has ", len( pmci.stations ), " stations"

    m1 = memory( m0 )
    print " memory consumed for inventory: ", m1 / 1000000

    pmc = quakepy.pmc.pmc.PMC_SCSN( pmci )
    pmc.importAliases( aliasfile )
    pmc.preprocessInventory( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.inventory.stations), " stations"

    m2 = memory( m1 )
    print " memory consumed for PMC inventory / import aliases: ", m2 / 1000000
    
    # loop over catalog tranches
    for catfile in catfiles:
        print " ----- importing uncompressed STP phase catalog file: ", catfile

        qpc = quakepy.quakepy.QPCatalog()
        qpc.importSTPPhase( catfile )
        print " catalog has ", len( qpc.eventParameters.event ), " events"

        m_curr = memory( m0 )
        print " total memory consumed after catalog chunk: ", m_curr / 1000000
    
        pmc.setCatalog( qpc )
        pmc.preprocessCatalog( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
        print " after catalog preprocess: catalog has ", len( pmc.catalog.eventParameters.event ), " events"

        print " assigning picks ... "
        pmc.assignPicks( reduce_catalog = False, verbose_level = 2 )

        m_curr = memory( m0 )
        print " total memory consumed after assign picks: ", m_curr / 1000000

    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.inventory.stations), " stations"
    
    m_end = memory( m0 )
    print " total memory consumed after assign picks: ", m_end / 1000000

    # -------------------------------------------------------------------------------
    
    print " pickling raw inventory with pickInfo / without distro ..."
    pmc.saveInventory( raw_inv_pickle )
    print "ready"

    print " call fillup() for all stations and delete pickInfo"
    for curr_sta in pmc.inventory.stations:
        print " processing station %s" % ( curr_sta.stationCode )
        curr_sta.distribution.fillup( pmc._calcMagnitudeFromDistance )
        del curr_sta.distribution.pickInfo
        
    print " pickling inventory without pickInfo ..."
    pmc.saveInventory( inv_pickle )
    print "ready"

def task_2():
    """
    for all catalog chunks we have: read in as STP, write out as pickle
    """
    exit()
    catnames       =  ( '../data/phase2001',
                        '../data/phase2002',
                        '../data/phase2003',
                        '../data/phase2004',
                        '../data/phase2005',
                        '../data/phase2006',
                        '../data/phase200701-06',
                        '../data/phase2001-200101',
                        '../data/phase.mini' )

    for catname in catnames:

        catfile    = catname + '.dat'
        picklefile = catname + '.pickle'
        
        print " ----- importing uncompressed STP phase catalog file: ", catfile

        qpc = QPCatalog()
        qpc.importSTPPhase( getQPDataSource( catfile ) )
        print " catalog has ", len( qpc.eventParameters.event ), " events"

        print " ----- pickling catalog chunk: ", picklefile
        qpc.save( picklefile )

        del qpc

def task_3():
    """
    for full Southern California catalog: read in as STP, write out as pickle
    NOTE: WILL TAKE LONG TIME & A LOT OF MEMORY
    CRASHES with MemoryError ON 4GB machine (othello)
    FIND OUT: does it crash while reading STP or when pickling?
    """
    catnames       =  ( '../data/phase2001-200706', )

    for catname in catnames:

        catfile    = catname + '.dat'
        picklefile = catname + '.pickle'
        
        print " ----- importing uncompressed STP phase catalog file: ", catfile

        qpc = QPCatalog()
        qpc.importSTPPhase( getQPDataSource( catfile ) )
        print " catalog has ", len( qpc.eventParameters.event ), " events"

        print " ----- pickling catalog chunk: ", picklefile
        qpc.save( picklefile )

        del qpc

def task_4():
    """
    compute pickInfo for all stations
    create PMC_SCSN with full station list, distStyle 0 so that we have no <distro> element in XML
    loop over catalog chunks, read in catalog chunk from STP file

    loop over stations:
     - assign picks (= fill pickInfo for station)
     - write station as XML snippet NO.NC.STA.pickInfo.xml
     - write station as pickle      NO.NC.STA.pickInfo.pickle
     - write XML station file name to stationlist file

    write PMC object with pickInfo as XML
    write PMC object with pickInfo as pickle
    """
    
    #catfiles       =  ( '../data/phase2001.pickle',
                        #'../data/phase2002.pickle',
                        #'../data/phase2003.pickle',
                        #'../data/phase2004.pickle',
                        #'../data/phase2005.pickle',
                        #'../data/phase2006.pickle',
                        #'../data/phase200701-06.pickle' )

    #catfiles       =  ( '../data/phase2001-200101.pickle', )

    
    catfiles       =  ( '../data/phase2001.dat',
                        '../data/phase2002.dat',
                        '../data/phase2003.dat',
                        '../data/phase2004.dat',
                        '../data/phase2005.dat',
                        '../data/phase2006.dat',
                        '../data/phase200701-06.dat' )

    stationfile   = '../data/stationlist.dat'
    #stationfile   = '../data/stationlist.test.dat' # POB, PSP
    
    aliasfile     = '../data/sc_alias.txt'

    pmc_xml_file    = '../data/pickInfo/pmc.sc.pickInfo.xml'
    pmc_pickle_file = '../data/pickInfo/pmc.sc.pickInfo.pickle'

    stationlist_file = '../data/pickInfo/pmc.sc.pickInfo.filelist.dat'

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, 0, smooth = True, encoding = 'html' )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"

    pmc.importAliases( aliasfile )
    pmc.preprocessInventory( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    # loop over catalog chunk pickles
    for catfile in catfiles:
        print " ----- importing catalog chunk from file: ", catfile

        #qpc = unpickleObj( catfile )
        qpc = QPCatalog()
        qpc.importSTPPhase( getQPDataSource( catfile ) )
        print " catalog has ", len( qpc.eventParameters.event ), " events"

        pmc.preprocessCatalog( qpc, [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
        print " after catalog preprocess: catalog has ", len( qpc.eventParameters.event ), " events"

        print " assigning picks ... "
        pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

        del qpc


    print " merging aliased stations ... "
    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.stations), " stations"

    # over stations - write station with pickInfo separately as XML and pickle
    print " +++++ NOW WRITING STATION DATA AS XML / PICKLE +++++ "

    # open stationlist file
    flist = writeQPData( stationlist_file )
    
    for curr_sta_idx, curr_sta in enumerate( pmc.stations ):

        sta_filename_xml    = str(curr_sta_idx) + '.' + curr_sta.networkCode + '.' + curr_sta.stationCode + '.pickInfo.xml'
        sta_file_xml        = '../data/pickInfo/' + sta_filename_xml
        sta_filename_pickle = str(curr_sta_idx) + '.' + curr_sta.networkCode + '.' + curr_sta.stationCode + '.pickInfo.pickle'
        sta_file_pickle     = '../data/pickInfo/' + sta_filename_pickle
    
        print " processing station %s: %s %s" % ( curr_sta_idx, curr_sta.networkCode, curr_sta.stationCode )

        print "  writing to XML"
        fh = writeQPData( sta_file_xml )
        fh.write( '<?xml version="1.0" encoding="utf-8"?>' )
        fh.write( '<PMC>' )
        curr_sta.toXML( 'PMCStation', fh )
        fh.write( '</PMC>' )
        fh.close()

        # write station filename to list file
        flist.write( ''.join( ( sta_filename_xml, '\n' ) ) )

        print "  writing to pickle"
        pickleObj( curr_sta, sta_file_pickle )

    flist.close()
        
    # write PMC as XML
    print "  writing PMC to XML"
    pmc.writeXML( writeQPData( pmc_xml_file ) )

    # write PMC as pickle
    print "  writing PMC to pickle"
    #pmc.save( pmc_pickle_file )
    pickleObj( pmc, pmc_pickle_file )

def task_5():
    """
    test for computing distros on cluster
    1) load file list of pickInfo XML files
    2) load PMC for one station with pickInfo from XML file / get filename from file list (idx.NC.STA.pickInfo.xml)
    3) reset distStyle to 3
    4) call fillup on station
    5) write PMC (with only one station) as XML: idx.NC.STA.pickInfo.distro.xml
    6) delete pickInfo
    7) write PMC (with only one station) as XML: idx.NC.STA.distro.xml
    """

    useStationIdx = 280 # CI / POB
    useDistStyle  = 3

    pickinfo_filelist = '../data/pickInfo/pmc.sc.pickInfo.filelist.dat'
    pickinfo_path     = '../data/pickInfo/'

    pickinfo_distro_xml_path = '../data/distro/'
    distro_xml_path          = '../data/distro/'

    # read pickinfo filelist
    pickinfo_file_arr = getQPDataSource( pickinfo_filelist ).read().split()
    print " i could find %s pickinfo files" % ( len(pickinfo_file_arr) )

    #pmcfile_xml = '../data/pickInfo/' + str(useStationIdx) + '.' + use_sta.networkCode + '.' + use_sta.stationCode + '.pickInfo.xml'
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

def task_6():
    """
    get time stamps of all onTimes of all stations
    write to file
    for these times we later compute the probability grids
    """

    stationfile   = '../data/stationlist.dat'
    #stationfile   = '../data/stationlist.test.dat' # POB, PSP
    
    timestamp_file = '../data/grid/pmc.sc.ontimes.dat'

    startDate = Date( 2001, 1, 1 )
    endDate   = Date( 2007, 7, 1 ) 

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, 0, smooth = True, encoding = 'html' )
    print " PMC inventory has ", len(pmc.stations), " stations"

    pmc.preprocessInventory( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    time_arr = [ mxDateTime2ISO( [ startDate ] ) ]
    
    # over stations - write timestamp file
    print " +++++ NOW WRITING TIMESTAMP FILE +++++ "

    # open timestamp file
    flist = writeQPData( timestamp_file )
    
    for curr_sta in pmc.stations:

        for curr_cha in curr_sta.channels:

            for curr_time in curr_cha.onTime:

                for idx in xrange( 2 ):

                    # add string representation of time, since addUnique
                    # does not work with DateTime objects
                    if (     curr_time[idx] >= startDate
                         and curr_time[idx] <= endDate ):
                        
                        datetime_str = mxDateTime2ISO( [ curr_time[idx] ] )
                        #print " checking %s -> %s" % ( curr_time[idx], datetime_str )
                        addUnique( time_arr, [ datetime_str ] )

    print " found %s unique time stamps" % ( len( time_arr ) )

    flist.write( '\n'.join( reversed( sort( time_arr ) ) ) )
    flist.write( '\n' )

    flist.close()

def task_7():
    """
    combine all one-station distro XML files in one big PMC with all stations
    """
    basepath = '/home/geovault-00/fabian/prog/pyprog/quakepy/pmc/'

    pickinfo_path     = basepath + 'data/pickInfo/'
    distro_xml_path   = basepath + 'data/distro/'

    pickinfo_filelist = pickinfo_path + 'pmc.sc.pickInfo.filelist.dat'
    pmc_xml_out       = distro_xml_path + 'scsn.distro.xml'

    # get XML distro files from ASCII pickinfo filelist
    print " opening %s ..." % ( pickinfo_filelist )
    pickinfo_file_arr = getQPDataSource( pickinfo_filelist ).read().split()
    print " i could find %s pickinfo files" % ( len(pickinfo_file_arr) )

    print " ----- get PMC from XML -----  "

    # create master PMC object with first XML file in list (Index 0)
    distro_xml_filename = distro_xml_path + pickinfo_file_arr[0][:-12] + 'distro.xml'

    pmc = PMC_SCSN()
    print " importing PMC from XML file: ", distro_xml_filename
    pmc.readXML( getQPDataSource( distro_xml_filename ) )
    print " loaded PMC from XML: %s stations" % ( len( pmc.stations ) )

    # over other XML files, start from 1
    for idx in xrange( 1, len(pickinfo_file_arr) ):
    #for idx in xrange( 2, 3 ):
        
        distro_xml_filename = distro_xml_path + pickinfo_file_arr[idx][:-12] + 'distro.xml'
        
        pmc_add = PMC_SCSN()
        print " adding PMC from XML file: ", distro_xml_filename
        pmc_add.readXML( getQPDataSource( distro_xml_filename ) )

        pmc.merge( pmc_add )
        print " merged station into PMC: now %s stations" % ( len( pmc.stations ) )

    # write overall PMC as XML
    print " writing XML file: ", pmc_xml_out
    fh = writeQPData( pmc_xml_out )
    pmc.writeXML( fh )
    fh.close()

def task_8():
    """
    compute combination list
    NOTE: this is impossible on 2GB memory, execution takes longer and longer the higher num and sta_ctr
    """
    max_bignumber   = 500
    max_smallnumber = 4

    list_file = 'combinations-' + str(max_smallnumber) + '-' + str(max_bignumber) + '.pickle'
    
    combo = []
    
    for num in xrange( 1, max_smallnumber+1, 1 ):

        combo.append( [] )
        
        for sta_ctr in xrange( max_smallnumber, max_bignumber+1, 1 ):
            combo[num-1].append( [] )
            combo[num-1][sta_ctr-max_smallnumber] =  list( xuniqueCombinations( range(sta_ctr), num ) )

            #print " computed %s, %s, %s" % ( num, sta_ctr, list( xuniqueCombinations( range(sta_ctr), num ) ) )
            print " computed %s, %s " % ( num, sta_ctr )
            
    pickleObj( combo, list_file )

def task_9():
    """
    compute max number of stations for all timestamps
    """
    basepath = '/home/fab/prog/pyprog/quakepy/pmc/'

    distro_xml_path   = basepath + 'data/distro/'
    grid_xml_path     = basepath + 'data/grid/'

    pmc_xml_in       = distro_xml_path + 'scsn.distro.xml'
    date_filelist    = grid_xml_path + 'pmc.sc.ontimes.dat'

    # get timestamps
    timestamps_arr = getQPDataSource( date_filelist ).read().split()

    sta_no_arr = []
    
    for timestamp in timestamps_arr:
        
        t = ParseDateTimeUTC(timestamp)

        pmc = PMC_SCSN()

        print " importing PMC from XML file: ", pmc_xml_in
        pmc.readXML( getQPDataSource( pmc_xml_in ) )
        print " loaded PMC from XML: %s stations" % ( len( pmc.stations ) )

        pmc._selectStationsInOperation( t )
        sta_no = len( pmc.stations )
        
        sta_no_arr.append( sta_no )
        print timestamp, sta_no

    print "maximum: ", max( sta_no_arr )
        
def task_10():
    """
    compute combination list class
    use NumPy version
    class has changed from revision 57, old pickles won't work any more
    
    Southern California: M=3, N=375
    Northern California: M=3, N=630
    """
    
    max_N = 630   # 200, 375, 630
    max_M = 3

    pickle_file_py   = 'combinations-' + str(max_M) + '-' + str(max_N) + '.py.pickle'
    pickle_file_numpy = 'combinations-' + str(max_M) + '-' + str(max_N) + '.numpy.pickle'

    # use native python
    print " Native Python: compute combinations for M=%s, N=%s" % ( max_M, max_N )
    c = Combinations( max_M, max_N )

    print " pickling combinations object"
    pickleObj( c, pickle_file_py )

    del c
    
    # use numpy - MEMORY ERROR ON othello (4 / 3.25 GB) FOR 3 / 630 when creating numpy.ones for m=3 !!
    print " NumPy: compute combinations for M=%s, N=%s" % ( max_M, max_N )
    c = CombinationsNumPy( max_M, max_N )

    print " pickling combinations object"
    pickleObj( c, pickle_file_numpy )
    
    del c
    

def task_11():
    """
    convert combinations pickle from native Python to numpy
    """
    
    max_N = 630   # 200, 375, 630
    max_M = 3

    pickle_file_py   = '../data/combinations-' + str(max_M) + '-' + str(max_N) + '.pickle'
    pickle_file_numpy = 'combinations-' + str(max_M) + '-' + str(max_N) + '.numpy.pickle'

    # use native python
    print " read native Python pickle for M=%s, N=%s" % ( max_M, max_N )
    c_nat = unpickleObj( pickle_file_py )

    # use numpy
    print " create NumPy: compute combinations for M=%s, N=%s" % ( max_M, max_N )
    c_num = CombinationsNumPy( max_M, max_N, nocompute = True )

    # transfer data from native python to numpy
    
    # copy index list
    c_num.index = copy.deepcopy( c_nat.index )

    # copy combinations list
    # over m: 2, ..., M
    for m in xrange( 2, c_nat.M+1, 1):

        dim = len( c_nat.combinations[m] )

        c_num.combinations[m] = numpy.ones( ( dim, m ), dtype=int ) * numpy.nan

        for ctr in xrange( dim ):

            print " transfer m %s, line %s: %s " % ( m, dim, c_nat.combinations[m][ctr] )
            c_num.combinations[m][ctr] = c_nat.combinations[m][ctr]


    print " pickling combinations object"
    pickleObj( c_num, pickle_file_numpy )

    #print " COMPARE"
    #print " native: index %s, combi %s" % ( c_nat.index, c_nat.combinations )
    #print " numpy: index %s, combi %s" % ( c_num.index, c_num.combinations )
    
def task_12():
    """
    compare memory consumption of native Python combinations to numpy combinations
    """
    
    max_N = 630   # 200, 375, 630
    max_M = 3

    pickle_file_py   = '../data/combinations-' + str(max_M) + '-' + str(max_N) + '.pickle'
    pickle_file_numpy = '../data/combinations-' + str(max_M) + '-' + str(max_N) + '.numpy.pickle'

    m0 = memory()
    
    # use native python
    print " read native Python pickle for M=%s, N=%s" % ( max_M, max_N )
    c_nat = unpickleObj( pickle_file_py )

    m1 = memory( m0 )
    print " memory used for native Python: ", m1 / 1000000
    del c_nat

    # use numpy python
    print " read NumPy pickle for M=%s, N=%s" % ( max_M, max_N )
    c_num = unpickleObj( pickle_file_numpy )

    m2 = memory( m1 )
    print " memory used for NumPy: ", m2 / 1000000

def task_13():
    """
    create list of active stations for all timestamps for SCSN network

    format:
    network_code | station_code | longitude | latitude | elevation

    filenames:
    scsn.YYYY-MM-DDTHH:MM:SS.station.dat
    """

    basepath = '/home/fab/prog/pyprog/quakepy/pmc/'

    distro_xml_path   = basepath + 'data/distro/'
    grid_xml_path     = basepath + 'data/grid/'
    station_dat_path  = basepath + 'data/station/'
    
    pmc_xml_in       = distro_xml_path + 'scsn.distro.xml'
    date_filelist    = grid_xml_path + 'pmc.sc.ontimes.dat'

    # get timestamps
    timestamps_arr = getQPDataSource( date_filelist ).read().split()
    
    for timestamp in timestamps_arr:
        
        t = ParseDateTimeUTC(timestamp)

        pmc = PMC_SCSN()

        print " importing PMC from XML file: ", pmc_xml_in
        pmc.readXML( getQPDataSource( pmc_xml_in ) )
        print " loaded PMC from XML: %s stations" % ( len( pmc.stations ) )

        # no need to preprocessInventory, since the XML files we read are already preprocessed

        pmc._selectStationsInOperation( t )

        # write station file
        sta_filename = station_dat_path + 'scsn.' + mxDateTime2ISO( [ timestamp ] ) + '.station.dat'

        print " writing file %s with %s stations" % ( sta_filename, len( pmc.stations ) )
        
        sta_fh = open( sta_filename, 'w' )

        # over stations
        for curr_sta in pmc.stations:

            sta_out = ( curr_sta.networkCode,
                        curr_sta.stationCode,
                        '%9.4f' % curr_sta.longitude,
                        '%9.4f' % curr_sta.latitude,
                        '%6.1f' % curr_sta.elevation
                      )

            sta_fh.writelines( ( '\t'.join( sta_out ), os.linesep ) )

        sta_fh.close()

def task_14():
    """
    compute completeness on given grid for specific point in time
    write grid with annotation as XML
    """
    
    useDate       = ( 2007, 06, 01 )
    useDistStyle  = 3
    
    maggrid       = frange( 0.0, 4.0, 1.0 )
    mp_prob_list  = [ 0.9, 0.999 ]
    
    # set geographical region for completeness evaluation
    geobox = { 'lonmin': -122.0, 'lonmax': -113.5, 
               'latmin': 31.5, 'latmax': 38.0,
               'londelta': 3.0, 'latdelta': 3.0,
               'depthmin': 7.5, 'depthmax': 7.5 }
               
    basepath = '/home/fab/prog/pyprog/quakepy/pmc/scsn/'
    
    distro_xml_path   = basepath + 'data/distro/'
    grid_xml_path     = basepath
    combi_path        = basepath + 'data/'
    
    pmc_xml_in        = distro_xml_path + 'scsn.distro.xml'
    combi_pickle      = combi_path + 'combinations-3-375.numpy.pickle.r56'
    
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
    mapdate_str = mxDateTime2ISO( mapdate )
    
    # set grid annotation
    g.annotation.setProperty( date = utc(),
                              starttime = mapdate,
                              endtime = mapdate,
                              latmin = geobox['latmin'],
                              latmax = geobox['latmax'],
                              lonmin = geobox['lonmin'],
                              lonmax = geobox['lonmax'] )
                              
    print " compute probability map for Date: %s" % ( mapdate_str )
    pmc.getProbabilityMap( g, maggrid, mapdate, mp_prob_list, verbose = False )
    
    # add station data to grid - delete not required fields
    g.stations = []
    for curr_sta in pmc.stations:
        del curr_sta.channels
        del curr_sta.distribution
        g.stations.append( curr_sta )
    
    grid_file_out = grid_xml_path + 'scsn.test.' + mapdate_str + '.grid.xml'
    print " write computed grid file ", grid_file_out
    g.writeXML( writeQPData( grid_file_out ) )


    
## ---------------------------------------------------------------------------

    
def task_99():
    """
    create PMC for Southern California (PMC_SCSN) with catalog for all years (in chunks)
    import full station list, alias list
    preprocess inventory & catalog, assign picks, merge alias stations
    export PMC object without distro / with pickInfo as XML
    compute distro with style 3
    delete pickInfo
    export PMC object with distro / without pickInfo as XML
    """
    catfiles       =  ( '../data/phase2001.dat',
                        '../data/phase2002.dat',
                        '../data/phase2003.dat',
                        '../data/phase2004.dat',
                        '../data/phase2005.dat',
                        '../data/phase2006.dat',
                        '../data/phase200701-06.dat' )
    
    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'

    pmc_raw_xml_file    = 'pmc.scsn.raw.xml.bz2'
    pmc_nopick_xml_file = 'pmc.scsn.nopickinfo.xml.bz2'

    m0 = memory() 

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, 3 )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"

    m1 = memory( m0 )
    print " memory consumed for PMC inventory: ", m1 / 1000000
    
    pmc.importAliases( aliasfile )
    pmc.preprocessInventory( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    m2 = memory( m1 )
    print " memory consumed for PMC inventory / import aliases: ", m2 / 1000000

    # loop over catalog tranches
    for catfile in catfiles:
        print " ----- importing uncompressed STP phase catalog file: ", catfile

        qpc = QPCatalog()
        qpc.importSTPPhase( catfile )
        print " catalog has ", len( qpc.eventParameters.event ), " events"

        m_curr = memory( m0 )
        print " total memory consumed after catalog chunk: ", m_curr / 1000000
    
        pmc.preprocessCatalog( qpc, [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
        print " after catalog preprocess: catalog has ", len( qpc.eventParameters.event ), " events"

        print " assigning picks ... "
        pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

        m_curr = memory( m0 )
        print " total memory consumed after assign picks: ", m_curr / 1000000


    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.stations), " stations"
    
    m_end = memory( m0 )
    print " total memory consumed after assign picks: ", m_end / 1000000

    # -------------------------------------------------------------------------------

    # write PMC with pickInfo to XML
    # CRASHES! maybe too big?
    #print " write PMC with pickInfo as XML file", pmc_raw_xml_file
    #pmc.writeXML( pmc_raw_xml_file, compression = 'bz2'  )
    #print "ready"

    # fillup distro, delete pickInfo
    print " call fillup() for all stations and delete pickInfo"
    for curr_sta_no, curr_sta in enumerate( pmc.stations ):
        print " processing station no. %s: %s" % ( curr_sta_no+1, curr_sta.stationCode )
        curr_sta.distribution.fillup( pmc._calcMagnitudeFromDistance )
        del curr_sta.distribution.pickInfo
        
    print " write PMC with distro / without pickInfo as XML file", pmc_nopick_xml_file
    pmc.writeXML( pmc_nopick_xml_file, compression = 'bz2'  )
    print "ready"

    
def main():
    
    #task_1()
    #task_2()
    #task_3()
    #task_4()
    #task_5()
    #task_6()
    #task_7()
    #task_8()
    #task_9()
    #task_10()
    #task_11()
    #task_12()
    #task_13()
    task_14()
        
main()