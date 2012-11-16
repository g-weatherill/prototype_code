#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# quakepy/pmc/scsn/examples/pmc_scsn_tests.py
# $Id: pmc_scsn_tests.py 65 2008-05-27 16:08:53Z fab $
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

__version__  = '$Id: pmc_scsn_tests.py 65 2008-05-27 16:08:53Z fab $'
__revision__ = '$Revision: 65 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import random
import datetime

import pickle
import cPickle

import matplotlib
from matplotlib import rcParams
from matplotlib.dates import MO
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, WeekdayLocator, DateFormatter
matplotlib.use('PS')

from pylab import *

import Polygon
from   mx.DateTime     import DateTime
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
from qpplot           import *
from qpseismicityplot import *
 

def test_1():
    """
    read catalog from "phase file", pickle catalog, and plot events (do not use phase data)
    """
    #infile   = '../pmc/data/phase2001-200706.dat'
    #infile      = '../data/phase2001.dat'
    infile           = '../data/phase2001-200101.dat'
    picklefile       = 'phase2001-200101.pickle'
    picklefile_check = 'phase2001-200101.pickle.check'
    outfile_mtl      = 'phase2001_mtl.eps'
    outfile_gmt      = 'phase2001_gmt.eps'
    
    qpc = QPCatalog()
    qpc.importSTPPhase( infile )
    print " read uncompressed STP phase catalog with ", len( qpc.eventParameters.event ), " events" 

    # qpc.cut( minmag=5.0 )
    plot = qpseismicityplot.QPSeismicityPlot().plot_matplotlib( qpc, outfile_mtl, mpl_coastlines=True )
    
    #plot = qpseismicityplot.QPSeismicityPlot().plot_gmt( qpc, 
                                                         #outfile_gmt, 
                                                         #gmt_proj='merc', 
                                                         #gmt_event_symbol_type='c', 
                                                         #gmt_event_symbol_color=(0,0,255), 
                                                         #gmt_basemaptype='fancy' )
                                                         
    # pickle catalog
    print " pickling catalog ..."
    pickleObj( qpc, picklefile )
    print "ready"

    # check pickle - load catalog from pickle
    print " loading catalog pickle: ", picklefile
    qpc_check = unpickleObj( picklefile )
    print " checked catalog has ", len( qpc_check.eventParameters.event ), " events"

    # write new pickle from checked catalog
    print " pickling checked catalog ..."
    pickleObj( qpc_check, picklefile_check )
    print "ready"
    
    try:
        assert qpc_check == qpc, "catalog instances are not equal"
    except AssertionError, args:
        print '%s: %s' % ( args.__class__.__name__, args )
        exit()
        
    print " SUCCESS: catalog instances are equal "

def test_2():
    """
    create PMC from full station list, no catalog, look for specific station
    """
    stationfile = '../data/stationlist.dat'

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile )
    print " PMC_SCSN: PMC inventory has ", len( pmc.stations ), " stations"
    
    for curr_sta in pmc.stations:
        for curr_cha in curr_sta.channels:
            if curr_cha.waveformID.stationCode == 'POB':
                print curr_sta.__dict__
                print curr_cha.__dict__
                print curr_cha.waveformID.__dict__
    
def test_3():
    """
    create PMC with 2001 catalog and POB station, assignPicks, pickle inventory, plot raw distribution
    
    """
    catfile       = '../data/phase2001.dat'
    stationfile   = '../data/station_pollybutte.dat'
    inventoryfile = 'pmc-2001-POB.pickle'
    imgfile       = 'rawdistribution.eps'
    
    qpc = QPCatalog()
    qpc.importSTPPhase( catfile )
    print " read uncompressed STP phase catalog with ", len( qpc.eventParameters.event ), " events" 
    
    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile )
    print " PMC_SCSN: PMC inventory has ", len( pmc.stations ), " stations"
    
    pmc.assignPicks( qpc )
    print " wir haben ", len(pmc.stations), " stationen mit daten"
    print " wir haben ", len(pmc.stations[0].distribution.pickInfo), " picks an station 0"
    
    print " pickling PMC ..."
    pickleObj( pmc, inventoryfile )
    print "ready"
    
    print " plotting raw distribution ..."
    pmc.plotRawDistribution( pmc.stations[0], imgfile )
    
def test_4():
    """
    load PMC (POB, with assigned picks) from pickle, plot raw distribution
    """
    pmcpickle     = "pmc-2001-POB.pickle"
    imgfile       = "rawdistribution_frompickle.eps"

    print " loading PMC pickle: ", pmcpickle

    pmc = unpickleObj( pmcpickle )

    print " wir haben ", len(pmc.stations), " stationen mit daten"
    print " wir haben ", len(pmc.stations[0].distribution.pickInfo), " picks an station 0"
    
    print " plotting raw distribution ..."
    pmc.plotRawDistribution( pmc.stations[0], imgfile )
    
def test_5():
    """
    create PMC with 2001 catalog, full stations and preprocess inventory & catalog
    """
    catfile      = '../data/phase2001.dat'
    stationfile  = "../data/stationlist.dat"
    
    qpc = QPCatalog()
    
    print " importing uncompressed STP phase catalog file: ", catfile
    qpc.importSTPPhase( catfile )
    print " catalog has ", len( qpc.eventParameters.event ), " events" 

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"
    
    pmc.preprocess( qpc, [ Date(2001, 1, 1), Date(2007, 7, 1) ] )
    print " after preprocess: inventory has ", len(pmc.stations), " stations"
    
def test_6():
    """
    create PMC with 2001 catalog, full stations, alias list
    preprocess inventory & catalog, assign picks, merge alias stations, postprocess
    pickle inventory, read pickled inventory in new object, compare objects
    """
    #catfile       = '../data/phase2001-200706.dat'
    catfile       = '../data/phase2001.dat'
    #catfile       = '../data/phase2001-200101.dat'
    
    # stationfile   = '../data/stationlist.dat'
    stationfile   = '../data/station_pollybutte.dat'
    aliasfile     = '../data/sc_alias.txt'
    
    inventoryfile       = 'inventory2001.pickle'
    check_inventoryfile = 'check_inventory2001.pickle'
    
    qpc = QPCatalog()
    
    print " importing uncompressed STP phase catalog file: ", catfile
    qpc.importSTPPhase( catfile )
    print " catalog has ", len( qpc.eventParameters.event ), " events" 
    
    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile )
    print " before preprocess: inventory has ", len(pmc.stations), " stations"
    
    pmc.importAliases( aliasfile )
    pmc.preprocess( qpc, [ Date(2001, 1, 1), Date(2007, 7, 1) ] )
    print " after preprocess: inventory has ", len(pmc.stations), " stations"
    
    print " assigning picks ... "
    pmc.assignPicks( qpc )
    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.stations), " stations"
    
    print " pickling inventory ..."
    pickleObj( pmc, inventoryfile )
    print "ready"
    
    # check pickle
    print " loading PMC pickle: ", inventoryfile

    pmc_check = unpickleObj( inventoryfile )
    print " checked PMC inventory has ", len(pmc_check.stations), " stations"
    
    try:
        assert pmc_check == pmc, "PMC instances are not equal"
    except AssertionError, args:
        print '%s: %s' % ( args.__class__.__name__, args )
        exit()
        
    print " SUCCESS: inventories are equal "
    
def test_7():
    """
    create PMC with full catalog, full stations, alias list
    preprocess inventory & catalog, assign picks, merge alias stations, postprocess
    pickle inventory
    """
    catfile       = '../data/phase2001-200706.dat'
    catfile_zmap  = '../data/phase2001-200706.processed.zmap.dat'
    
    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'
    
    inventoryfile = 'inventory2001-200706.pickle'
    
    qpc = QPCatalog()
    
    print " importing uncompressed STP phase catalog file: ", catfile
    qpc.importSTPPhase( catfile )
    print " catalog has ", len( qpc.eventParameters.event ), " events" 
    
    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile )
    print " before preprocess: inventory has ", len(pmc.stations), " stations"
    
    pmc.importAliases( aliasfile )
    pmc.preprocess( qpc, [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after preprocess: inventory has ", len(pmc.stations), " stations"
    
    print " writing processed catalog as ZMAP, no of events = ", len(pmc.catalog.eventParameters.event)
    pmc.catalog.exportZMAP( catfile_zmap )
    print "ready"
    
    print " assigning picks ... "
    pmc.assignPicks( qpc )
    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.stations), " stations"
    
    print " pickling inventory ..."
    pickleObj( pmc, inventoryfile )
    print "ready"
    
def test_8():
    """
    create PMC with full catalog, full stations, alias list
    preprocess inventory & catalog
    pickle catalog
    CRASHES ON 4GB MACHINE, MemoryError
    """

    print " this test needs too much memory"
    exit()
    
    catfile       = '../data/phase2001-200706.dat'
    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'
    
    pickled_catfile = '../data/catalog2001-200706.pickle'
    
    qpc = QPCatalog()
    
    print " importing uncompressed STP phase catalog file: ", catfile
    qpc.importSTPPhase( catfile )
    print " catalog has ", len( qpc.eventParameters.event ), " events"

    print " pickling catalog ... "
    pickleObj( qpc, pickled_catfile )
    print "ready"
    
    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile )
    print " before preprocess: inventory has ", len(pmc.stations), " stations"
    
    pmc.importAliases( aliasfile )
    pmc.preprocess( qpc, [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after preprocess: inventory has ", len(pmc.stations), " stations"
    print " after preprocess: catalog has ", len(pmc.catalog.eventParameters.event), " events"
    

def test_9():
    """
    read in full SCSN catalog WITHOUT phases
    write catalog as XML
    pickle catalog
    """
    catfile         = '../data/phase2001-200706.dat'
    pickled_catfile = '../data/catalog2001-200706.nophases.pickle'
    xml_catfile     = '../data/catalog2001-200706.nophases.xml'
    
    qpc = QPCatalog()
    
    print " importing uncompressed STP phase catalog file: ", catfile
    qpc.importSTPPhase( getQPDataSource( catfile ), nopicks = True  )
    print " catalog has ", len( qpc.eventParameters.event ), " events" 

    print " exporting catalog as XML ... "
    qpc.writeXML( writeQPData( xml_catfile ) )
    print "ready"
    
    print " pickling catalog ... "
    pickleObj( qpc, pickled_catfile )
    print "ready"


def test_10():
    """
    create PMC with full catalog, full stations, alias list
    preprocess inventory & catalog, assign picks, merge alias stations, postprocess
    """
    catfile       = '../data/phase2001-200706.dat'
    #catfile       = '../data/phase2001.dat'
    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'
    
    m0 = memory() 
    
    qpc = QPCatalog()
    
    print " importing uncompressed STP phase catalog file: ", catfile
    qpc.importSTPPhase( catfile )
    print " catalog has ", len( qpc.eventParameters.event ), " events" 
    
    m1 = memory( m0 )
    print "memory consumed for catalog: ", m1 / 1000000
    
    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile )
    print " before preprocess: inventory has ", len(pmc.stations), " stations"

    m2 = memory( m1 )
    print "memory consumed for PMC inventory: ", m2 / 1000000
    
    pmc.importAliases( aliasfile )
    pmc.preprocess( qpc, [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after preprocess: inventory has ", len(pmc.stations), " stations"
    
    m3 = memory( m2 )
    print "memory consumed for PMC / import aliases: ", m3 / 1000000
    
    print " assigning picks ... " 
    pmc.assignPicks( qpc, reduce_catalog = True )
    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.stations), " stations"
    
    m4 = memory( m3 )
    m4_t = memory( m0 )
    print "total memory consumed after assign picks: ", m4_t / 1000000
    

def test_11():
    """
    create PMC with full catalog, in tranches for years
    import full stations, alias list
    preprocess inventory & catalog, assign picks, merge alias stations, postprocess
    """
    catfiles       =  ( '../data/phase2001.dat',
                        '../data/phase2002.dat',
                        '../data/phase2003.dat',
                        '../data/phase2004.dat',
                        '../data/phase2005.dat',
                        '../data/phase2006.dat',
                        '../data/phase200701-06.dat' )
                        
    #catfile       =  '../data/phase2001-200706.dat'
    #catfile       = '../data/phase2001.dat'
    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'
    inventoryfile = '../data/inventory2001-200706.pickle'
    
    m0 = memory() 
    
    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"

    m1 = memory( m0 )
    print " memory consumed for inventory: ", m1 / 1000000
    
    pmc.importAliases( aliasfile )
    pmc.preprocessInventory( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after preprocess: inventory has ", len(pmc.stations), " stations"

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
    
        pmc.setCatalog( qpc )
        pmc.preprocessCatalog( qpc, [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )

        print " assigning picks ... "
        pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

        m_curr = memory( m0 )
        print " total memory consumed after assign picks: ", m_curr / 1000000

        del qpc

    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.stations), " stations"
    
    m_end = memory( m0 )
    print " total memory consumed after assign picks: ", m_end / 1000000

    print " pickling PMC inventory ..."
    pickleObj( pmc, inventoryfile )
    print "ready"
    
def test_12():
    """
    create PMC with catalog for 2002 (one chunk)
    import test station list, alias list
    preprocess inventory & catalog, assign picks, merge alias stations, postprocess
    pickle inventory
    """
    catfiles      = ( '../data/phase2002.dat', )
    
    stationfile   = '../data/stationlist.test.dat'
    aliasfile     = '../data/sc_alias.txt'

    inv_pickle    = 'inventory.phase2002.2st.nopickinfo.pickle'

    zmapfile = 'test.zmap.dat'
    #zmapfile = 'test.zmap.special.dat'
    
    m0 = memory() 

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, encoding = 'html' )
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
    
        pmc.setCatalog( qpc )
        pmc.preprocessCatalog( qpc, [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
        print " after catalog preprocess: catalog has ", len( pmc.catalog.eventParameters.event ), " events"

        print " assigning picks ... "
        pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

        m_curr = memory( m0 )
        print " total memory consumed after assign picks: ", m_curr / 1000000

    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.stations), " stations"
    
    m_end = memory( m0 )
    print " total memory consumed after assign picks: ", m_end / 1000000

    print " call fillup() for all stations and delete pickInfo"
    for curr_sta in pmc.stations:
        print " processing station %s" % ( curr_sta.stationCode )
        curr_sta.distribution.fillup()
        del curr_sta.distribution.pickInfo
        
    print " pickling PMC inventory without pickInfo ..."
    pickleObj( pmc, inv_pickle )
    print "ready"

    print " write distro for station 1 as ascii"
    pmc.stations[1].distribution.dump( 'fillup-test-PSP-nopickinfo.dat' )
    print "ready"
    del pmc
    
    print " ----- importing PMC inventory pickle without pick info: ", inv_pickle

    pmc = unpickleObj( inv_pickle )
    print "ready"

    print " write distro for station 1 from pickle"
    pmc.stations[1].distribution.dump( 'fillup-test-PSP-frompickle.dat' )
    print "ready"

def test_13():
    """
    """
    inv_pickle = 'inventory.phase2002-test.pickle'

    check_values = ( (1.0, 10.0),
                     (1.5, 20.0),
                     (1.2, 25.0),
                     
                     (3.0, 40.0),
                     (4.0, 100.0),
                     (4.0, 100.0) )

    pmc = PMC_SCSN()

    print " importing PMC inventory pickle: ", inv_pickle
    pmc = unpickleObj( inv_pickle )
    print "ready"

    for val in check_values:
        value = pmc.stations[0].distribution.calcProbability( val[0], val[1] )
        print " prob. for input %s, %s : %s" % ( val[0], val[1], value )

def test_14():
    """
    create PMC from full station list, no catalog
    write PMC as XML
    """
    stationfile  = '../data/stationlist.dat'
    aliasfile    = '../data/sc_alias.txt'

    xmlout       = 'pmc.xml'
    xmlout_zip   = 'pmc.xml.bz2'
    xmlout_check = 'pmc.check.xml'

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, 1, encoding = 'html' )
    print " PMC_SCSN: PMC inventory has ", len( pmc.stations ), " stations"

    pmc.importAliases( aliasfile )
    pmc.preprocessInventory( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    print " write PMC to XML file: ", xmlout
    pmc.writeXML( xmlout )
    print " ready"

    print " write compressed PMC to XML file: ", xmlout_zip
    pmc.writeXML( xmlout_zip, compression = 'bz2' )
    print " ready"
    
    print " read compressed PMC from XML file: ", xmlout
    pmc_check = PMC_SCSN()
    # pmc_check.setDistroStyle( 1 )
    pmc_check.readXML( xmlout_zip, compression = 'bz2' )
    print " ready"

    print " write check PMC to XML file: ", xmlout_check
    pmc_check.writeXML( xmlout_check  )
    print " ready"

    pmc._preprocessStationScenario()
    print " no of stations = ", len(pmc.stations)

    pmc._selectStationsInOperation( DateTime(2007, 7, 1) )
    print " no of stations = ", len(pmc.stations)

    pmc._selectStationsInOperation( DateTime(1927, 7, 1) )
    print " no of stations = ", len(pmc.stations)
    
    #assert pmc_check == pmc, " PMCs are not equal"

def test_15():
    """
    create PMC with catalog for 2002 (one chunk)
    import test station list, alias list
    preprocess inventory & catalog, assign picks, merge alias stations
    compute distro, do not delete pickInfo
    write as XML
    """
    #catfiles      = ( '../data/phase2002.dat', )
    catfiles      = ( '../data/phase2001-200101.dat', )
    
    stationfile   = '../data/stationlist.test.dat'
    aliasfile     = '../data/sc_alias.txt'

    pmc_xml_file  = 'pmc.phase2002.2st.pickinfo.xml.bz2'

    m0 = memory() 

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, 1, encoding = 'html' )
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

    print " call fillup() for all stations and delete pickInfo"
    for curr_sta in pmc.stations:
        print " processing station %s" % ( curr_sta.stationCode )
        curr_sta.distribution.fillup( pmc._calcMagnitudeFromDistance )
        del curr_sta.distribution.pickInfo
        
    print " write PMC as XML file", pmc_xml_file
    pmc.writeXML( pmc_xml_file, compression = 'bz2'  )
    print "ready"

def test_16():
    """
    create PMC with full stations, import alias list
    preprocess stations: select only Anza network
    import catalog for 2002 (one chunk)
    preprocess catalog,
    assign picks, merge alias stations
    """
    catfiles      = ( '../data/phase2002.dat', )
    
    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'

    use_networks  = ( 'AZ', )
    
    pmc_xml_file  = 'pmc.test_16.xml.bz2'

    m0 = memory() 

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, 1, encoding = 'html' )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"

    m1 = memory( m0 )
    print " memory consumed for PMC inventory: ", m1 / 1000000
    
    pmc.importAliases( aliasfile )
    pmc.preprocessInventory( [ DateTime(2002, 1, 1, 0, 0, 0), DateTime(2003, 1, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    pmc.selectStationsByNetworks( use_networks )

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

def test_17():
    """
    create PMC with full stations, import alias list
    preprocess stations: select only Anza network
    import mini catalog chunk from 2002
    preprocess catalog,
    assign picks, merge alias stations
    ----
    do fillup for all stations (smoothed & not smoothed)
    dump distros
    """
    catfiles      = ( '../data/phase.mini.dat', )

    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'

    use_networks  = ( 'AZ', )

    pmc_xml_file  = 'pmc.test_17.xml.bz2'

    m0 = memory()

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, 3, smooth = True, encoding = 'html' )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"

    m1 = memory( m0 )
    print " memory consumed for PMC inventory: ", m1 / 1000000

    pmc.importAliases( aliasfile )
    pmc.preprocessInventory( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    pmc.selectStationsByNetworks( use_networks )

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

    # over stations - fillup / dump - start with smoothed distro
    print " NOW COMPUTING SMOOTHED DISTROS"
    distrodump_basename = 'fab.pmcdistro.smooth.'
    
    for curr_sta in pmc.stations:
        print " processing station %s" % ( curr_sta.stationCode )
        curr_sta.distribution.setSmooth( True )
        curr_sta.distribution.fillup( pmc._calcMagnitudeFromDistance )
        file_out = distrodump_basename + curr_sta.stationCode + ".dat"
        curr_sta.distribution.dump( file_out )
        print " wrote station distro file %s" % ( file_out )

    # same again for NON-SMOOTHED distro
    print " NOW COMPUTING NON-SMOOTHED DISTROS"
    override = True
    
    distrodump_basename = 'fab.pmcdistro.raw.'
    for curr_sta in pmc.stations:
        print " processing station %s" % ( curr_sta.stationCode )
        curr_sta.distribution.setSmooth( False )
        curr_sta.distribution.fillup( pmc._calcMagnitudeFromDistance, override )
        file_out = distrodump_basename + curr_sta.stationCode + ".dat"
        curr_sta.distribution.dump( file_out )
        print " wrote station distro file %s" % ( file_out )

def test_18():
    """
    test add PMCData to QPGrid
    """

    xml_out = 'testpmcgrid.xml'
    
    pmcdata_str = '<?xml version="1.0" encoding="utf-8"?>\
                   <PMCData>\
                     <probability magnitude="1.3">0.4456</probability>\
                     <probability magnitude="1.4">0.4456</probability>\
                     <probability magnitude="1.5">0.4456</probability>\
                     <mp probability="0.99">2.4</mp>\
                     <mp probability="0.9">2.0</mp>\
                   </PMCData>'

    # get pyRXP object
    tree = pyRXP.Parser().parse( pmcdata_str )

    pmcdata = PMCData()
    pmcdata.fromXML( tree )

    g = QPGrid()
    g.setupBox( -118, -115, 32.5, 34.5, 7.5, 7.5, 0.05, 0.05 )

    tp = ( 'PMCData', 'PMCData', 'element', PMCData, 'complex' )
    g.grid.depthLayer[0].cell[0].addObject( pmcdata, tp )

    g.writeXML( writeQPData( xml_out ) )

def test_19():
    """
    create PMC with full stations, import alias list
    preprocess stations: select only Anza network
    import mini catalog chunk from 2002
    preprocess catalog,
    assign picks, merge alias stations
    ----
    do fillup for all stations (smoothed)
    DO NOT DELETE pickInfo
    write as XML
    write pickle
    """
    catfiles      = ( '../data/phase.mini.dat', )

    stationfile   = '../data/stationlist.dat'
    aliasfile     = '../data/sc_alias.txt'

    use_networks  = ( 'AZ', )

    pmc_xml_file    = 'pmc.anza.test.xml'
    pmc_pickle_file = 'pmc.anza.test.pickle'

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_SCSN( stationfile, 3, smooth = True, encoding = 'html' )
    print " before preprocess: PMC inventory has ", len(pmc.stations), " stations"

    pmc.importAliases( aliasfile )
    pmc.preprocessInventory( [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
    print " after inventory preprocess: inventory has ", len(pmc.stations), " stations"

    pmc.selectStationsByNetworks( use_networks )

    # loop over catalog tranches
    for catfile in catfiles:
        print " ----- importing uncompressed STP phase catalog file: ", catfile

        qpc = QPCatalog()
        qpc.importSTPPhase( getQPDataSource( catfile ) )
        print " catalog has ", len( qpc.eventParameters.event ), " events"

        pmc.preprocessCatalog( qpc, [ DateTime(2001, 1, 1, 0, 0, 0), DateTime(2007, 7, 1, 0, 0, 0) ] )
        print " after catalog preprocess: catalog has ", len( qpc.eventParameters.event ), " events"

        print " assigning picks ... "
        pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )


    pmc.mergeAliasStations()
    print " after assign picks/merge aliases: inventory has ", len(pmc.stations), " stations"

    # over stations - fillup / smoothed distro
    print " NOW COMPUTING SMOOTHED DISTROS"
    
    for curr_sta in pmc.stations:
        print " processing station %s" % ( curr_sta.stationCode )
        curr_sta.distribution.setSmooth( True )
        curr_sta.fillup( pmc._calcMagnitudeFromDistance )

    # write PMC as XML
    fh = writeQPData( pmc_xml_file )
    pmc.writeXML( fh )
    fh.close()

    # write PMC as pickle
    pmc.save( pmc_pickle_file )

def test_20():
    """
    test pickled PMC
    test PMC loaded from XML
    """
    pmc_pickle   = 'pmc.anza.test.pickle'
    pmc_xml_in   = 'pmc.anza.test.xml'

    pmc_frompickle_check = 'pmc.anza.test.frompickle.xml'
    pmc_fromxml_check    = 'pmc.anza.test.fromxml.xml'

    print " ----- get PMC from pickle -----  "
    
    fh = open( pmc_pickle )
    pmc = cPickle.load( fh )
    pmc.writeXML( writeQPData( pmc_frompickle_check ) )

    print " ----- get PMC from XML -----  "
    
    pmc = PMC_SCSN()
    pmc.readXML( getQPDataSource( pmc_xml_in ) )
    pmc.writeXML( writeQPData( pmc_fromxml_check ) )

def test_21():
    """
    test getProbabilityMap
    get PMC from pickle / from XML

    correct grid file to compare with: anza.testgrid.pretty.check.xml
    """
    pmc_xml_in    = 'pmc.anza.test.xml'
    pmc_pickle_in = 'pmc.anza.test.pickle'
    
    grid_out_fromxml     = 'anza.testgrid.fromxml.xml'
    grid_out_frompickle  = 'anza.testgrid.frompickle.xml'

    maggrid      = frange( 0.0, 4.0, 0.1 )
    mp_prob_list = [ 0.9, 0.95, 0.99, 0.999, 0.99999 ]

    print " ----- get PMC from pickle -----  "
    
    fh = open( pmc_pickle_in )
    pmc = cPickle.load( fh )
    print " loaded PMC from pickle: %s stations" % ( len( pmc.stations ) )

    g = QPGrid()
    g.setupBox( -118, -115, 32.5, 34.5, 7.5, 7.5, 0.1, 0.1 )
    #g.setupBox( -118, -115, 32.5, 34.5, 7.5, 7.5, 1.0, 1.0 )

    print " compute prob. map"
    pmc.getProbabilityMap( g, maggrid, Date(2007, 7, 1), mp_prob_list, verbose = False )

    print " write computed grid file ", grid_out_frompickle
    g.writeXML( writeQPData( grid_out_frompickle ) )
    
    print " ----- get PMC from XML -----  "

    pmc = PMC_SCSN()
    pmc.readXML( getQPDataSource( pmc_xml_in ) )
    print " loaded PMC from XML: %s stations" % ( len( pmc.stations ) )
    
    g = QPGrid()
    g.setupBox( -118, -115, 32.5, 34.5, 7.5, 7.5, 0.1, 0.1 )
    #g.setupBox( -118, -115, 32.5, 34.5, 7.5, 7.5, 1.0, 1.0 )
    
    print " compute prob. map"
    pmc.getProbabilityMap( g, maggrid, Date(2007, 7, 1), mp_prob_list, verbose = True  )

    print " write computed grid file ", grid_out_fromxml
    g.writeXML( writeQPData( grid_out_fromxml ) )

def test_22():
    """
    test plotDistribution
    get PMC from XML for POB station (idx 208)
    """
    #pmc_xml_in  = '../data/distro/208.CI.POB.distro.xml'
    #pmc_eps_out = '208.CI.POB.distro.eps'

    pmc_xml_in  = '../data/distro/102.CI.CTW.distro.xml'
    pmc_eps_out = '102.CI.CTW.distro.eps'

    colortable_file = '../data/colortables/rastafari.cpt'
    
    pmc = PMC_SCSN()
    pmc.readXML( getQPDataSource( pmc_xml_in ) )

    pmc.plotDistribution( pmc.stations[0], pmc_eps_out, colortable_file )


    
def test_all():
    test_1()
    test_2()
    test_3()
    test_4()
    test_5()
    test_6()
    test_7()
    test_8()
    test_9()
    test_10()
    test_11()
    test_12()
    test_13()
    test_14()
    test_15()
    test_16()
    test_17()
    
def main():
    
    #test_all()
    #exit()
    
    #test_1()
    #test_2()
    #test_3()
    #test_4()
    #test_5()
    #test_6()
    #test_7()
    #test_8()
    #test_9()
    #test_10()
    #test_11()
    #test_12()
    #test_13()
    #test_14()
    #test_15()
    #test_16()
    #test_17()
    #test_18()
    #test_19()
    #test_20()
    #test_21()
    test_22()
    
main()
