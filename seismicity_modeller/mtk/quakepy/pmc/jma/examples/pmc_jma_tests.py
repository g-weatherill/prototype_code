#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# quakepy/pmc/jma/examples/pmc_jma_tests.py
# $Id: pmc_jma_tests.py 98 2008-07-30 07:12:44Z fab $
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

__version__  = '$Id: pmc_jma_tests.py 98 2008-07-30 07:12:44Z fab $'
__revision__ = '$Revision: 98 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import random
import datetime
import glob

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
from pmc.PMC_JMA  import *
from pmc.PMCData  import *

from QPGrid           import *
from qpplot           import *
from qpseismicityplot import *

def test_1():
    """ 
    read PMC from XML
    compute distro for first station, delete pickInfo
    write as XML
    plot distro
    """

    pmc_xml_file    = '../data/pickInfo/pmc.jma.2000-200003.pickInfo.xml.bz2'
    pmc_distro_file = '../data/distro/pmc.jma.00.JMA.ABASH2.2000-200003.distro.xml.bz2'

    colortable_file = '../../scsn/data/colortables/rastafari.cpt'
    pmc_eps_out     = '00.JMA.ABASH2.2000-200003.distro.eps'
    
    pmc = PMC_JMA()
    pmc.readXML( getQPDataSource( pmc_xml_file, compression='bz2' ) )

    pmc.setDistroStyle( 4 )

    print " call fillup() for one station and delete pickInfo"
    curr_sta = pmc.stations[0]
    
    print " processing station %s" % ( curr_sta.stationCode )
    
    curr_sta.distribution.fillup( pmc._calcMagnitudeFromDistance )
    del curr_sta.distribution.pickInfo
        
    print " write PMC as XML file", pmc_distro_file
    pmc.writeXML( writeQPData( pmc_distro_file, compression='bz2' ) )
    print "ready"

    pmc.plotDistribution( pmc.stations[0], pmc_eps_out, colortable_file )

def test_2():
    """ 
    check reading of catalog with minimal dataset and numeric ids
    """
    catdir   = '../data/catalog-bz2/'
    catfiles = ( catdir + '200007_A.deck.Z.bz2', )

    cat_xml_file = 'jma.catalog.minimal.qml'

    m0 = memory()
    
    qpc = QPCatalog( idstyle='numeric' )
    
    # loop over catalog chunk pickles
    for catfile in catfiles:

        print " ----- importing catalog chunk from file: ", catfile

        qpc.importJMADeck( catfile, compression='bz2', minimumDataset=True )
        m_curr = memory( m0 )
        
        print " catalog has now %s events" % qpc.size()
        print " total memory consumed after catalog chunk: ", m_curr / 1000000


    print " number of public QP objects: %s" % QPPublicObject.getPublicObjectCtr()
    
    print " ----- exporting catalog as XML to file: ", cat_xml_file
    qpc.writeXML( cat_xml_file )
    
def test_3():
    """ 
    check if full catalog fits into memory (2000-200803)
    write as XML
    """
    catdir   = '../data/catalog-bz2/'
    catfiles = sort( glob.glob( catdir + '*.Z.bz2' ) )

    cat_xml_file = '../data/catalog-xml/catalog.jma.2000-200803.xml.bz2'

    m0 = memory()
    
    qpc = QPCatalog()
    
    # loop over catalog chunk pickles
    for catfile in catfiles:

        print " ----- importing catalog chunk from file: ", catfile

        qpc.importJMADeck( catfile, compression='bz2', minimumDataset=True )
        m_curr = memory( m0 )
        
        print " catalog has now %s events" % qpc.size()
        print " total memory consumed after catalog chunk: ", m_curr / 1000000


    print " ----- exporting catalog as XML to file: ", cat_xml_file
    qpc.writeXML( cat_xml_file, compression='bz2' )

def test_4():
    """ 
    read catalog w/o phases, in chunks
    write each chunk as ZMAP
    """
    catdir   = '../data/catalog-bz2/'
    catfiles = sort( glob.glob( catdir + '*.Z.bz2' ) )

    cat_zmap_dir = '../data/catalog-zmap/'

    m0 = memory()

    chunk_idx = 0

    # loop over catalog chunk pickles
    for catfile in catfiles:

        chunk_idx += 1
        
        qpc = QPCatalog( idstyle='numeric' )
        
        print " ----- importing catalog chunk from file: ", catfile

        qpc.importJMADeck( catfile, compression='bz2', nopicks=True, minimumDataset=True )
        m_curr = memory( m0 )

        cat_zmap_file = cat_zmap_dir + "catalog.jma.%03d.zmap.dat.bz2" % chunk_idx

        print " ----- exporting catalog chunk as b2zipped ZMAP to file: ", cat_zmap_file
        qpc.exportZMAP( cat_zmap_file, compression='bz2' )
    
        print " chunk has %s events" % qpc.size()
        print " memory consumed for catalog chunk: ", m_curr / 1000000

        del qpc


    
def test_5():
    """
    write station list as ascii ( lon | lat )
    """

    stationfile       = '../data/jma_station_list_with_network_code.dat'
    stationlist_ascii = '../data/jma_station_positions.dat'

    print " importing PMC inventory file: ", stationfile
    pmc = PMC_JMA( stationfile, 0, smooth = True, encoding = 'html' )
    print " full inventory has ", len(pmc.stations), " stations"
    
    # over stations - write station with pickInfo separately as XML
    print " +++++ NOW WRITING STATION POSITIONS +++++ "

    # open stationlist file
    flist = writeQPData( stationlist_ascii )
    
    for curr_sta_idx in xrange( len( pmc.stations ) ):

        curr_sta = pmc.stations[curr_sta_idx]
        flist.write( "%8.4f\t%8.4f\n" % ( curr_sta.longitude, curr_sta.latitude ) )

    flist.close()
    
def main():
    
    #test_1()
    #test_2()
    #test_3()
    test_4()
    #test_5()
    
main()