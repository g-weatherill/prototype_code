#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

#
# jma_compute_all.py
# $Id: jma_compute_all.py 147 2008-12-08 19:36:40Z fab $
#

############################################################################
#    Copyright (C) 2008 by Fabian Euchner and Danijel Schorlemmer          #
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
import glob
import re
import random
import datetime

import pickle
import cPickle

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

# this gets the absolute base path  == path where this script resides
# subsequent directory hierarchy will be below this directory
scriptname = os.path.basename( os.path.abspath( sys.argv[0] ) )
basepath   = os.path.dirname( os.path.abspath( sys.argv[0] ) ) # no trailing slash

# directory where catalog (deck) files in bz2 format reside
# deck file names must end with '.Z.bz2'
# works only for deck files with A,B,C part in names (from 1997-10 on)
catdir = os.path.join( basepath, 'data/catalog-bz2/' )

# location of station list file
stationfile = os.path.join( basepath, 'data/jma_station_list_with_network_code.dat' )

# directory where monthly pick info files will be written to
pickInfoDir = os.path.join( basepath, 'data/pickInfo-monthly/' )

# directory where yearly probability distributions will be written to
distroDir   = os.path.join( basepath, 'data/distro-yearly/' )

# start/end years for pickInfo
time_year_start = 2000
time_year_end   = 2008
    
# time intervals for which distros/grid are computed

# this is the one-year interval from which the distros are taken
# there has to be a matching sub-directory of the distro dir.
grid_year_start  = '200704'
grid_year_end    = '200803'

time_intervals = ( ( grid_year_start, grid_year_end ), )

files_per_interval = 12
useDistStyle       = 5

# data for PMC grid

gridfile_prefix  = 'jma'

# compute grid for this date
useDate          = ( 2008, 4, 1 )
useMaxStationCnt = 200                                        # 200

maggrid       = frange( -1.0, 5.0, 0.1 )                       # -1.0, -5.0, 0.1
mp_prob_list  = [ 0.9, 0.95, 0.99, 0.999 ]

# set geographical region for completeness evaluation
geobox = { 'lonmin': 120.0, 'lonmax': 150.0,                  # 120.0-150.0
           'latmin': 20.0,   'latmax': 49.0,                  # 20.0-49.0
           'londelta': 0.1,  'latdelta': 0.1,                 # 0.1-0.1
           'depthmin': 30.0, 'depthmax': 30.0 }

distro_xml_path   = os.path.join( distroDir, ''.join( ( grid_year_start, '-', grid_year_end ) ) )
grid_xml_path     = os.path.join( basepath, 'data/grid', ''.join( ( grid_year_start, '-', grid_year_end ) ) )

# location of combinations pickle
combi_path        = os.path.join( basepath, 'data/combinations/' )
combi_pickle      = os.path.join( combi_path, 'combinations-3-200.numpy.pickle' )


def main():

    compPickInfos()
    compYearlyDistros()
    compGrids()
    
def compPickInfos():

    global basepath
    global catdir
    global stationfile
    global pickInfoDir

    global time_year_start
    global time_year_end
    
    print " ========== computing pickInfos =========="

    catfiles_full = sorted( glob.glob( os.path.join( catdir, '*.Z.bz2' ) ) )
    
    # loop over years
    for time_year in xrange( time_year_start, time_year_end+1 ):
        
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
            print " importing PMC inventory file: %s" % stationfile
            pmc = PMC_JMA( stationfile, 0, smooth = True, encoding = 'html' )
            print " full inventory has %s stations" % len(pmc.stations)
            
            print " before preprocess: PMC inventory has %s stations" % len(pmc.stations)
            pmc.preprocessInventory( [ DateTime(time_year, 1, 1, 0, 0, 0), DateTime(time_year, 12, 31, 23, 59, 59.99) ] )
            print " after inventory preprocess: inventory has %s stations" % len(pmc.stations)
            
            for catfile in catfiles[curr_month]['files']:
                print "   ----- importing catalog chunk from file: ", catfile

                qpc = QPCatalog( idstyle='numeric' )
                qpc.importJMADeck( catfile, compression='bz2', minimumDataset=True )
                print "   catalog chunk has %s events" % qpc.size()

                pmc.preprocessCatalog( qpc, [ DateTime(time_year, 1, 1, 0, 0, 0), DateTime(time_year, 12, 31, 23, 59, 59.99) ] )
                print "   after catalog preprocess: catalog chunk has %s events" % qpc.size()

                print "   assigning picks ... "
                pmc.assignPicks( qpc, reduce_catalog = False, verbose_level = 2 )

                del qpc

            print " now writing station data"
            for curr_sta_idx in reversed( range( len( pmc.stations ) ) ):

                curr_sta = pmc.stations[curr_sta_idx]

                print "   ----- processing station %s: %s %s" % ( curr_sta_idx, curr_sta.locationCode, curr_sta.stationCode )

                # see if directory is already there
                sta_file_dir = os.path.join( pickInfoDir, '.'.join( ( curr_sta.locationCode, curr_sta.stationCode ) ) )

                if not os.path.isdir( sta_file_dir ):
                    print "   creating path %s" % sta_file_dir
                    os.mkdir( sta_file_dir )

                sta_filename_xml = '.'.join( ( str(curr_sta_idx), curr_sta.locationCode, curr_sta.stationCode,
                                              catfiles[curr_month]['name'], 'pickInfo.xml.bz2' ) )
                sta_file_xml     = os.path.join( sta_file_dir, sta_filename_xml )

                print "   ----- writing to XML"
                fh = writeQPData( sta_file_xml, compression='bz2' )
                fh.write( '<?xml version="1.0" encoding="utf-8"?>' )
                fh.write( '<PMC>' )
                curr_sta.toXML( 'PMCStation', fh )
                fh.write( '</PMC>' )
                fh.close()

                # delete station just processed
                del pmc.stations[curr_sta_idx]

def compYearlyDistros():

    global basepath
    global catdir
    global stationfile
    global pickInfoDir
    global distroDir

    global time_year_start
    global time_year_end
    
    global time_intervals

    global files_per_interval
    global useDistStyle

    print " ========== computing yearly distros =========="
    
    print " importing PMC inventory file: ", stationfile
    pmc = PMC_JMA( stationfile )
    print " full inventory has %s stations" % len(pmc.stations)
    
    print " before preprocess: PMC inventory has %s stations" % len(pmc.stations)
    pmc.preprocessInventory( [ DateTime(time_year_start, 1, 1, 0, 0, 0),
                               DateTime(time_year_end, 12, 31, 23, 59, 59.99) ] )
    print " after inventory preprocess: inventory has %s stations" % len(pmc.stations)

    for curr_sta_idx in reversed( xrange( len( pmc.stations ) ) ):

        sta_number = station_chunk[0] + curr_sta_idx
        curr_sta   = pmc.stations[curr_sta_idx]

        print " processing station %s: %s %s" % ( sta_number, curr_sta.locationCode, curr_sta.stationCode )
        
        sta_dirname = '.'.join( ( curr_sta.locationCode, curr_sta.stationCode ) )

        # get list of all filenames
        pickfiles = sorted( glob.glob( os.path.join( pickInfoDir, sta_dirname, '*.pickInfo.xml.bz2' ) ) )
          
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
            distro_dir = os.path.join( distroDir, ''.join( ( curr_ti[0], '-', curr_ti[1] ) ) )

            if not os.path.isdir( distro_dir ):
                os.mkdir( distro_dir )
            
            # get filename for probdistro XML
            distro_filename_xml = ''.join( ( str(sta_number), '.', curr_sta.locationCode, '.', curr_sta.stationCode,
                                             '.', curr_ti[0], '-', curr_ti[1], '.distro.xml.bz2' ) )
            
            distro_file_out = os.path.join( distro_dir, distro_filename_xml )

            # write PMC as XML
            print " writing distro XML file: %s" % distro_file_out
            fh = writeQPData( distro_file_out, compression='bz2' )
            curr_pmc.writeXML( fh )
            fh.close()

            # delete station just processed 
            del curr_pmc

def compGrids():
    """
    Note: output grid file is NOT prettyfied, is b2zipped
          has name <prefix>.YYYY-MM-DDTHH-MM-SS.grid.xml.bz2
    """

    global basepath
    global catdir
    global stationfile
    global pickInfoDir
    global distroDir

    global time_year_start
    global time_year_end
    
    global time_intervals

    global files_per_interval
    global useDistStyle

    global gridfile_prefix
    
    global useDate
    global grid_year_start
    global grid_year_end

    global maggrid
    global mp_prob_list

    global geobox

    global distro_xml_path
    global grid_xml_path

    global combi_path
    global combi_pickle

    print " ========== computing PMC grid =========="
    
    print " ----- get PMC from XML -----  "

    print " load Combinations from pickle"
    combi = unpickleObj( combi_pickle )
        
    pmc = PMC_JMA( combinations = combi )

    # import PMC distro data: loop over distro files for stations
    distrofiles = sorted( glob.glob( os.path.join( distro_xml_path,
                                                 ''.join( ( '*.',
                                                            grid_year_start, '-', grid_year_end,
                                                            '.distro.xml.bz2' ) ) ) ) )

    for curr_sta_file in distrofiles:
        print " importing PMC from XML file: %s" % curr_sta_file
        curr_pmc = PMC_JMA()
        curr_pmc.readXML( getQPDataSource( curr_sta_file, compression='bz2' ) )

        # copy over PMC for station to master PMC object
        pmc.stations.append( curr_pmc.stations[0] )

        del curr_pmc

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
    mapdate_str = mxDateTime2ISO( [ mapdate ], timesepreplacechar = '-', showsecfrac = False )

    # set grid annotation
    g.annotation.setProperty( date = utc(),
                              starttime = mapdate,
                              endtime = mapdate,
                              latmin = geobox['latmin'],
                              latmax = geobox['latmax'],
                              lonmin = geobox['lonmin'],
                              lonmax = geobox['lonmax'] )
                                
    print " compute probability map for Date: %s" % ( mapdate_str )
    pmc.getProbabilityMap( g, maggrid, mapdate, useMaxStationCnt, mp_prob_list, verbose = False )

    # add station data to grid - delete not required fields
    g.stations = []
    for curr_sta in pmc.stations:
        del curr_sta.channels
        del curr_sta.distribution
        g.stations.append( curr_sta )
            
    grid_file_out = os.path.join( grid_xml_path, '.'.join( ( gridfile_prefix, mapdate_str, 'grid.xml' ) ) )

    print " write computed grid file %s" % grid_file_out
    g.writeXML( grid_file_out, compression='gz', prettyPrint = True  )

if __name__ == "__main__":
    main()