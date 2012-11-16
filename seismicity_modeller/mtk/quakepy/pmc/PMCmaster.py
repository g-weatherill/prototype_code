#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# PMCmaster.py
# $Id: PMCmaster.py 336 2012-07-09 14:10:49Z tyrone $
#

############################################################################
#    Copyright (C) 2008-2011 by Fabian Euchner and Danijel Schorlemmer     #
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

import cPickle
import cProfile
import datetime
import getopt
import glob
import os
import pickle
import random
import re
import sys

from copy import deepcopy

from mx.DateTime     import Date, DateTime, utc
from mx.DateTime.ISO import ParseDateTimeUTC

## set a few paths that do not go into metadata

# this gets the absolute base path  == path where this script resides
# subsequent directory hierarchy will be below this directory
scriptpath = os.path.abspath( sys.argv[0] )
scriptname = os.path.basename( scriptpath )
scriptdir = os.path.dirname( scriptpath )

# the place of all the following QuakePy modules is (usually) here:
sys.path.append(os.path.join( scriptdir, '..'))

from QPCore    import *
from QPUtils   import *
from QPCatalog import *
from QPPolygon import *

from pmc.PMC         import *
from pmc.PMCFactory  import *
from pmc.PMCMetadata import *

from pmc.PMCGridN    import *

# ugly, but still needed by now
import PMC_NZ

# metadata for this computation run
metadata = PMCMetadata()

# get working directory - subsequent directory hierarchy will be below
basepath = os.getcwd()

### location of station list and station alias files
stationdir = os.path.join( basepath, 'data/station/' )

### directory where data from this run will be written to
rundir = os.path.join( basepath, 'data/runs/' )

### location of catalog files
catalogdir = os.path.join( basepath, 'data/catalog/' )

### location of stationpicks-file (one file per station, containing precomputed values for every event) 
stationpicksdir = os.path.join( basepath, 'data/stationpicks/' )

### location of catalog files
colormapdir = os.path.join( scriptdir, '../data/colormaps/pov/' )

# ----------------------------------------------------------------------------------

metadata.smoothDist   = True

# metadata file prefix
metadata_file_prefix = 'pmcinfo'

# verbosity level
verbose_level = 0

# profiling
profileFlag   = False
profileFile   = 'profiler.result'

# ----------------------------------------------------------------------------------

VALID_NETWORK_IDENTIFIERS = ( 'SC', 'JP', 'NZ', 'NC' )

# ----------------------------------------------------------------------------------


def runParts(**kwargs):

    global metadata
    verbose_level = kwargs.get('verbose_level', 0)

    # TODO(tb) grid merge
    if metadata.mergebuckets:
        try:
            mergeGridBuckets( metadata, verbose_level=verbose_level )
        except Exception, e:
            raise ValueError, e
    elif metadata.createStationPicks:
        try:
            createStationPicks( metadata, verbose_level=verbose_level )
        except Exception, e:
            raise ValueError, e
    else:

        if metadata.skipPickInfo is False:
            compPickInfos( verbose_level=verbose_level )

        if metadata.skipProbDistro is False:
            compProbDistros( verbose_level=verbose_level )

        if metadata.skipGrid is False:
            compGrids( verbose_level=verbose_level )

# ---------------------------------------------------------------------------

def PrintHelp():
    global scriptname
    print 'computes PMC grid'
    print 'Usage: %s [OPTION]' % scriptname
    print '  Options'
    print '   -a FILE, --aliasfile=<file>       Alias file name'
    print '   -A DIR, --stationdir=<dir>        Station directory'
    print '   -b VALUE, --currentbucket=<value> current bucket of grid computation (default: 1)'
    print '   -B VALUE, --bucketcount=<value>   number of buckets for grid computation (default: 1)'
    print '   -c DIR, --catalogdir=<dir>        Catalog directory'
    print '   -C, --combipickprob               Compute prob distro directly after pickInfo for each station'
    print '   -d VALUE, --date=<value>          Date for which PMC grid is computed (YYYY-MM-DD)'
    print '   -D, --skipprobdistro              Skip computation of probability distributions'
    print '   -e, --stationpicksdir=<dir>       Directory containing precomputed data for every event per station'
    print '   -E, --createstationpicks          TODO(tb) Create one file per station with precomputed values for every event'
    print '   -f, --plotpickinfo                Plot pickInfos'
    print '   -F, --plotdistro                  Plot probability distros'
    print '   -g VALUE, --griddir=<value>       Directory to place PMC grid file'
    print '   -G, --skipgrid                    Skip computation of PMC Grid'
    print '   -h, --help                        Show this short help'
    print '   -H, --usestationpicks             Using precomputed station events file'
    print '   -l DIR, --colormapdir=<dir>       Color map directory'
    print '   -n VALUE, --network=<value>       Network identifier [SC,JP,NZ,NC]'
    print '   -m VALUE, --subnetwork=<value>    Subnetwork identifier (e.g., for NZ)'
    print '   -M, --mergebuckets                Merge computed grib buckets'
    print '   -O, --overwrite                   Write to already existing dirs, possibly overwrite existing files'
    print '   -P, --skippickinfo                Skip computation of pickInfo'
    print '   -p FILE, --polygonfile=<file>     Polygon boundary file name'
    print '   -r FILE, --regionfile=<file>      Region settings file name'
    print '   -R, --refine                      Refine pick info'
    print '   -s FILE, --stationfile=<file>     Station file name'
    print '   -S, --profile                     Profile run'
    print '   -t VALUE, --timebefore=<value>    Days before date'
    print '   -T VALUE, --timeafter=<value>     Days after date'
    print '   -v VALUE, --verbosity=<value>     Verbosity 0 .. 2 (default: 0)'
    print '   -w VALUE, --workingdir=<value>    Working directory (absolute path)'
    print '   -X, --discardpickinfo             Do not store pickInfo files, only together with -C option'
    print '   -h, --help                        print this information'

# ---------------------------------------------------------------------------

def setUp():

    global rundir
    global scriptname

    global metadata
    
    global stationdir
    global metadata_file_prefix

    global profileFlag
    global profileFile

    global verbose_level

    # command line variables
    in_skipdistro       = False
    in_skipgrid         = False
    in_skippick         = False
    in_combipickprob    = False
    in_discardpickinfo  = False
    in_refinepickinfo   = False
    in_overwrite        = False

    in_plotpickinfo     = False
    in_plotdistro       = False

    in_createstationpicks = False
    in_usestationpicks    = False
    
    in_aliasfile_name   = None
    in_stationdir       = None
    in_polygonfile_name = None
    in_regionfile_name  = None
    in_stationpicksdir  = None
    in_catalogdir       = None
    in_colormapdir      = None
    
    in_date             = None
    in_griddir          = None
    
    in_network          = None
    in_subnetwork       = None
    in_stationfile_name = None
    in_earlier          = None
    in_later            = None
    in_workingdir       = None

    in_current_bucket   = 1
    in_bucket_count     = 1
    in_mergebuckets     = False

    in_verbosity        = 0
               
    # Read commandline arguments
    cmdParams = sys.argv[1:]
    if len( cmdParams ) == 0:
        PrintHelp()
        sys.exit()
            
    opts, args = getopt.gnu_getopt( 
        cmdParams,
        'hSDEGhHPCXfFROMa:A:b:B:c:d:e:g:l:n:m:p:r:s:t:T:v:w:',
        [ 'help',
            'profile',
            'skipprobdistro',
            'createstationpicks',
            'skipgrid',
            'help',
            'usestationpicks',
            'skippickinfo',
            'combipickprob',
            'discardpickinfo',
            'plotpickinfo',
            'plotdistro',
            'refine',
            'overwrite',
            'mergebuckets',
            'aliasfile=',
            'stationdir=',
            'currentbucket=',
            'bucketcount=',
            'catalogdir=',
            'date=',
            'stationpicksdir=',
            'griddir=',
            'colormapdir=',
            'network=',
            'subnetwork=',
            'polygonfile=',
            'regionfile=',
            'stationfile=',
            'timebefore=',
            'timeafter=',
            'verbosity=',
            'workingdir=' ] )

    for option, parameter in opts:

        if option == '-h' or option == '--help':
            PrintHelp()
            sys.exit(0)

        if option == '-S' or option == '--profile':
            profileFlag = True

        if option == '-E' or option == '--createstationpicks':
            in_createstationpicks = True

        if option == '-H' or option == '--usestationpicks':
            in_usestationpicks = True

        if option == '-D' or option == '--skipprobdistro':
            in_skipdistro = True

        if option == '-G' or option == '--skipgrid':
            in_skipgrid = True

        if option == '-P' or option == '--skippickinfo':
            in_skippick = True

        if option == '-C' or option == '--combipickprob':
            in_combipickprob = True

        if option == '-X' or option == '--discardpickinfo':
            in_discardpickinfo = True

        if option == '-f' or option == '--plotpickinfo':
            in_plotpickinfo = True

        if option == '-F' or option == '--plotdistro':
            in_plotdistro = True

        if option == '-R' or option == '--refine':
            in_refinepickinfo = True

        if option == '-O' or option == '--overwrite':
            in_overwrite = True

        if option == '-a' or option == '--aliasfile':
            in_aliasfile_name = parameter

        if option == '-A' or option == '--stationdir':
            in_stationdir = parameter

        if option == '-b' or option == '--currentbucket':
            in_current_bucket = int(parameter)

        if option == '-B' or option == '--bucketcount':
            in_bucket_count = int(parameter)

        if option == '-M' or option == '--mergebuckets':
            in_mergebuckets = True

        if option == '-c' or option == '--catalogdir':
            in_catalogdir = parameter

        if option == '-e' or option == '--stationpicksdir':
            in_stationpicksdir = parameter

        if option == '-d' or option == '--date':
            in_date = parameter

        if option == '-g' or option == '--griddir':
            in_griddir = parameter

        if option == '-l' or option == '--colormapdir':
            in_colormapdir = parameter

        if option == '-n' or option == '--network':
            in_network = parameter

        if option == '-m' or option == '--subnetwork':
            in_subnetwork = parameter

        if option == '-p' or option == '--polygonfile':
            in_polygonfile_name = parameter

        if option == '-r' or option == '--regionfile':
            in_regionfile_name = parameter

        if option == '-s' or option == '--stationfile':
            in_stationfile_name = parameter

        if option == '-t' or option == '--timebefore':
            in_earlier = parameter

        if option == '-T' or option == '--timeafter':
            in_later = parameter

        if option == '-v' or option == '--verbosity':
            in_verbosity = int(parameter)

        if option == '-w' or option == '--workingdir':
            in_workingdir = parameter
            
    verbose_level = in_verbosity

    # check if valid network identifier has been specified
    if in_network is not None and in_network in VALID_NETWORK_IDENTIFIERS:
        metadata.network = in_network
    else:
        error_str = "%s - no valid network identifier has been specified" % ( scriptname )
        raise ValueError, error_str

    ## evaluate directory settings from commandline options

    if in_stationdir is not None:
        metadata.stationDir = in_stationdir
    else:
        metadata.stationDir = stationdir
    
    if in_catalogdir is not None:
        metadata.catalogBaseDir = in_catalogdir
    else:
        metadata.catalogBaseDir = catalogdir

    if in_stationpicksdir is not None:
        metadata.stationPicksDir = in_stationpicksdir
    else:
        metadata.stationPicksDir = stationpicksdir

    if in_colormapdir is not None:
        metadata.colormapDir = in_colormapdir
    else:
        metadata.colormapDir = colormapdir

    # color map file for prob distros
    metadata.probDistroColorMapFile = os.path.join( metadata.colormapDir, 
        'rastafari.pov' )

    # bucket computation:
    if in_bucket_count < 1:
        raise ValueError, "%s - we must define at least one bucket" % ( scriptname )
    elif not 1 <= in_current_bucket <= in_bucket_count:
        raise ValueError, "%s - current bucket must be 1 <= current_bucket (%d) <= bucket_count (%d)" % ( scriptname, in_current_bucket, in_bucket_count )
    else:
        metadata.currentBucket = in_current_bucket
        metadata.bucketCount = in_bucket_count

    metadata.mergebuckets = in_mergebuckets
    metadata.createStationPicks = in_createstationpicks
    metadata.useStationPicks = in_usestationpicks

    # handle overwrite when we have to create stationPicksFiles
    if metadata.createStationPicks and os.path.isdir(metadata.stationPicksDir) and not in_overwrite:
        errMsg = 'createStationsPicksFiles requested, and stationPicksDirectory %s already exists,' % metadata.stationPicksDir
        errMsg += ' but overwrite (-O, --overwrite) is not requested'
        raise IOError, errMsg

    # if we're going to use stationPicksFiles, than at least the directory where we supposed to find them must exist:
    if metadata.useStationPicks and not os.path.isdir(metadata.stationPicksDir):
        errMsg = 'useStationPicks is requested, but the stationPicksDirectory %s does not exist' % metadata.stationPicksDir
        raise IOError, errMsg

    # network-specific settings, hard-coded for the moment
    if metadata.network == 'SC':

        ### NETWORK SETTINGS
        
        metadata.catalogFilesPattern = '*.dat.gz'

        # CHECK:
        metadata.minimumPickCount = 10

        metadata.catalogStartDate = DateTime( 1980, 1, 1 )
        metadata.catalogEndDate   = DateTime( 2008, 1, 1 )

        # distStyle for prob distros
        metadata.useDistStyle = 3

        # stuff for resulting grid XML
        metadata.targetMagArray  = frange( 0.0, 4.0, 0.1 )
        metadata.targetProbArray = [ 0.9, 0.95, 0.99, 0.999 ]

        # coordinate limits: lon (-122.0,-113.5), lat (31.5,38.0)
        metadata.areaPolygon = { 'vertices': ( ( -122.0, 31.5 ),
                                               ( -113.5, 31.5 ),
                                               ( -113.5, 38.0 ),
                                               ( -122.0, 38.0 ),
                                               ( -122.0, 31.5 ) ),
                                               
                                 'metadata': { 'lonDelta'              : 0.1,
                                               'latDelta'              : 0.1,
                                               'lonAlign'              : 0.0,
                                               'latAlign'              : 0.0,
                                               'includePointOnBoundary': True,
                                               'shift'                 : False  }
                               } 

        metadata.areaDepth = 7.5

    elif metadata.network == 'JP':

        ### NETWORK SETTINGS

        # deck file names must end with '.Z.bz2'
        # works only for deck files with A,B,C part in names (starting from October 1997)
        # metadata.catalogFilesPattern = '*.Z.bz2'
        metadata.catalogFilesPattern = '*.bz2'

        metadata.minimumPickCount = 100
        
        # catalog start and end date
        # Note: this is hard-coded at the moment, should be computed on the fly later
        metadata.catalogStartDate = DateTime( 2000, 1, 1 )
        metadata.catalogEndDate   = DateTime( 2009, 2, 1 )

        # distStyle for prob distros
        metadata.useDistStyle = 5

        # stuff for resulting grid XML
        metadata.targetMagArray  = frange( -1.0, 5.0, 0.1 )
        metadata.targetProbArray = [ 0.9, 0.95, 0.99, 0.999 ]

        # coordinate limits: lon (120,150), lat (20,49)
        metadata.areaPolygon = { 'vertices': ( ( 120.0, 20.0 ),
                                               ( 120.0, 36.0 ),
                                               ( 134.0, 49.0 ),
                                               ( 150.0, 49.0 ),
                                               ( 150.0, 20.0 ),
                                               ( 120.0, 20.0 ) ),
                                               
                                 'metadata': { 'lonDelta'              : 0.1,
                                               'latDelta'              : 0.1,
                                               'lonAlign'              : 0.0,
                                               'latAlign'              : 0.0,
                                               'includePointOnBoundary': True,
                                               'shift'                 : False  }
                               } 

        metadata.areaDepth = 30.0

    
    elif metadata.network == 'NZ':

        ### NETWORK SETTINGS

        # check if valid subnetwork identifier has been specified
        if in_subnetwork is not None:
            if in_subnetwork in PMC_NZ.SUBNETWORKS:
                metadata.subnetwork = in_subnetwork
            else:
                error_str = "no valid subnetwork identifier has been specified for network NZ"
                raise ValueError, error_str

        metadata.catalogFilesPattern = 'nz.????-??.xml.gz'

        # CHECK:
        metadata.minimumPickCount = 10

        metadata.catalogStartDate = DateTime( 1950, 1, 1 )
        metadata.catalogEndDate   = DateTime( 2007, 12, 1 )

        # distStyle for prob distros
        metadata.useDistStyle = 3

        # NZ subnetworks are installed from 1987-01-01 on
        if metadata.useDate >= DateTime( 1987, 1, 1 ):
            metadata.separateRegionalNetworks = True
        else:
            # no subnetworks
            metadata.separateRegionalNetworks = False

        # stuff for resulting grid XML
        metadata.targetMagArray  = frange( 0.0, 5.0, 0.1 )
        metadata.targetProbArray = [ 0.9, 0.95, 0.99, 0.999 ]

        # coordinate limits: lon (163.0, 182.0), lat (-50.0, -33.0)
        metadata.areaPolygon = { 'vertices': ( ( 163.0, -50.0 ),
                                               ( 163.0, -33.0 ),
                                               ( 182.0, -33.0 ),
                                               ( 182.0, -50.0 ) ),

                                 'metadata': { 'lonDelta'              : 0.1,
                                               'latDelta'              : 0.1,
                                               'lonAlign'              : 0.0,
                                               'latAlign'              : 0.0,
                                               'includePointOnBoundary': True,
                                               'shift'                 : False  }
                               }

        metadata.areaDepth = 7.5

    elif metadata.network == 'NC':

        ### NETWORK SETTINGS
        
        metadata.catalogFilesPattern = 'Phase*'

        # CHECK:
        metadata.minimumPickCount = 10

        metadata.catalogStartDate = DateTime( 1968, 1, 1 )
        metadata.catalogEndDate   = DateTime( 2012, 1, 1 )

        # distStyle for prob distros
        metadata.useDistStyle = 3

        # stuff for resulting grid XML
        metadata.targetMagArray  = frange( 0.0, 4.0, 0.1 )
        metadata.targetProbArray = [ 0.9, 0.95, 0.99, 0.999 ]

        # coordinate limits: lon (-127.0,-115.0), lat (33.0,44.0)
        metadata.areaPolygon = { 'vertices': ( ( -127.0, 33.0 ),
                                               ( -115.0, 33.0 ),
                                               ( -115.0, 44.0 ),
                                               ( -127.0, 44.0 ),
                                               ( -127.0, 33.0 ) ),
                                               
                                 'metadata': { 'lonDelta'              : 0.1,
                                               'latDelta'              : 0.1,
                                               'lonAlign'              : 0.0,
                                               'latAlign'              : 0.0,
                                               'includePointOnBoundary': True,
                                               'shift'                 : False  }
                               } 

        metadata.areaDepth = 7.5
    else:
        raise ValueError, 'Unkown network identifier'

    # PMC grid file prefix
    metadata.gridFilePrefix = metadata.network.lower()

    ## computation datetime
    if not metadata.createStationPicks:
        try:
            # useDate as mxDateTime instance
            metadata.useDate = ParseDateTimeUTC( in_date )
        except:
            error_str = "%s - illegal date format %s" % ( scriptname, in_date )
            raise ValueError, error_str

        ## surrounding interval of computation datetime

        # get start and end dates
        if in_earlier is not None:
            metadata.earlierDays = int( in_earlier )
            metadata.startDate   = metadata.useDate - DateTimeDelta( float( metadata.earlierDays ) )
        else:
            # use data from beginning of catalog
            metadata.startDate   = metadata.catalogStartDate
            metadata.earlierDays = ( metadata.useDate - metadata.startDate ).days

        if in_earlier is not None:
            metadata.laterDays = int( in_later )
            metadata.endDate   = metadata.useDate + DateTimeDelta( float( metadata.laterDays ) )
        else:
            # use data until end of catalog
            metadata.endDate   = metadata.catalogEndDate
            metadata.laterDays = ( metadata.endDate - metadata.useDate ).days
    else:
        if verbose_level > 1:
            print 'we create stationPicksFiles - no defined time(s)'

    # refine pickInfo using Michael's method
    metadata.refinePickInfo = in_refinepickinfo

    # set working dircetories
    if in_workingdir is not None:

        # working directory has explicitly been given
        metadata.runPath = in_workingdir

        if not os.path.isdir( metadata.runPath ):
            os.makedirs( metadata.runPath )

        # if grid dir is given, create directly below output dir (one level above pickInfo and distro)
        # no date component in dir name
        # this is to collect all grid files of a sequence of runs
        # if no explicit grid dir is given, create std grid dir below working dir
        if in_griddir is not None:
            metadata.gridDir = os.path.join( metadata.runPath, in_griddir )

        # we know where to find pickInfo and probdistro
        # i.e., skipping of computation is possible
        # set flags for skipping computation
        metadata.skipPickInfo   = in_skippick
        metadata.skipProbDistro = in_skipdistro
        metadata.skipGrid       = in_skipgrid

    else:

        # no explicit working directory, use PID
        metadata.runPath = os.path.join( rundir, str( os.getpid() ) )

        if os.path.isdir( metadata.runPath ):
            error_str = "%s - output path already exists, %s" % ( scriptname, metadata.runPath )
            raise IOError, error_str
        else:
            os.makedirs( metadata.runPath )

    # skipping/combining options
    metadata.combinePickInfoProbDistro = in_combipickprob
    metadata.discardPickInfo           = in_discardpickinfo

    # setting combined pickInfo and prob distro switch overrides individual skip switches
    if metadata.combinePickInfoProbDistro is True:
        metadata.skipPickInfo   = False
        metadata.skipProbDistro = True
            
    # plotting of results
    metadata.plotPickInfo   = in_plotpickinfo
    metadata.plotProbDistro = in_plotdistro

    if not metadata.createStationPicks:
        # set instance dir: runPath + date
        instanceDir = os.path.join( metadata.runPath,
                                    mxDateTime2ISO( metadata.useDate, showtime=False ) )

        # if new pickInfos or distros are computed:
        # if instance path exists, stop program execution (don't overwrite)
        # if overwrite option is not explicitly set
        if (     os.path.isdir( instanceDir )
                 and ( metadata.skipPickInfo is False or metadata.skipProbDistro is False )
                 and metadata.mergebuckets is False
                 and in_overwrite is False ):
            error_str = "%s - instance directory already exists, %s" % ( scriptname, instanceDir )
            raise IOError, error_str
        
        # create pickInfo and distro dirs
        metadata.pickInfoDir = os.path.join( instanceDir, 'pickInfo/' )
        metadata.distroDir   = os.path.join( instanceDir, 'distro/' )

        metadata.pickInfoPlotDir = os.path.join( instanceDir, 'pickInfo-plot/' )
        metadata.distroPlotDir   = os.path.join( instanceDir, 'distro-plot/' )
                                 
        # create grid dir, if not explicitly set
        if metadata.gridDir is None:
            metadata.gridDir = os.path.join( instanceDir, 'grid/' )

    # TODO(tb) we don't need a station file wenn metadata.mergebuckets:
    if in_stationfile_name is not None:
        metadata.stationfile = os.path.join( metadata.stationDir, in_stationfile_name )
    else:
        error_str = "%s - no station file name has been specified" % ( scriptname )
        raise ValueError, error_str

    if in_aliasfile_name is not None:
        metadata.aliasfile = os.path.join( metadata.stationDir, in_aliasfile_name )

    if in_polygonfile_name is not None:
        metadata.polygonfile = os.path.join( metadata.stationDir, in_polygonfile_name )

    # set metadata file
    if not metadata.createStationPicks:
        metadataFilename = '.'.join( ( metadata_file_prefix,
                                   mxDateTime2ISO( metadata.useDate, showtime=False ),
                                   '00.dat' ) )
    else:
        metadataFilename = '.'.join( ( metadata_file_prefix, 'CSP', '00.dat') )
    metadata.metadataFile = os.path.join( metadata.runPath, metadataFilename )
    
    # check if metadata file exists, if yes, try another name
    metadataFileCtr = 0
    while ( os.path.isfile( metadata.metadataFile ) ):

        print " - skipping metadata file %s" % metadata.metadataFile

        metadataFileCtr += 1
        if not metadata.createStationPicks:
            metadataFilename = '.'.join( ( metadata_file_prefix,
                                       mxDateTime2ISO( metadata.useDate, showtime=False ),
                                       "%02d" % metadataFileCtr,
                                       'dat' ) )
        else:
            metadataFilename = '.'.join( ( metadata_file_prefix, 'CSP',
                                           '%02d' % metadataFileCtr, 'dat' ) )
        metadata.metadataFile = os.path.join( metadata.runPath, metadataFilename )

    # write metadata to file and screen
    print "write run metadata to file %s" % metadata.metadataFile
    writeMetadata( writeQPData( metadata.metadataFile ) )
    writeMetadata( sys.stdout )

    # set name of profiling file
    if profileFlag is True:
        if not metadata.createStationPicks:
            profileFile = profileFile + '.' + mxDateTime2ISO(
                        metadata.useDate, showtime=False ) + ".%02d" % ( metadataFileCtr ) +'.' + mxDateTime2ISO(
                        metadata.creationDate, showsecfrac=False ) + '.dat'
        else:
            profileFile = profileFile + '.' + 'CSP' + '.%02d' % ( metadataFileCtr ) + '.' + 'CSP' + '.dat'
    if verbose_level > 1:
        print 'setUp() done'

def writeMetadata( stream ):
    """
    metadata is a dictionary with relevant metadata of computation run
    """

    global metadata

    # TODO(tb) add output for:
    # - bucket computation

    stream.write( ''.join( ( "=========================================================================================", '\n' ) ) )
    if not metadata.createStationPicks:
        stream.write( ''.join( ( "PMC computation for network %s, date %s" % (
                  metadata.network, mxDateTime2ISO( metadata.useDate, showtime=False ) ), '\n' ) ) )
        stream.write( ''.join( ( " phase data starts at %s (%s days earlier)" % \
                  ( mxDateTime2ISO( metadata.startDate, showtime=False ),
                    metadata.earlierDays ), '\n' ) ) )
        stream.write( ''.join( ( " phase data ends at %s (%s days later)" % \
                  ( mxDateTime2ISO( metadata.endDate, showtime=False ),
                    metadata.laterDays ), '\n' ) ) )
    else:
        stream.write( ''.join( ( "Create station picks files for network %s" % \
                  metadata.network, '\n' ) ) )
    stream.write( ''.join( ( " run started: %s" % \
                  mxDateTime2ISO( metadata.creationDate, showsecfrac=False ), '\n' ) ) )
    stream.write( ''.join( ( "=========================================================================================", '\n' ) ) )

    if not metadata.createStationPicks:
        if metadata.combinePickInfoProbDistro is True:
            stream.write( ''.join( ( " COMPUTING combined pickInfo and probability distributions", '\n' ) ) )

            if metadata.discardPickInfo is True:
                stream.write( ''.join( ( " -> DISCARDING pickInfo files", '\n' ) ) )
        else:
            if metadata.skipPickInfo is True:
                stream.write( ''.join( ( " SKIPPING pickInfo", '\n' ) ) )
            else:
                stream.write( ''.join( ( " COMPUTING pickInfo", '\n' ) ) )
            if metadata.skipProbDistro is True:
                stream.write( ''.join( ( " SKIPPING probability distributions", '\n' ) ) )
            else:
                stream.write( ''.join( ( " COMPUTING probability distributions", '\n' ) ) )

        if metadata.skipGrid is True:
            stream.write( ''.join( ( " SKIPPING PMC Grid", '\n' ) ) )
        else:
            stream.write( ''.join( ( " COMPUTING PMC Grid", '\n' ) ) )

        if metadata.useStationPicks is True:
            stream.write( ''.join( ( "  using stationPicksFiles in: %s" % metadata.stationPicksDir, '\n' ) ) )
        else:
            stream.write( ''.join( ( "  using catalog data in: %s" % metadata.catalogBaseDir, '\n' ) ) )

    stream.write( ''.join( ( "=========================================================================================", '\n' ) ) )
    stream.write( ''.join( ( " input files:", '\n' ) ) )
    stream.write( ''.join( ( "  station file: %s" % metadata.stationfile, '\n' ) ) )
    stream.write( ''.join( ( "  alias file: %s" % metadata.aliasfile, '\n' ) ) )
    stream.write( ''.join( ( "  polygon file: %s" % metadata.polygonfile, '\n' ) ) )
    stream.write( ''.join( ( "=========================================================================================", '\n' ) ) )
    stream.write( ''.join( ( " operating on paths:", '\n' ) ) )
    if not metadata.createStationPicks:
        stream.write( ''.join( ( "  pickInfo dir: %s" % metadata.pickInfoDir, '\n' ) ) )
        stream.write( ''.join( ( "  prob distro dir: %s" % metadata.distroDir, '\n' ) ) )
        stream.write( ''.join( ( "  grid dir: %s" % metadata.gridDir, '\n' ) ) )
        stream.write( ''.join( ( "=========================================================================================", '\n' ) ) )
        stream.write( ''.join( ( " computation parameters:", '\n' ) ) )
        stream.write( ''.join( ( "  refine pickInfo: %s" % metadata.refinePickInfo, '\n' ) ) )
        stream.write( ''.join( ( "  minimumPickCount (ignore stations with fewer picks): %s" % metadata.minimumPickCount, '\n' ) ) )
    else:
        stream.write( ''.join( ( "  catalog dir: %s" % metadata.catalogBaseDir, '\n' ) ) )
        stream.write( ''.join( ( "  stationPicksFiles dir: %s" % metadata.stationPicksDir, '\n' ) ) )
    stream.write( ''.join( ( "=========================================================================================", '\n' ) ) )
    stream.write( ''.join( ( "the command line was:", '\n' ) ) )
    stream.write( ' '.join( sys.argv ) )
    stream.write( '\n' )
    stream.write( ''.join( ( "=========================================================================================", '\n' ) ) )


def compPickInfos( **kwargs ):

    global metadata

    verbose_level = kwargs.get('verbose_level', 0)

    print " ========== computing pickInfos =========="

    useDateStr = mxDateTime2ISO( metadata.useDate, showtime=False )

    # new inventory object for useDate

    # create output dir
    if not os.path.isdir( metadata.pickInfoDir ):
        if verbose_level > 1:
            print "   creating path %s" % metadata.pickInfoDir
        os.makedirs( metadata.pickInfoDir )
    else:
        if verbose_level > 1:
            print '   path %s already exists' % metadata.pickInfoDir

    if (     ( not os.path.isdir( metadata.pickInfoPlotDir ) )
         and metadata.plotPickInfo is True ):
        if verbose_level > 1:
            print "   creating path %s" % metadata.pickInfoPlotDir
        os.makedirs( metadata.pickInfoPlotDir )
    else:
        if verbose_level > 1:
            print '   path %s already exists or not needed' % metadata.pickInfoPlotDir
        
    if metadata.combinePickInfoProbDistro is True:
        if not os.path.isdir( metadata.distroDir ):
            if verbose_level > 1:
                print "   creating path %s" % metadata.distroDir
            os.makedirs( metadata.distroDir )
        else:
            if verbose_level > 1:
                print '   path %s already exists' % metadata.distroDir

        if (     ( not os.path.isdir( metadata.distroPlotDir ) )
              and metadata.plotProbDistro is True ):
            if verbose_level > 1:
                print "   creating path %s" % metadata.distroPlotDir
            os.makedirs( metadata.distroPlotDir )
        else:
            if verbose_level > 1:
                print '   path %s already exists or not needed' % metadata.distroPlotDir

    print " importing PMC station file: %s" % metadata.stationfile
    pmc = PMCFactory().createPMC( metadata.network, dist_style=0, 
        useSubnetworks=metadata.separateRegionalNetworks, **kwargs )
    pmc.importStations( metadata.stationfile, encoding='html' )
    
    # set timeZoneShift in metadata
    metadata.timeZoneShift = pmc.timeZoneShift
    
    if metadata.aliasfile is not None:
        print " importing PMC alias file: %s" % metadata.aliasfile
        pmc.importAliases( metadata.aliasfile )

    print " before preprocess: PMC inventory has %s stations" % len( 
        pmc.stations )

    pmc.preprocessInventory( [ metadata.startDate, metadata.endDate ] )
    print " after inventory preprocess: inventory has %s stations" %  len( 
        pmc.stations )

    if verbose_level > 1:
        print ' - assign picks for %s stations' % len( pmc.stations )

    # TODO(tb) we need this only if we don't use station picks files
    if metadata.useStationPicks == False:
        pmc.getCatalogFiles( metadata )
        print ' -- using %s catalog files' % len( metadata.catalogFiles )
    else:
        print ' -- using precomputed station picks files in %s' % metadata.stationPicksDir

    pmc.assignPicks( metadata=metadata, compression='bz2', 
        verbose_level=verbose_level )

    if metadata.aliasfile is not None:
        pmc.mergeAliasStations()
        print " after merging aliases: inventory has %s stations" %  len( 
            pmc.stations )


def compProbDistros( **kwargs ):

    global metadata

    verbose_level = kwargs.get('verbose_level', 0)

    print " ========== computing probability distros =========="

    # create output dir
    if not os.path.isdir( metadata.distroDir ):
        if verbose_level > 1:
            print "   creating path %s" % metadata.distroDir
        os.makedirs( metadata.distroDir )
    else:
        if verbose_level > 1:
            print '   path %s already exists' % metadata.distroDir

    if (     ( not os.path.isdir( metadata.distroPlotDir ) )
         and metadata.plotProbDistro is True ):
        if verbose_level > 1:
            print "   creating path %s" % metadata.distroPlotDir
        os.makedirs( metadata.distroPlotDir )
    else:
        if verbose_level > 1:
            print '   path %s already exists' % metadata.distroPlotDir

    # get pick info files
    pickinfofiles = sorted( glob.glob( os.path.join( metadata.pickInfoDir, 
        "*.%s" % PICKINFO_FILE_SUFFIX ) ) )
    
    # loop over stations
    for pickinfofile in pickinfofiles:

        print " processing pick info file %s" % ( pickinfofile )
        if verbose_level > 0:
            myStartTime = utc()

        # create PMC for current station
        curr_pmc = PMCFactory().createPMC( metadata.network, 
            dist_style=metadata.useDistStyle, 
            useSubnetworks=metadata.separateRegionalNetworks )
        
        curr_pmc.readXML( getQPDataSource( pickinfofile, compression='bz2' ) )

        # fillup
        curr_pmc.stations[0].distribution.restoreDistStyle( 
            metadata.useDistStyle )
        curr_pmc.stations[0].distribution.setSmooth( metadata.smoothDist )

        print " ----- now calling fillup ..."
        curr_pmc.stations[0].fillup( curr_pmc._calcMagnitudeFromDistance )

        # delete pick info
        if hasattr( curr_pmc.stations[0].distribution, 'pickInfo' ):
            del curr_pmc.stations[0].distribution.pickInfo

        # write distro to disk
        curr_pmc.writeProbDistros( metadata, 
            verbose_level=verbose_level )

        # plot prob distro
        if metadata.plotProbDistro is True:

            if verbose_level > 1:
                print " --- creating plot of prob distro"

            curr_pmc.plotDistribution( curr_pmc.stations[0],
                metadata, colortable=metadata.probDistroColorMapFile )
        
        # delete station just processed
        del curr_pmc
        if verbose_level > 0:
            print ' --- time used: %s' % ( str( utc() - myStartTime ) )

        
def compGrids( **kwargs ):
    """
    Note: output grid file is gzipped
          do not use built-in pretty-print (uses too much memory)
          use 'xmllint --format' for pretty-print
          has name <prefix>.YYYY-MM-DDTHH-MM-SS.xml.gz
    """

    global metadata
    verbose_level = kwargs.get('verbose_level', 0)
    if verbose_level > 1:
        verbose_grid = True
    else:
        verbose_grid = False

    print " ========== computing PMC grid =========="
    
    print " ----- get PMC from XML -----  "

    # create output dir
    if not os.path.isdir( metadata.gridDir ):
        if verbose_level > 1:
            print "   creating path %s" % metadata.gridDir
        os.makedirs( metadata.gridDir )
    else:
        if verbose_level > 1:
            print '   path %s already exists' % metadata.gridDir

    pmc = PMCFactory().createPMC( metadata.network, 
        useSubnetworks=metadata.separateRegionalNetworks )
    
    # import PMC distro data: loop over distro files for stations
    distrofiles = sorted( glob.glob( os.path.join( metadata.distroDir, 
        "*.%s" % DISTRO_FILE_SUFFIX ) ) )
    _fill_stations_from_distrofiles( distrofiles, pmc, metadata )
    print " loaded PMC from XML: %s stations" % len( pmc.stations )

    print " ----- compute QPGrid -----  "

    ## compute grid

    # create polygon object for border of region
    if metadata.polygonfile is None:
        polygon = QPPolygon( metadata.areaPolygon['vertices'] )
    else:
        polygon = QPPolygon( metadata.polygonfile )
        metadata.areaPolygon['vertices'] = polygon.vertices
    
    # result grid
    g = getGridN( metadata, pmc, polygon )
    mapdate_str = mxDateTime2ISO( [ metadata.useDate ], 
        timesepreplacechar='-', showsecfrac=False )

    if pmc.subnetworks is not None:

        # compute grids for local networks

        # this is an 3-dim array which contains the probability vectors from each
        # subnetwork grid as layers in first dimension
        # axis 0 (layers): probabilities per subnetwork, dimension: len(pmc.subnetworks)
        # axis 1 (rows): cells, dimension: grid.cells.shape[0]
        # axis 2 (cols): probability per magnitude bin, dimension: grid.pmcDataProbCols
        
        local_grids = numpy.zeros( 
            ( len(pmc.subnetworks), g.cells.shape[0], g.pmcDataProbCols ), 
            dtype=float )

        for subnetwork_ctr, curr_subnetwork in enumerate(
            pmc.subnetworks.keys() ):

            # copy master PMC object, since inventory preprocess removes
            # stations
            local_pmc = deepcopy( pmc )

            # set trigger condition
            local_pmc.setSubNetwork( curr_subnetwork )
            metadata.subnetwork = curr_subnetwork

            # create temporary grid object for local network
            local_g = getGridN( metadata, local_pmc, polygon )

            print " compute probability map for date: %s, network %s,"\
                " subnet %s, no. of stations: %s, triggering stations: %s" % ( metadata.useDate, 
                metadata.network, curr_subnetwork, 
                len(local_pmc.subnetworks[curr_subnetwork]['stations']),
                local_pmc.stationsRequiredForTriggering )

            local_pmc.getProbabilityMap( local_g, metadata, verbose=verbose_grid )

            # write grid file for subnet 
            grid_add_station_data( local_g, local_pmc )
            write_gridfile( local_g, 
                grid_filename( metadata, mapdate_str, temp=True, 
                    subnet=curr_subnetwork ),
                grid_filename( metadata, mapdate_str, temp=False, 
                    subnet=curr_subnetwork ) )

            # write to array that collects all subnet data
            local_grids[subnetwork_ctr, :, :] = local_g.pmcDataProb.copy()

            del local_pmc
            del local_g

        # for each grid point, select largest probability value of the 
        # subnet data sets
        # get maximum value of stacked prob arrays, in axis 0
        # and assign to pmcDataProb of result grid
        g.pmcDataProb = numpy.nanmax( local_grids, axis=0 )
        del local_grids

        # re-compute array of Mc values per grid cell
        pmc.computeDataMP( g, metadata )

    else:
        print " compute probability map for date: %s, network %s" % ( 
            metadata.useDate, metadata.network)
        pmc.getProbabilityMap( g, metadata, verbose=verbose_grid )

    grid_add_station_data( g, pmc )
    write_gridfile( g, grid_filename( metadata, mapdate_str, temp=True ),
        grid_filename( metadata, mapdate_str ) )


def createStationPicks( metadata, **kwargs ):
    """
    creates stationPicksFiles for the whole catalog timespan
    for every station in the station file. These files are
    quick2read ascii files - one file for each station
    """
    verbose_level = kwargs.get('verbose_level', 0)
    if verbose_level > 0:
        print 'start in createStationPicks(), verbosity: %d' % verbose_level
    if not os.path.isdir(metadata.stationPicksDir):
        # TODO(tb) create directory
        if verbose_level > 1:
            print '  create station picks files directory: %s' % metadata.stationPicksDir
        try:
            os.makedirs( metadata.stationPicksDir )
        except Exception, e:
            print 'Error: ', e
            raise 'Cannot create station picks files directory: %s' % metadata.stationPicksDir
        else:
            if verbose_level > 1:
                print '  ...done'
    else:
        if verbose_level > 1:
            print 'using the already existent directory %s' % metadata.stationPicksDir
            
    pmc = PMCFactory().createPMC( metadata.network, dist_style=0, 
        useSubnetworks=metadata.separateRegionalNetworks, **kwargs )
    pmc.importStations( metadata.stationfile, encoding='html' )
    
    # set timeZoneShift in metadata
    metadata.timeZoneShift = pmc.timeZoneShift
    
    if metadata.aliasfile is not None:
        print " importing PMC alias file: %s" % metadata.aliasfile
        pmc.importAliases( metadata.aliasfile )

    print " before preprocess: PMC inventory has %s stations" % len( 
        pmc.stations )

    pmc.preprocessInventory()
    print " after inventory preprocess: inventory has %s stations" %  len( 
        pmc.stations )
    print ' - assign picks for %s stations' % len( pmc.stations )

    pmc.getCatalogFiles( metadata )

    event_counter = 0
    stationFileHandles = {}

    for catalog_file in metadata.catalogFilesAll:
        print '  processing file: %s' % catalog_file
        if verbose_level > 1:
            catalogStartTime = utc()
        catalog = pmc.getCatalog( catalog_file )
        if verbose_level > 1:
            print '    raw catalog size: %d events' % catalog.size()
        pmc.preprocessCatalog( catalog )
        if verbose_level > 1:
            print '    after preprocessing: %d events' % catalog.size()
        for curr_event_idx, curr_event in enumerate( catalog.eventParameters.event ):
            try:
                magnitude = curr_event.getPreferredMagnitude()
            except Exception, e:
                if verbose_level > 1:
                    print '    ... no preferred magnitude for event (idx: %d) found. Exception: ' % curr_event_idx, e
                continue
            else:
                event_counter += 1
            curr_ori = curr_event.getPreferredOrigin()
            for station_idx, station in enumerate(pmc.stations):
                if stationFileHandles.get(station_idx, None) == None:
                    resultFileName = '%s/picks_%s.%s.dat' % ( metadata.stationPicksDir, station.networkCode, station.stationCode )
                    # os.path.isfile(resultFileName)
                    try:
                        stationFileHandles[station_idx] = open(resultFileName, 'w')
                    except:
                        raise IOError, 'cannot open station result file: %s' % resultFileName
                distanceXYZ = distanceBetweenPointsWGS84( ( curr_ori.latitude.value,
                                                                    curr_ori.longitude.value,
                                                                    curr_ori.depth.value ),
                                                                  ( station.latitude,
                                                                    station.longitude,
                                                                    -(station.elevation/1000.0) ) )
                distanceXY = distanceBetweenPointsWGS84( ( curr_ori.latitude.value,
                                                                    curr_ori.longitude.value,
                                                                    0 ),
                                                                  ( station.latitude,
                                                                    station.longitude,
                                                                    0 ) )
                myData = {
                          'eventId': curr_ori.publicID,
                          'longitude': float(curr_ori.longitude.value),
                          'latitude': float(curr_ori.latitude.value),
                          'distanceXY': float(distanceXY),
                          'distanceXYZ': float(distanceXYZ),
                          'magnitude': float(magnitude.mag.value),
                          'depth': float(curr_ori.depth.value),
                          'stationNetworkCode': station.networkCode,
                          'stationStationCode': station.stationCode,
                          'picked': int(pmc._eventPicked(curr_event, station.networkCode, station.stationCode)),
                          'time': float(curr_ori.time.value.toDecimalYear()),
                          'magnitudeFromDistance': float(pmc._calcMagnitudeFromDistance(distanceXYZ))
                          }
                line = '{time} {magnitude} {magnitudeFromDistance} {picked} {distanceXY} {distanceXYZ} {eventId}\n'.format(**myData)
                stationFileHandles.get(station_idx).write(line)

        del catalog
        if verbose_level > 1:
            print '  catalog file %s done, time used: %s' % ( catalog_file, str( utc() - catalogStartTime ) )
    if verbose_level > 1:
        print ' processed %d events.' % event_counter
    # for fh_idx in xrange( len(stationFileHandles) ):
    for station_idx, station in enumerate(pmc.stations):
        try:
            stationFileHandles[station_idx].close()
        except:
            print 'OOPS - error while closing filehandle for station %s.%s ?!' % ( station.networkCode, station.stationCode )
        else:
            if verbose_level > 1:
                print 'closed station picks file for station %s.%s' % ( station.networkCode, station.stationCode )

    return True

def mergeGridBuckets( metadata, **kwargs ):

    import tempfile
    import gzip

    verbose_level = kwargs.get('verbose_level', 0)

    print " ========== merging grid buckets =========="

    print 'searching in grid directory: %s' % metadata.gridDir

    pmc = PMCFactory().createPMC( metadata.network, 
        useSubnetworks=metadata.separateRegionalNetworks )

    tempMetadata = metadata
    tempMetadata.mergebuckets = False
    mapdate_str = mxDateTime2ISO( [ metadata.useDate ], timesepreplacechar='-', showsecfrac=False )
    for idx in range(1, ( 1 + metadata.bucketCount ) ):
        tempMetadata.currentBucket = idx
        bucketFileName = grid_filename( tempMetadata, mapdate_str )
        if verbose_level > 1:
            print 'loading file: %s' % bucketFileName

        try:
            # this does not work: True == os.access(bucketFileName, os.R_OK)
            fp = gzip.open(bucketFileName, 'rb')
        except:
            raise IOError, 'cannot open bucket no. %d file: %s' % (idx, bucketFileName)
        # ToDo compressed read file
        tempFile = tempfile.TemporaryFile()
        for line in fp:
            tempFile.write(line)
        fp.close()
        print 'Import XML...'
        pmc.readXML(tempfile)
        print '    ...done'
        tempfile.close()

    pmc.writeXML()

    print " ========== DONE: merging grid buckets =========="
    return True
# ----------------------------------------------------------------------------------

def getGridN( metadata, pmc, polygon ):
    """Return grid object."""

    g = PMCGridN( metadata )
    g.setGridParameter( metadata.areaPolygon['metadata'] )
    g.setupPolygon( polygon, metadata )

    # add annotation object to grid
    g.annotation = pmc.annotation
    g.annotation.setProperty( date                 = utc(),
                              starttime            = metadata.useDate,
                              endtime              = metadata.useDate,
                              lonmin               = g.extent['lonMin'],
                              lonmax               = g.extent['lonMax'],
                              latmin               = g.extent['latMin'],
                              latmax               = g.extent['latMax'],
                              boundary             = polygon.vertices,
                              observationStartTime = metadata.startDate,
                              observationEndTime   = metadata.endDate
                            )
    return g

def grid_add_station_data( grid, pmc ):
    """Add station data to grid - delete fields not required."""
    grid.stations = []
    for curr_sta in pmc.stations:
        del curr_sta.distribution
        del curr_sta.channels
        del curr_sta.onTime
        del curr_sta.onTimeRefined
        del curr_sta.interPickTime

        grid.stations.append( curr_sta )

def grid_filename( metadata, date_str, temp=False, subnet=None ):
    """Set filename for grid file."""
    if temp is False:
        suffix = 'xml.gz'
    else:
        suffix = 'short.xml'

    if metadata.mergebuckets or ( metadata.bucketCount == 1):
        suffix = '%s.%s' % ( date_str, suffix)
    else:
        # from math import *
        bucketNumberLength = 1 + int(math.floor(log10(metadata.bucketCount)))
        suffixTemplateString = '%%0%dd.%%0%dd.%%s.%%s' % ( bucketNumberLength, bucketNumberLength)
        suffix = suffixTemplateString % ( metadata.bucketCount, metadata.currentBucket, date_str, suffix)

    if subnet is None:
        filename = "%s.%s" % ( metadata.gridFilePrefix, suffix )
    else:
        filename = "%s.%s.%s" % ( 
            metadata.gridFilePrefix, subnet, suffix )

    return os.path.join( metadata.gridDir, filename )


def write_gridfile( grid, file_temp, file_final ):
    """Write grid object to XML file."""

    #print " write uncompressed, compact grid file %s" % file_temp
    grid.writeXML( file_temp, prettyPrint=False, noStationChildElements=True )
    
    #print " creating pretty-printed, compressed grid file: %s" % file_final

    exec_str = "xmllint --format %s | gzip > %s" % ( file_temp, file_final )
    os.system( exec_str )

    # check for correct execution of pretty-print (resulting file length > 0)
    # remove compact file only if pretty-print succeeded
    try:
        if os.path.getsize( file_final ) > 0:
            # remove compact file
            os.remove( file_temp )
    except Exception:
        # pretty-printing failed, do nothing
        pass

def _fill_stations_from_distrofiles( distrofiles, pmc, metadata ):
    """Read files with distro information and add station objects
    to master PMC object."""

    for curr_sta_file in distrofiles:

        print " importing PMC from XML file: %s" % curr_sta_file

        curr_pmc = PMCFactory().createPMC( metadata.network,
            useSubnetworks=metadata.separateRegionalNetworks )
        curr_pmc.readXML( getQPDataSource( curr_sta_file, compression='bz2' ) )

        # copy over PMC for station to master PMC object
        pmc.stations.append( curr_pmc.stations[0] )
        del curr_pmc


def main(**kwargs):

    # TODO(tb) in **kwargs?
    global scriptname
    global metadata

    global verbose_level
    global profileFlag
    global profileFile

    programStartTime = utc()

    verbose_level = kwargs.get('verbose_level', 0)

    print " +++ program start time: %s" % ( 
        mxDateTime2ISO( programStartTime, showsecfrac=False ) )

    try:
        setUp()
    except Exception, e:
        raise SystemExit, 'ERROR: %s' % e

    if (     ( metadata.combinePickInfoProbDistro is False )
         and ( metadata.skipProbDistro is True )
         and (     metadata.skipPickInfo is False
               and metadata.skipGrid is False ) ):
        print "%s error: skipping only distro creation (-D) does not make "\
            "sense. Combine with switches -P and -G" % ( scriptname )
        sys.exit()

    if (     ( metadata.discardPickInfo is True )
         and ( metadata.combinePickInfoProbDistro is False ) ):
        print "%s error: you cannot discard pickInfo if not combined with "\
            "distro computation. Use switch -C" % ( scriptname )
        sys.exit()

    if profileFlag is False:
        try:
            runParts(verbose_level=verbose_level)
        except Exception, e:
            raise SystemError, e
    else:
        cProfile.run( 'runParts(verbose_level=verbose_level)', profileFile )

    programEndTime = utc()
    print " +++ program end time: %s" % mxDateTime2ISO( 
        programEndTime, showsecfrac=False )
    print " +++ execution time: %s" % str( programEndTime - programStartTime )


if __name__ == '__main__':
    main(verbose_level=verbose_level)