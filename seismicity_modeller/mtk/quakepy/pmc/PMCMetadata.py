# -*- coding: utf-8 -*-
#
# quakepy/pmc/PMCMetadata.py
# $Id: PMCMetadata.py 336 2012-07-09 14:10:49Z tyrone $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2007-2011 by Fabian Euchner and Danijel Schorlemmer     #
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

__version__  = '$Id: PMCMetadata.py 336 2012-07-09 14:10:49Z tyrone $'
__revision__ = '$Revision: 336 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import os
import glob
import re

from mx.DateTime import DateTimeType, TimeDelta, utc

class PMCMetadata( object ):
    """
    all times used in this class are UTC
    """
    
    def __init__( self, **kwargs ):
        '''
        all the configuration values predefined with some reasonable
        default values
        '''

        # init all variables to None or False
        self.useDate           = None
        self.startDate         = None
        self.endDate           = None

        self.network           = None

        self.timeZoneShift     = None

        self.runPath           = None
        self.pickInfoDir       = None
        self.pickTimeDiffDir   = None
        self.distroDir         = None
        self.gridDir           = None

        self.pickInfoPlotDir   = None
        self.distroPlotDir     = None
        
        self.stationDir              = None
        self.combiDir                = None
        self.colormapDir             = None

        self.catalogBaseDir          = None
        self.catalogFilesPattern     = None
        self.catalogFilesPerInterval = None
        self.stationPicksDir         = None
        
        self.catalogFilesAll     = []
        self.catalogFiles        = []

        self.gridFilePrefix    = None

        self.catalogStartDate    = None
        self.catalogEndDate      = None

        self.stationfile         = None
        self.aliasfile           = None
        self.polygonfile         = None

        self.useDistStyle           = None
        self.smoothDist             = None
        self.useMaxStationCnt       = None
        self.refinePickInfo         = None
        self.minimumPickCount       = None

        self.separateRegionalNetworks = False
        self.subnetwork               = None

        # TODO(fab): do we need this list?
        self.subnetworkCodes          = []

        self.targetMagArray           = []
        self.targetProbArray          = []

        self.areaPolygon              = {}
        self.areaDepth                = []

        self.createStationPicks = False
        self.useStationPicks    = False
        self.skipPickInfo       = False
        self.skipProbDistro     = False
        self.skipGrid           = False

        self.combinePickInfoProbDistro = False
        self.discardPickInfo           = False

        self.plotPickInfo              = False
        self.plotProbDistro            = False

        self.currentBucket             = 1
        self.bucketCount               = 1
        self.mergebuckets              = False

        self.probDistroColorMapFile    = None

        self.creationDate              = utc()


    def getCatalogFilesAll( self ):
        if (     self.catalogBaseDir is not None 
             and len( self.catalogBaseDir ) > 0
             and self.catalogFilesPattern is not None 
             and len( self.catalogFilesPattern ) > 0 ):

            self.catalogFilesAll = sorted( glob.glob( os.path.join( 
                self.catalogBaseDir, self.catalogFilesPattern ) ) )

            if len(self.catalogFilesAll) == 0:
                print 'OOPS: no catalog file in the directory %s matches the pattern "%s"' % (self.catalogBaseDir, self.catalogFilesPattern)
        else:
            print 'OOPS: list of all catalogfiles is empty because of empty catalogBaseDir (%s) and/or catalogFilesPattern (%s)' % (self.catalogBaseDir, self.catalogFilesPattern)

#    # TODO(tb) add all the metadata properties...
#    def __str__(self):
#        me = 'This representation of metadata:\n'
#        me += '\tDate: ', self.useDate.toISO()
#        me += '\tStartDate ', self.startDate.toISO()
#        me += '\tEndDate ', self.endDate.toISO()
#        me += '\n'
#        return me