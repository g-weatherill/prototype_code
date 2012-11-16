#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# quakepy/pmc/jma/jma_create_full_jobs.py
# $Id: jma_create_full_jobs.py 195 2009-04-06 21:39:28Z fab $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2009 by Fabian Euchner and Danijel Schorlemmer          #
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

jma_create_full_jobs.py

prepare PMC/Japan computations for submission on dynamic cluster,
for several computation dates between startDate, endDate

1) create job file which calls jma_compGrid.py with commandline parameters (shell script)
2) create shell script that submits job files to 'quakeml' queue on cluster
"""

__version__  = '$Id: jma_create_full_jobs.py 195 2009-04-06 21:39:28Z fab $'
__revision__ = '$Revision: 195 $'
__author__   = "Danijel Schorlemmer <ds@usc.edu>, Fabian Euchner <fabian@fabian-euchner.de>"
__license__  = "GPL"

import sys
import os
import re
import datetime

from   mx.DateTime     import DateTime, RelativeDateTime
from   mx.DateTime.ISO import ParseDateTimeUTC

sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')

from QPUtils import *

# on dynamic.usc.edu: $ qsub -l nodes=1:ppn=2 -q <queuename> <scriptname.py>
# to submit to queue: 'quakeml'

# run computation for each month, on first day, with +/- 60 days of catalog data
# last computation: 2008-10-01
# first computation: 2003-07-01

startDate = DateTime( 2008, 10, 1 )
endDate   = DateTime( 2003, 7, 1 )

daysBefore = 60
daysAfter  = 60

workingDir = 'op'

# scriptname: name of the computation program to start
# note: use full path, for qsub system
scriptname = '/home/geovault-00/fabian/prog/pyprog/quakepy/pmc/jma/jma_compGrid.py'

# runfile: create shell script that submits all individual job files to queue 
runfile = 'qsub_jma_full_jobs.sh'

jobfile_header = '#!/bin/sh\n\n'

# open qsub submit file
qsub_fh = writeQPData( runfile )
qsub_fh.write( '#!/bin/sh\n\n' )

useDate = startDate
while ( useDate >= endDate ):

    # run command line, use default station file
    run_command = "python %s -C -X -d %s -t %s -T %s -w %s\n" % ( scriptname,
                                                                  mxDateTime2ISO(useDate, showtime = False),
                                                                  daysBefore,
                                                                  daysAfter,
                                                                  workingDir )

    new_jobfile = mxDateTime2ISO(useDate, showtime = False) + '.sh'
    
    # write jobfile to disk
    fh = writeQPData( new_jobfile )
   
    fh.write( jobfile_header )
    fh.write( run_command )
 
    # make jobfile executable
    os.chmod( new_jobfile, 0755 )

    # add jobfile to qsub file
    qsub_fh.writelines( ( 'qsub -l nodes=1:ppn=1 -q quakeml ', new_jobfile, '\n' ) )
    qsub_fh.writelines( ( 'sleep 3', '\n' ) )
    
    print " wrote jobfile ", new_jobfile

    # set useDate to one month earlier
    useDate = useDate - RelativeDateTime( months = 1, day = 01 )
    
qsub_fh.close()

# make qsub file executable
os.chmod( runfile, 0755 )
