#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# quakepy/pmc/scsn/create_grid_jobs.py
# $Id: create_grid_jobs.py 66 2008-05-28 12:32:00Z fab $
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

__version__  = '$Id: create_grid_jobs.py 66 2008-05-28 12:32:00Z fab $'
__revision__ = '$Revision: 66 $'
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

#from qpplot           import *
#from qpseismicityplot import *


# auf dynamic: $ qsub -l nodes=1:ppn=2 -q <queuename> <scriptname.py>
# queuename: 'quakeml'
# dann schiebt er's in die queue

# there are 453 stations with pickInfo / inventory originally has 462 (after preprocess)

jobfile_template = 'grid_job_template.py'
replace_pattern  = r'DATE_USED'

runfile          = 'qsub_grid_jobs.sh'

date_filelist    = '../data/grid/pmc.sc.ontimes.dat'

# read in template job file as string
jobfile_str = getQPDataSource( jobfile_template ).read()

# open qsub submit file
qsub_fh = writeQPData( runfile )
qsub_fh.write( '#!/bin/sh\n\n' )

# read timestamps
date_file_arr = getQPDataSource( date_filelist ).read().split()
print " i could find %s dates" % ( len(date_file_arr) )

# loop over indices of date_file_arr
for idx in xrange( len(date_file_arr) ):
#for idx in xrange( 0, 2, 1 ):

    curr_date = date_file_arr[idx]
    
    # create DateTime object
    dt = ParseDateTimeUTC( curr_date )

    replace_str = '(' + str(dt.year) + ',' + str(dt.month) + ',' + str(dt.day) + ')'
    
    # replace pattern with idx
    new_jobfile_str = re.sub( replace_pattern, replace_str, jobfile_str )

    # create new job filename
    # new_jobfile = "%04d" % ( idx ) + '.' + curr_date + '.py'
    new_jobfile = "%04d" % ( idx ) + '.py'
    
    # write new jobfile to disk
    fh = writeQPData( new_jobfile )
    fh.write( new_jobfile_str )
    fh.write( '\n' )
    fh.close()

    # make jobfile executable
    os.chmod( new_jobfile, 0755 )

    # add jobfile to qsub file
    qsub_fh.writelines( ( 'qsub -l nodes=1:ppn=2 -q quakeml ', new_jobfile, '\n' ) )
    qsub_fh.writelines( ( 'sleep 3', '\n' ) )
    
    print " wrote jobfile ", new_jobfile

qsub_fh.close()

# make qsub file executable
os.chmod( runfile, 0755 )
