#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

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

sys.path.append('../../../recipes')
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')

from recipes import *

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

jobfile_template = 'distro_job_template.py'
replace_pattern  = r'STATION_INDEX'

runfile          = 'qsub_distro_jobs.sh'

pickinfo_filelist = '../data/pickInfo/pmc.sc.pickInfo.filelist.dat'

# read in template job file as string
jobfile_str = getQPDataSource( jobfile_template ).read()

# open qsub submit file
qsub_fh = writeQPData( runfile )
qsub_fh.write( '#!/bin/sh\n\n' )

# read pickinfo filelist
pickinfo_file_arr = getQPDataSource( pickinfo_filelist ).read().split()
print " i could find %s pickinfo files" % ( len(pickinfo_file_arr) )

# loop over indices of pickinfo_file_arr
for idx in xrange( len( pickinfo_file_arr ) ):
#for idx in xrange( 280, 281, 1 ):
    
    # replace pattern with idx
    new_jobfile_str = re.sub( replace_pattern, str(idx), jobfile_str )

    # create new job filename
    new_jobfile = str(idx) + '.py'

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
