# -*- coding: utf-8 -*-
#
# quakepy/QPCatalog.py
# $Id: QPCatalog.py 334 2012-06-01 17:46:57Z fab $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2007-2010 by Fabian Euchner and Danijel Schorlemmer     #
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

__version__  = '$Id: QPCatalog.py 334 2012-06-01 17:46:57Z fab $'
__revision__ = '$Revision: 334 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import bz2
import cPickle
import cStringIO
import csv
import datetime
import gzip
import math
import os
import re
import string
import time

import numpy
import pyRXP
import shapely.geometry

from mx.DateTime     import DateTime
from mx.DateTime     import utc
from mx.DateTime     import DateTimeDeltaFromSeconds, DateTimeDelta, DateTimeDeltaFrom, TimeDelta
from mx.DateTime     import RelativeDateTime
from mx.DateTime.ISO import ParseDateTimeUTC

from   xml.sax import saxutils

# internal includes

from EventParameters            import EventParameters
from Event                      import Event
from Origin                     import Origin
from Magnitude                  import Magnitude
from FocalMechanism             import FocalMechanism
from MomentTensor               import MomentTensor
from Tensor                     import Tensor
from DataUsed                   import DataUsed
from SourceTimeFunction         import SourceTimeFunction
from PrincipalAxes              import PrincipalAxes
from Axis                       import Axis
from NodalPlanes                import NodalPlanes
from NodalPlane                 import NodalPlane
from OriginUncertainty          import OriginUncertainty
from StationMagnitude           import StationMagnitude
from StationMagnitudeReference  import StationMagnitudeReference
from Amplitude                  import Amplitude
from Pick                       import Pick
from Arrival                    import Arrival
from CreationInfo               import CreationInfo
from Comment                    import Comment

from Phase                      import Phase
from RealQuantity               import RealQuantity
from IntegerQuantity            import IntegerQuantity
from TimeQuantity               import TimeQuantity
from TimeWindow                 import TimeWindow
from CompositeTime              import CompositeTime
from WaveformStreamID           import WaveformStreamID
from EventDescription           import EventDescription
from OriginQuality              import OriginQuality

from QPCore import *
from QPUtils import *

import QPCatalogCompact
import QPPolygon
import QPGrid

import qpfmd
import cumuldist
import qpplot
import qpseismicityplot

POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)

class QPCatalog( QPObject ):
    """
    QuakePy: QPCatalog 
    represents an earthquake catalog
    holds an object of type EventParameters and provides methods for reading, writing, cutting, ...
    """

    # currently, QuakeML version 1.0.1 is supported
    # the namespace required by the schema, however, has 1.0
    xml_namespace = 'http://quakeml.org/xmlns/quakeml/1.0'
    root_attributes = {}
    
    def __init__( self, input=None, **kwargs ):
        """
        input can either be a iostream like a file handle or StringIO object,
        or a string, which is then interpreted as a filename
        """
        super( QPCatalog, self ).__init__( **kwargs )

        # set publicID style
        if 'idstyle' in kwargs and kwargs['idstyle'] is not None :
            QPPublicObject.setPublicIDStyle( kwargs['idstyle'] )

        # set element axis
        self.setElementAxis( '/quakeml' )

        # add eventParameters explicitly
        # child elements of eventParameters are defined via QPElementList
        self.eventParameters = EventParameters( parentAxis=self.elementAxis,
                                                elementName='eventParameters' )
        
        if input is not None:

            if isinstance( input, basestring ):
                istream = getQPDataSource( input, **kwargs )
            else:
                istream = input
                
            self.readXML( istream )

    def merge( self, T ):
        """
        merge contents of another QPCatalog object with self
        NOTE: loses publicID, creationInfo, and comment of merged catalog
        """
        self.eventParameters.event.extend( T.eventParameters.event )
        
    # -------------------------------------------------------------------------
    
    def __eq__( self, T ):
        """
        compare two catalogs
        call __eq__() method on eventParameter attribute
        """

        if hasattr( self, 'eventParameters' ) and self.eventParameters is not None:
          
            if hasattr( T, 'eventParameters' ) and T.eventParameters is not None:
                if not ( self.eventParameters == T.eventParameters ):
                    return False
            else:
                return False
            
        # if no unequal comparison so far, return True
        return True

    # -------------------------------------------------------------------------
    
    def save( self, filename ):
        """
        write catalog object to cPickle
        """
        fh = writeQPData( filename, binary=True )
        try:
            cPickle.dump( self, fh, 2 )
        except cPickle.PickleError:
            raise cPickle.PickleError, "QPCatalog::save - error pickling catalog"
        fh.close()
    
    # -------------------------------------------------------------------------
    
    def readXML( self, input, **kwargs ):
        """
        read catalog from QuakeML serialization
        """
        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input
                
        # get whole content of stream at once
        lines = istream.read()
        
        # check if it is XML
        if not lines.startswith('<?xml'):
            raise IOError, 'QPCatalog::readXML - input stream is not XML'
            
        tree = pyRXP.Parser().parse( lines )
        
        if tree[POS_TAGNAME] != 'quakeml':
            raise TypeError, 'QPCatalog::readXML - input stream is not QuakeML'

        # save attributes of quakeml element, make shallow copy of dict
        if tree[POS_ATTRS] is not None:
            self.root_attributes = tree[POS_ATTRS].copy()
        
        # look for eventParameters tag in children component of 4-tuple
        for child in tree[POS_CHILDREN]:
            
            # skip whitespace children, grab first eventParameters child
            if child[POS_TAGNAME] == 'eventParameters':
                self.eventParameters.fromXML( child )
                break
                    
        if not hasattr( self, 'eventParameters' ):
            raise TypeError, 'QPCatalog::readXML - no eventParameters tag found'

    def writeXML( self, output, **kwargs ):
        """
        serialize catalog to QuakeML

        kwargs:
            prettyPrint - pretty formatting of output XML
        """
        if isinstance( output, basestring ):
            ostream = writeQPData( output, **kwargs )
        else:
            ostream = output

        if 'prettyPrint' in kwargs and kwargs['prettyPrint'] is False:
            prettyPrint = False
        else:
            prettyPrint = True
            
        if prettyPrint is True:

            # serialize to string stream
            streamSuccess = False
            try:
                curr_stream = cStringIO.StringIO()
                self.toXML( curr_stream )
                streamSuccess = True
            except:
                print "QPCatalog::writeXML - error in self.toXML() to StringIO "

            if streamSuccess is True :
                try:
                    xmlPrettyPrint( curr_stream, ostream )
                    return
                except:
                    print "QPCatalog::writeXML - error in xmlPrettyPrint()"

                

        # write to output stream w/o pretty print
        # fallback if pretty-print has not succeeded
        try:
            self.toXML( ostream )
        except:
            raise IOError, "QPCatalog::writeXML - error in self.toXML()"

    # -------------------------------------------------------------------------
    
    def toXML( self, stream ):

        stream.write( '<?xml version="1.0" encoding="utf-8"?>')

        xmlns_set = False
        stream.write( '<quakeml' )

        # loop over root attributes from input document
        for curr_attr_name in self.root_attributes.keys():

            stream.write( " %s=\"%s\"" % ( curr_attr_name, self.root_attributes[curr_attr_name] ) )

            # check if attribute defines a namespace
            if curr_attr_name.find( 'xmlns' ) >= 0:
                xmlns_set = True

        # if no namespace definition  in input document, use standard QuakeML namespace
        if xmlns_set is False:
            stream.write( " xmlns=\"%s\"" % ( self.xml_namespace ) )

        stream.write( '>' )
        
        try:
            self.eventParameters.toXML( 'eventParameters', stream )
        except:
            print "QPCatalog::toXML - error in eventParameters.toXML()"

        stream.write( '</quakeml>' )

    ## -----------------------------------------------------------------------

    def importZMAP( self, input, **kwargs ):
        """ 
        input ZMAP stream has to provide the first 10 columns as specified below (classical ZMAP format)
        there can be additional columns
        CSEP uses three additional columns with uncertainties

        NOTE: all input columns are floats, including those who are in principle integers
              (month, day, hour, minute)
              class attributes, however, are integers

        extended ZMAP format as used in CSEP (www.cseptesting.org) is supported by using
        keyword argument withUncertainties=True (default: False)

        additional CSEP columns:
        11: horizontal error (km) 12: depth error (km) 13: magnitude error

        ZMAP format
        
         col   value                     type
         ---   -----                     ----
           1   longitude                 float
           2   latitude                  float
           3   decimal year              float
           4   month                     float (!)
           5   day                       float (!)
           6   magnitude                 float
           7   depth, in km              float
           8   hour                      float (!)
           9   minute                    float (!)
          10   second                    float
         ---   [CSEP extension] ----------------
          11   horizontal error, in km   float
          12   depth error, in km        float
          13   magnitude error           float

        correct illegal time components: YES
        """

        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        if 'withUncertainties' in kwargs and kwargs['withUncertainties'] is True:
            withUncertainties = True
        else:
            withUncertainties = False
            
        # over lines (= events) in ZMAP input stream
        for line_ctr, line in enumerate( istream ):

            # line entries need to be separated by whitespace
            # split input line into components 
            zmap_pars = line.strip().split()
            
            if len( zmap_pars ) < 10:
                error_str = "QPCatalog.importZMAP - format error in input file, line %s: %s" % ( line_ctr+1, line )
                raise ValueError, error_str
              
            # create event
            ev = Event()
            ev.add( self.eventParameters, 'event' )
              
            # create origin
            ori = Origin()
            ori.add( ev, 'origin' )
            
            ori.longitude = RealQuantity( float( zmap_pars[0] ) )
            ori.latitude  = RealQuantity( float( zmap_pars[1] ) )
            ori.depth     = RealQuantity( float( zmap_pars[6] ) )
            
            # get time components
            corrected_time = correctedDateTimeFromString( 
                float(zmap_pars[2]), float(zmap_pars[3]), float(zmap_pars[4]), 
                float(zmap_pars[7]), float(zmap_pars[8]), 
                float(zmap_pars[9]) )
            ori.time = TimeQuantity( corrected_time[0] )
                
            if corrected_time[1] is True:
                print "corrected illegal datetime components in line %s:"\
                    "%s, %s" % ( line_ctr+1, zmap_pars[2:5], zmap_pars[7:10] )
            
            # create magnitude
            mag = Magnitude()
            mag.add( ev, 'magnitude' )
            mag.mag = RealQuantity( float( zmap_pars[5] ) )
            mag.setOriginAssociation( ori.publicID )
            
            # set preferred origin and magnitude
            ev.preferredOriginID    = ori.publicID
            ev.preferredMagnitudeID = mag.publicID
            # self.eventParameters.event.append( ev )

            ## add uncertainties, if present

            # input line must have at least 13 whitespace-separated entries
            if ( withUncertainties is True ) and ( len( zmap_pars ) >= 13 ):

                # horizontal error
                try:
                    horizontal_error = float( zmap_pars[10] )

                    ou = OriginUncertainty()
                    ou.horizontalUncertainty = horizontal_error
                    ou.add( ori, 'originUncertainty' )
                except:
                    pass
            
                # depth error
                try:
                    ori.depth.uncertainty = float( zmap_pars[11] )
                except:
                    pass

                # magnitude error
                try:
                    mag.mag.uncertainty = float( zmap_pars[12] )
                except:
                    pass


    def exportZMAP( self, output, **kwargs ):
        """ 
        NOTE: all output columns are floats, including those which are in principle integers
              (month, day, hour, minute)
              class attributes, however, are integers

        extended ZMAP format as used in CSEP (www.cseptesting.org) is supported by using
        keyword argument withUncertainties=True (default: False)
        """
        
        if isinstance( output, basestring ):
            ostream = writeQPData( output, **kwargs )
        else:
            ostream = output

        if 'withUncertainties' in kwargs and kwargs['withUncertainties'] is True:
            withUncertainties = True
        else:
            withUncertainties = False
            
        outstring_arr = []
        for ev in self.eventParameters.event:

            # check if event has preferred origin and coordinates, otherwise skip
            try:
                ori = ev.getPreferredOrigin()
                curr_lon = ori.longitude.value
                curr_lat = ori.latitude.value
            except:
                continue

            try:
                mag = ev.getPreferredMagnitude()
            except:
                mag = None
            
            # if event has no associated magnitude, set magnitude column to 'NaN'
            if ( hasattr( mag, 'mag' ) and mag is not None ):
                mag_str = str( mag.mag.value )
            else:
                mag_str = 'NaN'

            # if origin has no depth, set depth column to 'NaN'
            if hasattr( ori, 'depth' ):
                depth_str = str( ori.depth.value )
            else:
                depth_str = 'NaN'

            if ( withUncertainties is True ):

                hzErrorFound = False
                
                # look if explicit horizontal error is given in OriginUncertainty object
                # this overrides possible separate lat/lon errors
                if len( ori.originUncertainty ) > 0:
                    ou = ori.originUncertainty[0]

                    if hasattr( ou, 'horizontalUncertainty' ):
                        try:
                            horizontal_uncertainty_str = str( ou.horizontalUncertainty )
                            hzErrorFound = True
                        except:
                            pass

                # if no explicit horizontal error is given, compute horizontal error from lat/lon errors
                if hzErrorFound is False:

                    if ( hasattr( ori.longitude, 'uncertainty' ) and hasattr( ori.latitude, 'uncertainty' ) ):

                        try:
                            curr_lon_err = ori.longitude.uncertainty
                            curr_lat_err = ori.latitude.uncertainty
                            horizontal_uncertainty_str = str( math.sqrt( math.pow( curr_lat_err * 111.0, 2 ) +
                                                              math.pow( curr_lon_err * math.cos(curr_lat * math.pi/180.0) *
                                                              111.0, 2 ) ) )
                            hzErrorFound = True
                        except:
                            pass


                if hzErrorFound is False:
                    horizontal_uncertainty_str = 'NaN'

                # depth error
                if ( ( depth_str != 'NaN' ) and hasattr( ori.depth, 'uncertainty' ) and ( ori.depth.uncertainty is not None ) ):
                    depth_uncertainty_str = str( ori.depth.uncertainty )
                else:
                    depth_uncertainty_str = 'NaN'

                # magnitude error
                if ( ( mag_str != 'NaN' ) and hasattr( mag.mag, 'uncertainty' ) and ( mag.mag.uncertainty is not None ) ):
                    magnitude_uncertainty_str = str( mag.mag.uncertainty )
                else:
                    magnitude_uncertainty_str = 'NaN'

                line_arr = ( '%10.6f' % ori.longitude.value,
                             '%10.6f' % ori.latitude.value,
                             '%18.12f' % ori.time.value.toDecimalYear(),
                             str( float( ori.time.value.datetime.month ) ),
                             str( float( ori.time.value.datetime.day ) ),
                             mag_str,
                             depth_str,
                             str( float( ori.time.value.datetime.hour ) ),
                             str( float( ori.time.value.datetime.minute ) ),
                             str( ori.time.value.datetime.second ),
                             horizontal_uncertainty_str,
                             depth_uncertainty_str,
                             magnitude_uncertainty_str )
            else:
                line_arr = ( '%10.6f' % ori.longitude.value,
                             '%10.6f' % ori.latitude.value,
                             '%18.12f' % ori.time.value.toDecimalYear(),
                             str( float( ori.time.value.datetime.month ) ),
                             str( float( ori.time.value.datetime.day ) ),
                             mag_str,
                             depth_str,
                             str( float( ori.time.value.datetime.hour ) ),
                             str( float( ori.time.value.datetime.minute ) ),
                             str( ori.time.value.datetime.second ) )
                        
            ostream.writelines( ( '\t'.join( line_arr ), '\n' ) )

    # -------------------------------------------------------------------------
    
    def importSTPPhase( self, input, **kwargs ):
        """
        import SCSN event/phase data as obtained via STP:
        http://www.data.scec.org/STP/stp.html
      
  14018180 le 2004/01/01,00:28:59.260   34.1630   -116.4237  13.14  1.52  l 1.0
    CI    BLA HHZ --   34.0695  -116.3890  1243.0 P d. w  1.0   10.87   2.760
    CI    BLA HHN --   34.0695  -116.3890  1243.0 S .. e  0.5   10.87   4.929
    CI    RMR EHZ --   34.2128  -116.5763  1663.0 P c. i  1.0   15.09   3.391
    CI    GTM EHZ --   34.2946  -116.3560   836.0 P .. e  0.5   15.89   3.362
    CI    GTM EHZ --   34.2946  -116.3560   836.0 S .. i  0.8   15.89   5.913
    CI    JVA HHZ --   34.3662  -116.6127   904.0 P c. i  1.0   28.49   5.207
    CI    JVA HHE --   34.3662  -116.6127   904.0 S .. e  0.5   28.49   8.937

        correct illegal time components: NO
        
        kwargs: nopicks = True - do not read in pick lines
                authorityID    - authority id used in public ids, default: 'SCSN'
        """

        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'authorityID' in kwargs and isinstance( kwargs['authorityID'], basestring ):
            auth_id = kwargs['authorityID']
        else:
            auth_id = 'SCSN'
            
        event_line  = False
        event_ctr   = 0
        
        for line in istream:
            
            scsn_pars = line.split()
            
            # check type of line
            if len( scsn_pars ) == 9:
                
                ## event line
                event_line = True
                pick_ctr   = 0
                
                # create event
                # TODO: event type
                curr_id = scsn_pars[0]
                ev = Event( ''.join( ( 'smi:', auth_id, '/event/', curr_id ) ) )
                ev.add( self.eventParameters, 'event' )
                # self.eventParameters.event.append( ev )
                
                # create origin
                ori = Origin( ''.join( ( 'smi:', auth_id, '/origin/', curr_id ) ) )
                ori.add( ev, 'origin' )
                
                ori.latitude  = RealQuantity( float( scsn_pars[3] ) )
                ori.longitude = RealQuantity( float( scsn_pars[4] ) )
                ori.depth     = RealQuantity( float( scsn_pars[5] ) )
                
                # get time components
                datetime_str = scsn_pars[2]
    
                year   = datetime_str[0:4]
                month  = datetime_str[5:7]
                day    = datetime_str[8:10]
                hour   = datetime_str[11:13]
                mins   = datetime_str[14:16]
                sec    = datetime_str[17:]
                ori.time = TimeQuantity( ( int(year), int(month), int(day), 
                                           int(hour), int(mins), float(sec) ) )
                
                # create magnitude
                # TODO: magnitude type
                mag = Magnitude( ''.join( ( 'smi:', auth_id, '/magnitude/', curr_id ) ) )
                mag.add( ev, 'magnitude' )
                mag.mag = RealQuantity( float( scsn_pars[6] ) )
                mag.setOriginAssociation( ori.publicID )
                
                # set preferred origin and magnitude
                ev.preferredOriginID    = ori.publicID
                ev.preferredMagnitudeID = mag.publicID
                
                event_ctr += 1
                
            elif len( scsn_pars ) == 13:
                
                # if nopicks set, skip line
                if 'nopicks' in kwargs and kwargs['nopicks'] is True:
                    continue
                    
                # if no event line defined: error
                if not event_line:
                    raise ValueError, "QPCatalog::importSTPPhase - phase line without event"
                
                ## phase line
                pick_ctr    = pick_ctr + 1
                
                # create pick
                curr_pickid = ''.join( ( 'smi:', auth_id, '/event/', curr_id, '/pick/', str(pick_ctr) ) )
                pick = Pick( curr_pickid )
                
                # get pick time (origin time plus time diff in column 13)
                pick.time = TimeQuantity( ori.time.value.datetime + DateTimeDeltaFromSeconds( float(scsn_pars[12]) ) )
                
                # get waveform id
                pick.waveformID = WaveformStreamID( scsn_pars[0], scsn_pars[1], scsn_pars[2] )
                if ( scsn_pars[3] != '--' ):
                    pick.waveformID.locationCode = scsn_pars[3]
                    
                # onset
                if ( scsn_pars[9] == 'e' or scsn_pars[9] == 'E' ):
                    pick.onset = 'emergent'
                elif ( scsn_pars[9] == 'i' or scsn_pars[9] == 'I' ):
                    pick.onset = 'impulsive'
                elif ( scsn_pars[9] == 'w' or scsn_pars[9] == 'W' ):
                    pick.onset = 'questionable'
                
                # TODO: polarity
                
                pick.add( ev, 'pick' )
                
                # create arrival
                arrv = Arrival()
                arrv.pickID   = curr_pickid
                arrv.phase    = Phase( scsn_pars[7] )

                # TODO: caution! what unit has epicentral distance in STP phase format?
                # there is nothing said in the documentation
                # QuakeML <= version 1.0.1 uses km, maybe a transformation from degrees is required
                arrv.distance = float( scsn_pars[11] )
                arrv.add( ori, 'arrival' )
                
            elif len( scsn_pars ) == 0:
                # skip empty line
                continue
            else:
                raise ValueError, "QPCatalog::importSTPPhase - format error in input stream"

    # -------------------------------------------------------------------------
    
    def importCMT( self, input, **kwargs ):
        """
        import data from Global CMT catalog in NDK format:
            http://www.globalcmt.org/
            http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/allorder.ndk_explained

        12345678901234567890123456789012345678901234567890123456789012345678901234567890

      
        PDE  2005/01/01 01:20:05.4  13.78  -88.78 193.1 5.0 0.0 EL SALVADOR             
        C200501010120A   B:  4    4  40 S: 27   33  50 M:  0    0   0 CMT: 1 TRIHD:  0.6
        CENTROID:     -0.3 0.9  13.76 0.06  -89.08 0.09 162.8 12.5 FREE S-20050322125201
        23  0.838 0.201 -0.005 0.231 -0.833 0.270  1.050 0.121 -0.369 0.161  0.044 0.240
        V10   1.581 56  12  -0.537 23 140  -1.044 24 241   1.312   9 29  142 133 72   66
        PDE  2005/01/01 01:42:24.9   7.29   93.92  30.0 5.1 0.0 NICOBAR ISLANDS, INDIA R
        C200501010142A   B: 17   27  40 S: 41   58  50 M:  0    0   0 CMT: 1 TRIHD:  0.7
        CENTROID:     -1.1 0.8   7.24 0.04   93.96 0.04  12.0  0.0 BDY  S-20050322125628
        23 -1.310 0.212  2.320 0.166 -1.010 0.241  0.013 0.535 -2.570 0.668  1.780 0.151
        V10   3.376 16 149   0.611 43  44  -3.987 43 254   3.681 282 48  -23  28 73 -136

        correct illegal time components: YES (correct seconds=60.0 by hand)
        """

        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input
            
        self.LINES_PER_EVENT = 5

        line_ctr    = 0
        target_line = 0
        
        for line in istream:
            # line  event_ctr  ev_line_ctr
            # 0     0          0
            # 1     0          1
            # 2     0          2
            # 3     0          3
            # 4     0          4
            # 5     1          0
            # 6     1          1
            # 7     1          2
            # 8     1          3
            # 9     1          4
            # 10    2          0
            
            # current event counter: use integer division
            # current line in event = line_ctr modulo 5
            ( event_ctr, ev_line_ctr ) = divmod( int(line_ctr), int(self.LINES_PER_EVENT) )

            # print " line %s, event %s, event_line %s " % ( str(line_ctr), str(event_ctr), str(ev_line_ctr) )
            
            # truncate line to 80 chars (remove trailing newline chars)
            line = line[0:80]
                
            if line_ctr < target_line:
                line_ctr   = line_ctr+1
                continue

            if ev_line_ctr == 0:

                # First line: Hypocenter line
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                
                # PDE  2005/01/01 01:20:05.4  13.78  -88.78 193.1 5.0 0.0 EL SALVADOR
                #[1-4]   Hypocenter reference catalog (e.g., PDE for USGS location, ISC for
                        #ISC catalog, SWE for surface-wave location, [Ekstrom, BSSA, 2006])
                #[6-15]  Date of reference event
                #[17-26] Time of reference event
                #[28-33] Latitude
                #[35-41] Longitude
                #[43-47] Depth
                #[49-55] Reported magnitudes, usually mb and MS
                #[57-80] Geographical location (24 characters)

                try:
                    ref_catalog       = line[0:4].strip()
                    curr_date_str     = line[5:15].strip()
                    curr_time_str     = line[15:26].strip()
                    curr_lat          = float( line[26:33].strip() )
                    curr_lon          = float( line[33:41].strip() )
                    curr_depth        = float( line[41:47].strip() )
                    curr_mag_str      = line[47:55].strip()
                    curr_location_str = line[55:80].strip()
                except:
                    print " importCMT - error in hypocenter input line %s: %s" % ( line_ctr, line )
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue
                
                ev = Event()
                ev.add( self.eventParameters, 'event' )

                # check if curr_location_str is set (is sometimes missing in NDK file)
                if curr_location_str != '':

                    # sanitize curr_location_str for XML
                    #curr_location_str_unicode = unicode( curr_location_str, 'ascii' )
                    
                    #curr_location_str_xml = saxutils.escape(curr_location_str_unicode).encode('UTF-8')
                    curr_location_str_xml = saxutils.escape(curr_location_str)
                    
                    descr = EventDescription( curr_location_str_xml, 'region name' )
                    ev.description.append( descr )

                # self.eventParameters.event.append( ev )

                ori = Origin()
                ori.add( ev, 'origin' )
                
                # write reference catalog code to origin comment, if not empty string
                # this may be changed in a later version of QuakeML
                if len( ref_catalog ) > 0:
                    ct = Comment( ''.join( ( 'CMT:catalog=', ref_catalog ) ) )
                    ori.comment.append( ct )
                
                ori.latitude  = RealQuantity( curr_lat )
                ori.longitude = RealQuantity( curr_lon )
                ori.depth     = RealQuantity( curr_depth )
                
                # get time components
                try:
                    year   = int( curr_date_str[0:4] )
                    month  = int( curr_date_str[5:7] )
                    day    = int( curr_date_str[8:10] )
                    hour   = int( curr_time_str[0:2] )
                    mins   = int( curr_time_str[3:5] )
                    sec    = float( curr_time_str[6:] )
                except:
                    print " importCMT - date/time error in input line ", line_ctr, "  ", curr_time_str
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue

                # look if seconds are 60.0 (can occur in NDK format)
                # - in that case, set seconds to 59, create time object, and add one second
                add_second = False
                if sec == 60.0:
                    sec = 59.0
                    add_second = True
                    
                try:
                    ori.time = TimeQuantity( ( year, month, day, hour, mins, sec ) )
                except ValueError:
                    print " importCMT - TimeQuantity error in input line ", line_ctr, " ", curr_time_str, " seconds: ", sec
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue

                # if required, add one second to ori.time
                if add_second:
                    ori.time = TimeQuantity( ori.time.value.datetime + DateTimeDeltaFromSeconds( 1.0 ) )
                                           
                # 2 magnitude values, 2nd can be 0.0 (means not set)
                try:
                    curr_mag_arr = [ float( curr_mag_str[0:4].strip() ),
                                     float( curr_mag_str[4:8].strip() ) ]
                except:
                    print " importCMT - magnitude error in input line ", line_ctr, "  ", curr_mag_str
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue
                    
                for curr_mag_idx, curr_mag in enumerate( curr_mag_arr ):

                    # if one of the magnitudes is 0.0, skip
                    if curr_mag > 0.0:
                        mag = Magnitude()
                        mag.mag = RealQuantity( curr_mag )

                        if curr_mag_idx == 0:
                            # first magnitude is always there, is preferred
                            mag.type = 'mb'
                            ev.preferredMagnitudeID = mag.publicID
                        elif curr_mag_idx == 1:
                            mag.type = 'MS'

                        mag.add( ev, 'magnitude' )
                        mag.setOriginAssociation( ori.publicID )

            elif ev_line_ctr == 1:
                # Second line: CMT info (1)
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                
                # C200501010120A   B:  4    4  40 S: 27   33  50 M:  0    0   0 CMT: 1 TRIHD:  0.6
                #[1-16]  CMT event name. This string is a unique CMT-event identifier. Older
                        #events have 8-character names, current ones have 14-character names.
                        #See note (1) below for the naming conventions used.
                #[18-61] Data used in the CMT inversion. Three data types may be used: 
                        #Long-period body waves (B), Intermediate-period surface waves (S),
                        #and long-period mantle waves (M). For each data type, three values
                        #are given: the number of stations used, the number of components 
                        #used, and the shortest period used.
                #[63-68] Type of source inverted for: "CMT: 0" - general moment tensor; 
                        #"CMT: 1" - moment tensor with constraint of zero trace (standard); 
                        #"CMT: 2" - double-couple source.
                #[70-80] Type and duration of moment-rate function assumed in the inversion. 
                        #"TRIHD" indicates a triangular moment-rate function, "BOXHD" indicates
                        #a boxcar moment-rate function. The value given is half the duration
                        #of the moment-rate function. This value is assumed in the inversion,
                        #following a standard scaling relationship (see note (2) below),
                        #and is not derived from the analysis.

                ev_name         = line[0:16].strip()
                data_used_str   = line[17:61].strip()
                source_type_str = line[62:68].strip()
                moment_rate_str = line[69:80].strip()

                fm = FocalMechanism()
                fm.triggeringOriginID = ori.publicID
                fm.add( ev, 'focalMechanism' )
                
                mt = MomentTensor()
                mt.triggeringOriginID = ori.publicID
                mt.cmtName = ev_name

                ## dataUsed
                du_fields = re.match( r'^B:(.+)S:(.+)M:(.+)$', data_used_str )

                # split dataUsed string into components
                du_B_part = du_fields.group(1)
                du_S_part = du_fields.group(2)
                du_M_part = du_fields.group(3)
                
                ## long-period body waves
                try:
                    du_B_1    = int( du_B_part[0:3].strip() )
                    du_B_2    = int( du_B_part[3:8].strip() )
                    du_B_3    = float( du_B_part[8:12].strip() )
                except:
                    print " importCMT - error (dataUsed B) in input line ", line_ctr, " ", du_B_part
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue
                
                # only create object if first value (number of stations) is not zero
                if du_B_1 > 0:
                    du_B = DataUsed( 'long-period body waves', du_B_1, du_B_2, du_B_3 )
                    mt.dataUsed.append( du_B )

                ## intermediate-period surface waves
                try:
                    du_S_1    = int( du_S_part[0:3].strip() )
                    du_S_2    = int( du_S_part[3:8].strip() )
                    du_S_3    = float( du_S_part[8:12].strip() )
                except:
                    print " importCMT - error (dataUsed S) in input line ", line_ctr, " ", du_S_part
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue
                
                # only create object if first value (number of stations) is not zero
                if du_S_1 > 0:
                    du_S = DataUsed( 'intermediate-period surface waves', du_S_1, du_S_2, du_S_3 )
                    mt.dataUsed.append( du_S )

                ## long-period mantle waves
                try:
                    du_M_1    = int( du_M_part[0:3].strip() )
                    du_M_2    = int( du_M_part[3:8].strip() )
                    du_M_3    = float( du_M_part[8:12].strip() )
                except:
                    print " importCMT - error (dataUsed M) in input line ", line_ctr, " ", du_M_part
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue
                
                # only create object if first value (number of stations) is not zero
                if du_M_1 > 0:
                    du_M = DataUsed( 'long-period mantle waves', du_M_1, du_M_2, du_M_3 )
                    mt.dataUsed.append( du_M )

                # method
                method_nr_fields = re.match( r'^CMT:(.+)$', source_type_str )
                method_nr_arr = method_nr_fields.group(1).split()

                if int( method_nr_arr[0] ) == 0:
                    method_nr_str = 'CMT - general moment tensor'
                elif int( method_nr_arr[0] ) == 1:
                    method_nr_str = 'CMT - moment tensor with zero trace'
                elif int( method_nr_arr[0] ) == 2:
                    method_nr_str = 'CMT - double-couple source'
                else:
                    raise ValueError, "importCMT: illegal MomentTensor method"
                mt.method = method_nr_str

                # sourceTimeFunction
                source_time_arr = moment_rate_str.split(':')
                source_time_type = source_time_arr[0].strip()
                source_time_duration = 2.0 * float( source_time_arr[1].strip() )

                if source_time_type == 'TRIHD':
                    source_time_type_str = 'triangle'
                elif source_time_type == 'BOXHD':
                    source_time_type_str = 'box car'
                else:
                    source_time_type_str = 'unknown'
                mt.sourceTimeFunction = SourceTimeFunction( source_time_type_str, source_time_duration )
                mt.add( fm, 'momentTensor' )
                ev.preferredFocalMechanismID = fm.publicID
                
            elif ev_line_ctr == 2:
                # Third line: CMT info (2)
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                
                # CENTROID:     -0.3 0.9  13.76 0.06  -89.08 0.09 162.8 12.5 FREE S-20050322125201
                #[1-58]  Centroid parameters determined in the inversion. Centroid time, given
                        #with respect to the reference time, centroid latitude, centroid
                        #longitude, and centroid depth. The value of each variable is followed
                        #by its estimated standard error. See note (3) below for cases in
                        #which the hypocentral coordinates are held fixed.
                #[60-63] Type of depth. "FREE" indicates that the depth was a result of the
                        #inversion; "FIX " that the depth was fixed and not inverted for;
                        #"BDY " that the depth was fixed based on modeling of broad-band 
                        #P waveforms.
                #[65-80] Timestamp. This 16-character string identifies the type of analysis that
                        #led to the given CMT results and, for recent events, the date and 
                        #time of the analysis. This is useful to distinguish Quick CMTs ("Q-"), 
                        #calculated within hours of an event, from Standard CMTs ("S-"), which 
                        #are calculated later. The format for this string should not be 
                        #considered fixed.
                        
                centroid_par_str = line[0:58].strip()
                depth_type       = line[59:63].strip()
                timestamp_str    = line[64:80].strip()

                try:
                    curr_time      = float( centroid_par_str[9:18].strip() )
                    curr_time_err  = float( centroid_par_str[18:22].strip() )
                    curr_lat       = float( centroid_par_str[22:29].strip() )
                    curr_lat_err   = float( centroid_par_str[29:34].strip() )
                    curr_lon       = float( centroid_par_str[34:42].strip() )
                    curr_lon_err   = float( centroid_par_str[42:47].strip() )
                    curr_depth     = float( centroid_par_str[47:53].strip() )
                    curr_depth_err = float( centroid_par_str[53:58].strip() )
                except:
                    print " importCMT - error (centroid) in input line ", line_ctr, " ", centroid_par_str
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue
                
                # origin from inversion
                ori_inv = Origin()
                mt.derivedOriginID = ori_inv.publicID
                ori_inv.add( ev, 'origin' )

                # set origin from mt inversion as preferred origin
                ev.preferredOriginID = ori_inv.publicID

                ori_inv.latitude  = RealQuantity( curr_lat, curr_lat_err )
                ori_inv.longitude = RealQuantity( curr_lon, curr_lon_err )
                ori_inv.depth     = RealQuantity( curr_depth, curr_depth_err )

                # add seconds of centroid correction to time from triggering origin
                inv_time = ori.time.value.datetime + DateTimeDeltaFromSeconds( curr_time )
                ori_inv.time = TimeQuantity( inv_time,  curr_time_err )
                                      
                # depth type
                if depth_type == 'FREE':
                    depth_type_str = 'from moment tensor inversion'
                elif depth_type == 'FIX':
                    depth_type_str = 'from location'
                elif depth_type == 'BDY':
                    depth_type_str = 'from modeling of broad-band P waveforms'
                else:
                    depth_type_str = 'other'

                ori_inv.depthType = depth_type_str

                # status
                if timestamp_str[0:1] == 'S':
                    mt.status = 'standard CMT solution'
                elif timestamp_str[0:1] == 'Q':
                    mt.status = 'quick CMT solution'

                # MomentTensor/FocalMechanism: creationInfo: creationTime

                # check if timestamp is valid (not valid if it starts with 0)
                if not timestamp_str[2:3] == '0':
                    mt.creationInfo = CreationInfo()
                    mt.creationInfo.creationTime = QPDateTime( ( int(timestamp_str[2:6]),
                                                                 int(timestamp_str[6:8]),
                                                                 int(timestamp_str[8:10]),
                                                                 int(timestamp_str[10:12]),
                                                                 int(timestamp_str[12:14]),
                                                                 float(timestamp_str[14:]) ))

                    fm.creationInfo = CreationInfo()
                    fm.creationInfo.creationTime = QPDateTime( ( int(timestamp_str[2:6]),
                                                                 int(timestamp_str[6:8]),
                                                                 int(timestamp_str[8:10]),
                                                                 int(timestamp_str[10:12]),
                                                                 int(timestamp_str[12:14]),
                                                                 float(timestamp_str[14:]) ) )
                 
            elif ev_line_ctr == 3:
                # Fourth line: CMT info (3)
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                
                # 23 -1.310 0.212  2.320 0.166 -1.010 0.241  0.013 0.535 -2.570 0.668  1.780 0.151
                #[1-2]   The exponent for all following moment values. For example, if the
                        #exponent is given as 24, the moment values that follow, expressed in 
                        #dyne-cm, should be multiplied by 10**24.
                #[3-80]  The six moment-tensor elements: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, where r
                        #is up, t is south, and p is east. See Aki and Richards for conversions
                        #to other coordinate systems. The value of each moment-tensor
                    #element is followed by its estimated standard error. See note (4)
                    #below for cases in which some elements are constrained in the inversion.

                # errors are non-negative, therefore use only 6 digits (values use 7 digits)

                try:

                    moment_exponent = int( line[0:2].strip() )
                    
                    mrr             = float( line[2:9].strip() )
                    mrr_err         = float( line[9:15].strip() )
                    mtt             = float( line[15:22].strip() )
                    mtt_err         = float( line[22:28].strip() )
                    mpp             = float( line[28:35].strip() )
                    mpp_err         = float( line[35:41].strip() )
                    mrt             = float( line[41:48].strip() )
                    mrt_err         = float( line[48:54].strip() )
                    mrp             = float( line[54:61].strip() )
                    mrp_err         = float( line[61:67].strip() )
                    mtp             = float( line[67:74].strip() )
                    mtp_err         = float( line[74:80].strip() )
                    
                except:
                    print " importCMT - error (tensor) in input line ", line_ctr, " ", line
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue

                mt.tensor = Tensor( RealQuantity( exponentialFloatFromString( mrr, moment_exponent ),
                                                  exponentialFloatFromString( mrr_err, moment_exponent ) ),
                                    RealQuantity( exponentialFloatFromString( mtt, moment_exponent ),
                                                  exponentialFloatFromString( mtt_err, moment_exponent ) ),
                                    RealQuantity( exponentialFloatFromString( mpp, moment_exponent ),
                                                  exponentialFloatFromString( mpp_err, moment_exponent ) ),
                                    RealQuantity( exponentialFloatFromString( mrt, moment_exponent ),
                                                  exponentialFloatFromString( mrt_err, moment_exponent ) ),
                                    RealQuantity( exponentialFloatFromString( mrp, moment_exponent ),
                                                  exponentialFloatFromString( mrp_err, moment_exponent ) ),
                                    RealQuantity( exponentialFloatFromString( mtp, moment_exponent ),
                                                  exponentialFloatFromString( mtp_err, moment_exponent ) ) )
                                    
            elif ev_line_ctr == 4:
                # Fifth line: CMT info (4)
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                
                # V10   1.581 56  12  -0.537 23 140  -1.044 24 241   1.312   9 29  142 133 72   66
                #[1-3]   Version code. This three-character string is used to track the version 
                        #of the program that generates the "ndk" file.
                #[5-48]  Moment tensor expressed in its principal-axis system: eigenvalue,
                        #plunge, and azimuth of the three eigenvectors. The eigenvalue should be
                        #multiplied by 10**(exponent) as given on line four.
                #[50-56] Scalar moment, to be multiplied by 10**(exponent) as given on line four.
                #[58-80] Strike, dip, and rake for first nodal plane of the best-double-couple 
                        #mechanism, repeated for the second nodal plane. The angles are defined
                        #as in Aki and Richards.

                # dip is only from 0-90 degrees, therefore uses only 2 digits
                        
                try:
                    version_str         = line[0:3].strip()

                    pa_ei_1             = float( line[3:11].strip() )
                    pa_pl_1             = float( line[11:14].strip() )
                    pa_az_1             = float( line[14:18].strip() )

                    pa_ei_2             = float( line[18:26].strip() )
                    pa_pl_2             = float( line[26:29].strip() )
                    pa_az_2             = float( line[29:33].strip() )

                    pa_ei_3             = float( line[33:41].strip() )
                    pa_pl_3             = float( line[41:44].strip() )
                    pa_az_3             = float( line[44:48].strip() )

                    scalar_moment       = float( line[48:56].strip() )

                    np_st_1             = float( line[56:60].strip() )
                    np_di_1             = float( line[60:63].strip() )
                    np_ra_1             = float( line[63:68].strip() )

                    np_st_2             = float( line[68:72].strip() )
                    np_di_2             = float( line[72:75].strip() )
                    np_ra_2             = float( line[75:80].strip() )
                    
                except:
                    print " importCMT - error (principal axes/nodal planes) in input line ", line_ctr, " ", line
                    target_line = line_ctr + self.LINES_PER_EVENT - ev_line_ctr
                    line_ctr = line_ctr + 1
                    continue

                mt.cmtVersion = version_str

                # scalar moment M0 (in dyne*cm) to moment magnitude Mw:
                # Kanamori (1977): Mw = (2/3)*(log10(M0) - 16.1)
                # see http://www.globalcmt.org/CMTsearch.html#MWnote
                mt.scalarMoment = RealQuantity( exponentialFloatFromString( scalar_moment, moment_exponent ) )

                fm.principalAxes = PrincipalAxes( Axis( RealQuantity( pa_az_1 ),
                                                        RealQuantity( pa_pl_1 ),
                                                        RealQuantity( exponentialFloatFromString( pa_ei_1, moment_exponent ) ) ),
                                                  Axis( RealQuantity( pa_az_2 ),
                                                        RealQuantity( pa_pl_2 ),
                                                        RealQuantity( exponentialFloatFromString( pa_ei_2, moment_exponent ) ) ),
                                                  Axis( RealQuantity( pa_az_3 ),
                                                        RealQuantity( pa_pl_3 ),
                                                        RealQuantity( exponentialFloatFromString( pa_ei_3, moment_exponent ) ) ) )

                fm.nodalPlanes = NodalPlanes( NodalPlane( RealQuantity( np_st_1 ),
                                                          RealQuantity( np_di_1 ),
                                                          RealQuantity( np_ra_1 ) ),
                                              NodalPlane( RealQuantity( np_st_2 ),
                                                          RealQuantity( np_di_2 ),
                                                          RealQuantity( np_ra_2 ) ) )

                
            else:
                raise ValueError, "importCMT: error - never get here"

            # everything fine, increase line_ctr and target_line
            line_ctr    = line_ctr+1
            target_line = target_line+1

    def exportCMT( self, output, **kwargs ):
        """
        output earthquake catalog in NDK format (Global CMT):
            http://www.globalcmt.org/
            http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/allorder.ndk_explained

        12345678901234567890123456789012345678901234567890123456789012345678901234567890

      
        PDE  2005/01/01 01:20:05.4  13.78  -88.78 193.1 5.0 0.0 EL SALVADOR             
        C200501010120A   B:  4    4  40 S: 27   33  50 M:  0    0   0 CMT: 1 TRIHD:  0.6
        CENTROID:     -0.3 0.9  13.76 0.06  -89.08 0.09 162.8 12.5 FREE S-20050322125201
        23  0.838 0.201 -0.005 0.231 -0.833 0.270  1.050 0.121 -0.369 0.161  0.044 0.240
        V10   1.581 56  12  -0.537 23 140  -1.044 24 241   1.312   9 29  142 133 72   66
        PDE  2005/01/01 01:42:24.9   7.29   93.92  30.0 5.1 0.0 NICOBAR ISLANDS, INDIA R
        C200501010142A   B: 17   27  40 S: 41   58  50 M:  0    0   0 CMT: 1 TRIHD:  0.7
        CENTROID:     -1.1 0.8   7.24 0.04   93.96 0.04  12.0  0.0 BDY  S-20050322125628
        23 -1.310 0.212  2.320 0.166 -1.010 0.241  0.013 0.535 -2.570 0.668  1.780 0.151
        V10   3.376 16 149   0.611 43  44  -3.987 43 254   3.681 282 48  -23  28 73 -136
        """

        if isinstance( output, basestring ):
            ostream = writeQPData( output, **kwargs )
        else:
            ostream = output
            
        self.LINES_PER_EVENT = 5

        ## loop over events - use preferred origin
        for curr_ev in self.eventParameters.event:

            # get shortcuts for required objects
            fm = curr_ev.getPreferredFocalMechanism()

            # NOTE: can there be more than one moment tensors?
            mt = fm.momentTensor[0]

            trig_ori    = curr_ev.getOrigin( fm.triggeringOriginID )
            derived_ori = curr_ev.getOrigin( mt.derivedOriginID )

            # First line: Hypocenter line
            # 12345678901234567890123456789012345678901234567890123456789012345678901234567890

            # PDE  2005/01/01 01:20:05.4  13.78  -88.78 193.1 5.0 0.0 EL SALVADOR
            
            #[1-4]   Hypocenter reference catalog (e.g., PDE for USGS location, ISC for
                    #ISC catalog, SWE for surface-wave location, [Ekstrom, BSSA, 2006])
            #[6-15]  Date of reference event
            #[17-26] Time of reference event
            #[28-33] Latitude
            #[35-41] Longitude
            #[43-47] Depth
            #[49-55] Reported magnitudes, usually mb and MS
            #[57-80] Geographical location (24 characters)

            # look if hypocenter reference is there
            # ('micro-format' in comment of triggering origin)
            match_str = r'CMT:catalog=(.*)'
            hyporef = None
            for curr_co in trig_ori.comment:
                matches = re.match( match_str, curr_co.text )
                if ( matches is not None and matches.group(1).strip() > 0 ):
                    hyporef = matches.group(1).strip()
                    break

            if hyporef is not None:
                # add a space char after max. 4 chars of hyporef
                ostream.write( '%-4s ' % string.upper(hyporef) )
            else:
                ostream.write( '     ' )

            # date/time of triggering origin
            # TODO: seconds fraction
            datetime_str = trig_ori.time.value.datetime.strftime( '%Y/%m/%d %H:%M:%S' )
            ostream.write( '%-21s' % datetime_str )

            # lat
            ostream.write( '%7.2f' % trig_ori.latitude.value )
            
            # lon
            ostream.write( '%8.2f' % trig_ori.longitude.value )

            # depth
            if trig_ori.depth.value is not None:
                ostream.write( '%6.1f' % trig_ori.depth.value )
            else:
                ostream.write( '%6.1f' % 0.0 )
            
            # mb
            mags = curr_ev.getMagnitudes( trig_ori )
            
            mb = None
            for curr_mag in mags:

                if curr_mag.type in ( 'mb', 'MB', 'Mb', 'mB' ):
                    mb = curr_mag.mag.value
                    break

            if mb is None:
                mb = 0.0

            ostream.write( '%4.1f' % float(mb) )

            # MS
            ms = None
            for curr_mag in mags:

                if curr_mag.type in ( 'MS', 'ms', 'Ms', 'mS' ):
                    ms = curr_mag.mag.value
                    break

            if ms is None:
                ms = 0.0

            ostream.write( '%4.1f' % float(ms) )
            
            # description (uppercase)
            desc_found = None
            for desc in curr_ev.description:

                if desc.type == 'region name':
                    desc_found = string.upper( desc.text.encode( 'ascii' ) )
                    break

            if desc_found is None:
                desc_found = ''

            # write leading space char before description
            ostream.write( ' %-24s' % desc_found )
            ostream.write( '\n' )
            
            ## 2nd line ()
            # Second line: CMT info (1)
            # 12345678901234567890123456789012345678901234567890123456789012345678901234567890

            # C200501010120A   B:  4    4  40 S: 27   33  50 M:  0    0   0 CMT: 1 TRIHD:  0.6
            #[1-16]  CMT event name. This string is a unique CMT-event identifier. Older
                    #events have 8-character names, current ones have 14-character names.
                    #See note (1) below for the naming conventions used.
            #[18-61] Data used in the CMT inversion. Three data types may be used:
                    #Long-period body waves (B), Intermediate-period surface waves (S),
                    #and long-period mantle waves (M). For each data type, three values
                    #are given: the number of stations used, the number of components
                    #used, and the shortest period used.
            #[63-68] Type of source inverted for: "CMT: 0" - general moment tensor;
                    #"CMT: 1" - moment tensor with constraint of zero trace (standard);
                    #"CMT: 2" - double-couple source.
            #[70-80] Type and duration of moment-rate function assumed in the inversion.
                    #"TRIHD" indicates a triangular moment-rate function, "BOXHD" indicates
                    #a boxcar moment-rate function. The value given is half the duration
                    #of the moment-rate function. This value is assumed in the inversion,
                    #following a standard scaling relationship (see note (2) below),
                    #and is not derived from the analysis.

            ostream.write( '%-16s' % mt.cmtName )

            # fill different types (B, S, M) of dataUsed
            # if not set, set to 0
            ( bStaCnt, bCompCnt, bShortestPeriod ) = ( None, None, None )
            ( sStaCnt, sCompCnt, sShortestPeriod ) = ( None, None, None )
            ( mStaCnt, mCompCnt, mShortestPeriod ) = ( None, None, None )

            for du in mt.dataUsed:
                if string.lower(du.waveType) == 'long-period body waves':
                    # ( val != None and [float(val)] or [None] )[0]
                    bStaCnt         = du.stationCount
                    bCompCnt        = du.componentCount
                    bShortestPeriod = du.shortestPeriod
                if string.lower(du.waveType) == 'intermediate-period surface waves':
                    sStaCnt         = du.stationCount
                    sCompCnt        = du.componentCount
                    sShortestPeriod = du.shortestPeriod
                if string.lower(du.waveType) == 'long-period mantle waves':
                    mStaCnt         = du.stationCount
                    mCompCnt        = du.componentCount
                    mShortestPeriod = du.shortestPeriod

                    
            ostream.write( ' B:' )
            ostream.write( '%3d' % ( bStaCnt is not None and [int(bStaCnt)] or [0] )[0] )
            ostream.write( '%5d' % ( bCompCnt is not None and [int(bCompCnt)] or [0] )[0] )
            ostream.write( '%4d' % ( bShortestPeriod is not None and [int(bShortestPeriod)] or [0] )[0] )

            ostream.write( ' S:' )
            ostream.write( '%3d' % ( sStaCnt is not None and [int(sStaCnt)] or [0] )[0] )
            ostream.write( '%5d' % ( sCompCnt is not None and [int(sCompCnt)] or [0] )[0] )
            ostream.write( '%4d' % ( sShortestPeriod is not None and [int(sShortestPeriod)] or [0] )[0] )

            ostream.write( ' M:' )
            ostream.write( '%3d' % ( mStaCnt is not None and [int(mStaCnt)] or [0] )[0] )
            ostream.write( '%5d' % ( mCompCnt is not None and [int(mCompCnt)] or [0] )[0] )
            ostream.write( '%4d' % ( mShortestPeriod is not None and [int(mShortestPeriod)] or [0] )[0] )

            ostream.write( ' CMT:' )

            # if moment tensor type is not set, set to default (1)
            if ( string.lower(mt.method) == 'cmt - general moment tensor' ):
                mtmethod = 0
            elif ( string.lower(mt.method) == 'cmt - double-couple source' ):
                mtmethod = 2
            else:
                mtmethod = 1

            ostream.write( ' %1d' % mtmethod )

            # if source time function not given, set type to TRI:
            # and value to 0.0
            if ( string.lower(mt.sourceTimeFunction.type) == 'triangle' ):
                stf = ' TRIHD:'
            elif ( string.lower(mt.sourceTimeFunction.type) == 'box car' ):
                stf = ' BOXHD:'
            else:
                stf = ' TRIHD:'

            ostream.write( stf )
            ostream.write( '%5.1f' % ( mt.sourceTimeFunction.duration is not None and [float(0.5*mt.sourceTimeFunction.duration)] or [0] )[0] )

            ostream.write( '\n' )
            
            ## 3rd line ()
            # Third line: CMT info (2)
            # 12345678901234567890123456789012345678901234567890123456789012345678901234567890

            # CENTROID:     -0.3 0.9  13.76 0.06  -89.08 0.09 162.8 12.5 FREE S-20050322125201
            #[1-58]  Centroid parameters determined in the inversion. Centroid time, given
                    #with respect to the reference time, centroid latitude, centroid
                    #longitude, and centroid depth. The value of each variable is followed
                    #by its estimated standard error. See note (3) below for cases in
                    #which the hypocentral coordinates are held fixed.
            #[60-63] Type of depth. "FREE" indicates that the depth was a result of the
                    #inversion; "FIX " that the depth was fixed and not inverted for;
                    #"BDY " that the depth was fixed based on modeling of broad-band
                    #P waveforms.
            #[65-80] Timestamp. This 16-character string identifies the type of analysis that
                    #led to the given CMT results and, for recent events, the date and
                    #time of the analysis. This is useful to distinguish Quick CMTs ("Q-"),
                    #calculated within hours of an event, from Standard CMTs ("S-"), which
                    #are calculated later. The format for this string should not be
                    #considered fixed.

            ostream.write( 'CENTROID:' )

            # time difference
            time_diff = derived_ori.time.value.datetime - trig_ori.time.value.datetime

            ostream.write( '%9.1f' % time_diff )
            ostream.write( '%4.1f' % derived_ori.time.uncertainty )

            # latitude & error
            ostream.write( '%7.2f' % derived_ori.latitude.value )
            ostream.write( '%5.2f' % derived_ori.latitude.uncertainty )

            # longitude & error
            ostream.write( '%8.2f' % derived_ori.longitude.value )
            ostream.write( '%5.2f' % derived_ori.longitude.uncertainty )

            # depth & error
            ostream.write( '%6.1f' % derived_ori.depth.value )
            ostream.write( '%5.1f' % derived_ori.depth.uncertainty )

            # if depth type not given, use FIX as default
            if ( string.lower(derived_ori.depthType) == 'from moment tensor inversion' ):
                dt = ' FREE'
            elif ( string.lower(derived_ori.depthType) == 'from modeling of broad-band p waveforms' ):
                dt = ' BDY '
            else:
                dt = ' FIX '
                
            ostream.write( dt )
            
            # if method not given, use standard 'S-' as default
            if ( string.lower(mt.status) == 'quick cmt solution' ):
                st = ' Q-'
            else:
                st = ' S-'
                
            ostream.write( st )
            if mt.creationInfo.creationTime is not None:
                ostream.write( mt.creationInfo.creationTime.datetime.strftime( '%Y%m%d%H%M%S' ) )
            else:
                ostream.write( '00000000000000' )
                
            ostream.write( '\n' )
            
            ## 4th line ()
            # Fourth line: CMT info (3)
            # 12345678901234567890123456789012345678901234567890123456789012345678901234567890

            # 23 -1.310 0.212  2.320 0.166 -1.010 0.241  0.013 0.535 -2.570 0.668  1.780 0.151
            #[1-2]   The exponent for all following moment values. For example, if the
                    #exponent is given as 24, the moment values that follow, expressed in
                    #dyne-cm, should be multiplied by 10**24.
            #[3-80]  The six moment-tensor elements: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, where r
                    #is up, t is south, and p is east. See Aki and Richards for conversions
                    #to other coordinate systems. The value of each moment-tensor
                #element is followed by its estimated standard error. See note (4)
                #below for cases in which some elements are constrained in the inversion.

            # errors are non-negative, therefore use only 6 digits (values use 7 digits)

            # get exponent from scalar moment value
            scm_mantissa, exponent = normalizeFloat( mt.scalarMoment.value )
            power = 10**exponent
            
            ostream.write( '%02d' % exponent )
            
            ostream.write( '%7.3f' % (mt.tensor.Mrr.value / power) )
            ostream.write( '%6.3f' % (mt.tensor.Mrr.uncertainty / power) )

            ostream.write( '%7.3f' % (mt.tensor.Mtt.value / power) )
            ostream.write( '%6.3f' % (mt.tensor.Mtt.uncertainty / power) )

            ostream.write( '%7.3f' % (mt.tensor.Mpp.value / power) )
            ostream.write( '%6.3f' % (mt.tensor.Mpp.uncertainty / power) )

            ostream.write( '%7.3f' % (mt.tensor.Mrt.value / power) )
            ostream.write( '%6.3f' % (mt.tensor.Mrt.uncertainty / power) )

            ostream.write( '%7.3f' % (mt.tensor.Mrp.value / power) )
            ostream.write( '%6.3f' % (mt.tensor.Mrp.uncertainty / power) )

            ostream.write( '%7.3f' % (mt.tensor.Mtp.value / power) )
            ostream.write( '%6.3f' % (mt.tensor.Mtp.uncertainty / power) )
            
            ostream.write( '\n' )
            
            ## 5th line ()
            # Fifth line: CMT info (4)
            # 12345678901234567890123456789012345678901234567890123456789012345678901234567890

            # V10   1.581 56  12  -0.537 23 140  -1.044 24 241   1.312   9 29  142 133 72   66
            #[1-3]   Version code. This three-character string is used to track the version
                    #of the program that generates the "ndk" file.
            #[5-48]  Moment tensor expressed in its principal-axis system: eigenvalue,
                    #plunge, and azimuth of the three eigenvectors. The eigenvalue should be
                    #multiplied by 10**(exponent) as given on line four.
            #[50-56] Scalar moment, to be multiplied by 10**(exponent) as given on line four.
            #[58-80] Strike, dip, and rake for first nodal plane of the best-double-couple
                    #mechanism, repeated for the second nodal plane. The angles are defined
                    #as in Aki and Richards.

            # dip is only from 0-90 degrees, therefore uses only 2 digits

            if mt.cmtVersion is not None:
                ostream.write( '%-3s ' % mt.cmtVersion )
            else:
                ostream.write( '    ' )

            ( tev, tpl, taz ) = ( None, None, None )
            ( pev, ppl, paz ) = ( None, None, None )
            ( nev, npl, naz ) = ( None, None, None )

            if fm.principalAxes is not None:
                if fm.principalAxes.tAxis is not None:
                    tev     = fm.principalAxes.tAxis.length.value / power
                    tpl     = fm.principalAxes.tAxis.plunge.value
                    taz     = fm.principalAxes.tAxis.azimuth.value
                if fm.principalAxes.pAxis is not None:
                    pev     = fm.principalAxes.pAxis.length.value / power
                    ppl     = fm.principalAxes.pAxis.plunge.value
                    paz     = fm.principalAxes.pAxis.azimuth.value
                if fm.principalAxes.nAxis is not None:
                    nev     = fm.principalAxes.nAxis.length.value / power
                    npl     = fm.principalAxes.nAxis.plunge.value
                    naz     = fm.principalAxes.nAxis.azimuth.value


            ostream.write( '%7.3f' % ( tev is not None and [float(tev)] or [0.0] )[0] )
            ostream.write( '%3d' % ( tpl is not None and [int(tpl)] or [0] )[0] )
            ostream.write( '%4d' % ( taz is not None and [int(taz)] or [0] )[0] )

            ostream.write( '%8.3f' % ( pev is not None and [float(pev)] or [0.0] )[0] )
            ostream.write( '%3d' % ( ppl is not None and [int(ppl)] or [0] )[0] )
            ostream.write( '%4d' % ( paz is not None and [int(paz)] or [0] )[0] )

            ostream.write( '%8.3f' % ( nev is not None and [float(nev)] or [0.0] )[0] )
            ostream.write( '%3d' % ( npl is not None and [int(npl)] or [0] )[0] )
            ostream.write( '%4d' % ( naz is not None and [int(naz)] or [0] )[0] )

            if mt.scalarMoment.value is not None:
                ostream.write( '%8.3f' % scm_mantissa )
            else:
                ostream.write( '%8.3f' % ( 0.0 ) )

            ( s1, d1, r1 ) = ( None, None, None )
            ( s2, d2, r2 ) = ( None, None, None )
            
            if fm.nodalPlanes is not None:
                if fm.nodalPlanes.nodalPlane1 is not None:
                    s1     = fm.nodalPlanes.nodalPlane1.strike.value
                    d1     = fm.nodalPlanes.nodalPlane1.dip.value
                    r1     = fm.nodalPlanes.nodalPlane1.rake.value
                if fm.nodalPlanes.nodalPlane2 is not None:
                    s2     = fm.nodalPlanes.nodalPlane2.strike.value
                    d2     = fm.nodalPlanes.nodalPlane2.dip.value
                    r2     = fm.nodalPlanes.nodalPlane2.rake.value

            ostream.write( '%4d' % ( s1 is not None and [int(s1)] or [0] )[0] )
            ostream.write( '%3d' % ( d1 is not None and [int(d1)] or [0] )[0] )
            ostream.write( '%5d' % ( r1 is not None and [int(r1)] or [0] )[0] )

            ostream.write( '%4d' % ( s2 is not None and [int(s2)] or [0] )[0] )
            ostream.write( '%3d' % ( d2 is not None and [int(d2)] or [0] )[0] )
            ostream.write( '%5d' % ( r2 is not None and [int(r2)] or [0] )[0] )
            
            ostream.write( '\n' )

    # -------------------------------------------------------------------------
    
    def importANSSUnified( self, input, **kwargs ):
        """
        import data from ANSS "reduced" unified catalog, one event per line
        
        get monthly chunks:
            ftp://www.ncedc.org/pub/catalogs/cnss/YYYY/YYYY.MM.cnss

        data example for first 10 events from 2001.01.cnss (note: ruler shown below is not part of data)
        
        123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012

        $loc 20010101000748.8000 34.28100-118.45000 17.7100H CI   22              0.1930        0.0000 4.6660L 20010101     9172296 $magP 1.20h CI              20010101     9172296 $add$loc                                                                                  9172296     9172296
        $loc 20010101001207.7600 38.81816-122.78683  1.7500H NC   10 50    1.0000 0.0500        0.2800 0.3700L 20070622             $magP 1.36d NC    9 0.17    20070622             $add$loc  10   0  1016418    0.2300 8015    0.290020865    0.4100                        21141544            
        $locP20010101003246.2300 45.61300  26.47900125.0000H NEI  17              0.6900                               200101014001                                                                                                                                                               
        $locP20010101005217.3600 55.57900-159.70401 56.2000H AEI  17              0.0000                               200101014002 $magP 3.30l AEI                     200101014002                                                                                                              
        $locP20010101011510.7200 18.18200 -67.13300 23.2000H RSP   4              0.0000                               200101014003 $magP 2.60l RSP                     200101014003                                                                                                              
        $loc 20010101011625.7640 35.12300-118.53800  5.1600H CI    8              0.0350        0.0000 0.4010L 20010101     9172298 $magP 1.00h CI              20010101     9172298 $add$loc                                                                                  9172298     9172298
        $loc 20010101012155.9550 35.04700-119.08300 10.9400H CI   18              0.1770        0.0000 1.7770L 20010101     9172300 $magP 1.60h CI              20010101     9172300 $add$loc                                                                                  9172300     9172300
        $loc 20010101013345.5100 36.23650-120.79283  8.0500H NC   10157    2.0000 0.0500        0.6200 0.3800L 20070622             $magP 1.37d NC    4 0.13    20070622             $add$loc  10   1   9316 6    0.290023862    0.4300 4226    0.6900                        21141564            
        $loc 20010101013549.5960 34.47000-116.04500  4.0700H CI   14              0.1350        0.0000 0.8970L 20010101     9172301 $magP 1.63c CI              20010101     9172301 $add$loc                                                                                  9172301     9172301
        $loc 20010101013808.8400 37.62933-118.98534  8.0100H NC   32 83    2.0000 0.0400        0.1900 0.3400L 20070622             $magP 1.54d NC   21 0.16    20070622             $add$loc  35   2  24265 0    0.140017512    0.190035277    0.3500                        21141567            

        correct illegal time components: NO
        
        kwargs:
            authorityID - authority id used in public ids, default: 'ANSS'
        """

        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'authorityID' in kwargs and isinstance( kwargs['authorityID'], basestring ):
            auth_id = kwargs['authorityID']
        else:
            auth_id = 'ANSS'
            
        line_ctr = 0
        for line in istream:

            line_ctr += 1
            
            # in order to accept a location, we require focal time, lat, lon, and depth to be present
            # ignore lines shorter than 51 characters
            if len( line.strip() ) < 51:
                continue

            ## get location information ($loc)

            #columns  frmt   description
            #------- -----   -------------
            #1-  4   4s    $loc
            #5-  5    s    identification for preferred location ("P") when 
                        #multiple entries present, otherwise blank
            #6-  9  *4d    year of origin time (all four digits required)
            #10- 11  *2d    month of origin time (1-12)
            #12- 13  *2d    day of origin time (1-31)
            #14- 15  *2d    hour of origin time (0-23)
            #16- 17  *2d    minutes of origin time (0-59)
            #18- 24  *7.4f  seconds of origin time (0-59.9999)
            #25- 33  *9.5f  latitude in decimal degrees (-90.00000 -  90.00000,  N = +)
            #34- 43  *10.5f longitude in decimal degrees (-180.00000 -  180.00000, E = +)
            #44- 51  *8.4f  depth in km (datum reference defined by method, + down)
            #52- 53   2s    type of location (Table 1)
            #54- 56  *3s    source code of location information (Table 2a)
            #57- 60  *4d    number of non-null weighted travel times (P & S) 
                        #used to compute this hypocenter
            #61- 63   3d    azimuthal gap in degrees
            #64- 73   10.4f distance to nearest station in km
            #74- 80   7.4f  RMS residual of phases used in location
            #81- 87   7.4f  origin time error (s)
            #88- 94   7.4f  horizontal error (km)
            #95-101   7.4f  depth error (km)
            #102-103  *2s    auxillary event remarks (Table 3)
            #104-111    8d   date solution created in the form YYYYMMDD (eg, 19960125)
            #112-123  *12d   data center id # (right justified)
            #------- -----   -------------
            
            try:
                # first get required fields
                # NOTE: "data center id" is claimed to be a required field,
                #       but is often only whitespace in the cnss catalog files
                #       therefore it is not treated as 'required' here
                #       all fields right from depth can be missing in cnss files, are not truly required 
                #       we do not process lines that are shorter that 51 characters
                curr_year   = int( line[5:9] )
                curr_month  = int( line[9:11] )
                curr_day    = int( line[11:13] )
                curr_hour   = int( line[13:15] )
                curr_minute = int( line[15:17] )
                curr_second = float( line[17:24] )
                
                curr_lat    = float( line[24:33].strip() )
                curr_lon    = float( line[33:43].strip() )
                curr_depth  = float( line[43:51].strip() )

            except:
                print " importANSSUnified - error in $loc block of line %s: %s" % ( line_ctr, line )
                continue

            ## define event id
            # - do not use 'data center id' since is is not always set
            # - use location datetime block instead
            # - use datetime block also for magnitude URI below
            #   (format provides only one location & magnitude,
            #    so we there is no need to discriminate several locations and magnitudes for same event)
            curr_id = line[5:24].strip()
            
            # create event
            ev = Event( ''.join( ( 'smi:', auth_id, '/event/', curr_id ) ) )
            ev.add( self.eventParameters, 'event' )
            # self.eventParameters.event.append( ev )

            # create origin
            ori = Origin( ''.join( ( 'smi:', auth_id, '/origin/', curr_id ) ) )
            ori.add( ev, 'origin' )
            
            ori.latitude  = RealQuantity( curr_lat )
            ori.longitude = RealQuantity( curr_lon )
            ori.depth     = RealQuantity( curr_depth )
            ori.time      = TimeQuantity( ( curr_year, curr_month, curr_day,
                                            curr_hour, curr_minute, curr_second ) )
            
            # set preferred origin
            ev.preferredOriginID = ori.publicID

            ## get optional fields
            
            try:
                curr_loctype = line[51:53].strip()
                if curr_loctype in ( 'H', 'h' ):
                    ori.type = 'hypocenter'
                elif curr_loctype in ( 'C', 'c' ):
                    ori.type = 'centroid'
                elif curr_loctype in ( 'A', 'a' ):
                    ori.type = 'amplitude'    
            except: 
                  pass
                
            # this can be a network code or a different agency/institution ID
            # map this to creationInfo.agencyID
            try:
                curr_locsource = line[53:56].strip()
                
                if len( curr_locsource ) > 0:
                    if ( not hasattr( ori, 'creationInfo' ) ) or ( not isinstance( ori.creationInfo, CreationInfo ) ):
                        ori.creationInfo = CreationInfo()
                    ori.creationInfo.agencyID = curr_locsource
            except:
                pass
                
            try:    
                curr_traveltime_cnt = int( line[56:60].strip() )
                if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                    ori.quality   = OriginQuality()
                ori.quality.usedPhaseCount = curr_traveltime_cnt
            except: 
                pass
            
            try:
                curr_azimuthalGap = float( line[60:63].strip() )
                if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                    ori.quality   = OriginQuality()
                ori.quality.azimuthalGap = curr_azimuthalGap
            except: 
                pass

            # NOTE: in ANSS format distance to nearest station is given in km
            # this is good for QuakeML version <= 1.0.1
            # for QuakeML version > 1.0.1 degrees will be used
            try:
                curr_nearestStation = float( line[63:73].strip() )
                if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                    ori.quality   = OriginQuality()
                ori.quality.minimumDistance = curr_nearestStation
            except: 
                pass
            
            try:
                curr_rms_residual_phases = float( line[73:80].strip() )
                if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                    ori.quality   = OriginQuality()
                ori.quality.standardError = curr_rms_residual_phases
            except: 
                pass
            
            try:    ori.time.uncertainty = float( line[80:87].strip() )
            except: pass
            
            try:
                che = float( line[87:94].strip() )
                ou = OriginUncertainty()
                ou.horizontalUncertainty = che
                ou.add( ori, 'originUncertainty' )
            except:
                pass
            
            try:    ori.depth.uncertainty = float( line[94:101].strip() )
            except: pass
            
            # QuakeML does not have EventType entries for all possible values of ANSS format
            try:
                curr_ev_remarks = line[101:103].strip()
                if curr_ev_remarks in ( 'N', 'n' ):
                    ev.type = 'nuclear explosion'
                elif curr_ev_remarks in ( 'Q', 'q' ):
                    ev.type = 'quarry blast'
            except:
                pass
            
            # date of creation of solution, has format YYYYMMDD
            try:
                curr_solution_date = line[103:111].strip()
                cict = QPDateTime( ( int( curr_solution_date[0:4] ), 
                                     int( curr_solution_date[4:6] ), 
                                     int( curr_solution_date[6:8] ) ) )
                                 
                if ( not hasattr( ori, 'creationInfo' ) ) or ( not isinstance( ori.creationInfo, CreationInfo ) ):
                    ori.creationInfo = CreationInfo()
                ori.creationInfo.creationTime = cict
            except:
                pass
            
            # data center id: there is no good way to preserve 'legacy ids' in QuakeML
            # ignore this field
            # curr_loc_datacenterid = line[111:123].strip()
            
            ## get magnitude information ($mag)

            #columns  frmt   description
            #------- -----   -------------
            #1-  4   4s    $mag
            #5-  5    s    identification for preferred magnitude ("P") when 
                        #multiple entries present, otherwise blank
            #6- 10  *5.2f  magnitude
            #11- 12  *2s    magnitude type (Table 4)
            #13- 15  *3s    source code of magnitude information (Table 2a)
            #16- 19  *4d    number of observations for magnitude determination
            #20- 24   5.2f  error in magnitude estimate 
                        #(type of error depends on magnitude definition)
            #25- 28   4.1f  total of magnitude weights.
            #29- 36    8d   date solution created in the form YYYYMMDD (eg, 19960125)
            #37- 48  *12d   data center id # (right justified)
            #------- -----   -------
            
            # if line length > 123 chars, look for magitude information
            if len( line.strip() ) > 123:

                # copy rest of line into new zero-offset string
                # $mag block goes potentially from column 125 to column 172 (48 columns, list indices 124-171)
                # but can contain fewer characters, therefore we cannot specify an end position explicitly
                mag_info = line[124:]

                ## NOTE: although documentation says that some fields are required,
                ## there lines in the cnss files with almost ALL entries missing (even magnitude value)
                ## 'source code' of magnitude information seems to be always there 
                ## input line can end after 'source code'
                ## -> we treat none of the fields as required
                ## NOTE: last field "data center id" is not padded with whitespace if missing
                
                # we need a magnitude value (is a required attribute for QuakeML)
                # skip magnitude block if magnitude is missing
                validMagBlock = True
                try:
                    curr_mag = float( mag_info[5:10].strip() )
                except:
                    validMagBlock = False
                    
                if validMagBlock is True:
                  
                    mag = Magnitude( ''.join( ( 'smi:', auth_id, '/magnitude/', curr_id ) ) )
                    mag.add( ev, 'magnitude' )
                    mag.mag = RealQuantity( curr_mag )
                    mag.setOriginAssociation( ori.publicID )
                    
                    # set preferred magnitude
                    ev.preferredMagnitudeID = mag.publicID
                    
                    # write ANSS magnitude code in comment
                    try:
                        curr_magtype = mag_info[10:12].strip()
                        mc = Comment( ''.join( ( 'ANSS:magnitude_type=', curr_magtype ) ) )
                        mag.comment.append( mc )
                        
                        if curr_magtype in ( 'b', 'B' ):
                            mag.type = 'mb'
                        elif curr_magtype in ( 'l', 'L', 'l1', 'L1', 'l2', 'L2' ):
                            mag.type = 'ML'
                        elif curr_magtype in ( 's', 'S' ):
                            mag.type = 'Ms'
                        elif curr_magtype in ( 'w', 'W' ):
                            mag.type = 'Mw'
                        elif curr_magtype in ( 'c', 'C' ):
                            mag.type = 'Mc'
                        elif curr_magtype in ( 'd', 'D' ):
                            mag.type = 'Md'    
                    except:
                        pass
                    
                    try:
                        curr_magsource = mag_info[12:15].strip()
                        if len( curr_magsource ) > 0:
                            if ( not hasattr( mag, 'creationInfo' ) ) or ( not isinstance( mag.creationInfo, CreationInfo ) ):
                                mag.creationInfo = CreationInfo()
                            mag.creationInfo.agencyID = curr_magsource
                    except:
                        pass
                    
                    try:    mag.stationCount = int( mag_info[15:19].strip() )
                    except: pass
                    
                    try:    mag.mag.uncertainty = float( mag_info[19:24].strip() )
                    except: pass
                    
                    # total of magnitude weights not provided in QuakeML -> ignore 
                    # curr_mag_weights = float( line[24:28].strip() )
                    
                    # date of creation of solution, has format YYYYMMDD
                    try:
                        curr_mag_solution_date = mag_info[28:36].strip()
                        cict = QPDateTime( ( int( curr_mag_solution_date[0:4] ), 
                                             int( curr_mag_solution_date[4:6] ), 
                                             int( curr_mag_solution_date[6:8] ) ) )
                                            
                        if ( not hasattr( mag, 'creationInfo' ) ) or ( not isinstance( mag.creationInfo, CreationInfo ) ):
                            mag.creationInfo = CreationInfo()
                        mag.creationInfo.creationTime = cict
                    except:
                        pass
                
                    # data center id: there is no good way to preserve 'legacy ids' in QuakeML
                    # ignore this field
                    #curr_mag_datacenterid = line[36:48].strip()

            ## get additional location information ($addloc)

            #columns  frmt   description
            #------- -----   -------------
            #1-  4   4s    $add
            #5-  8   4s    $loc
            #9- 12   4d    number of valid P & S readings with non-null weights
            #13- 16   4d    number of S readings with non-null weights
            #17- 20   4d    number of P first motions
            #21- 23   3d    azimuth of smallest principal error (deg E of N)
            #24- 25   2d    dip of smallest principal error (deg)
            #26- 35   10.4f magnitude of smallest principal error (km)
            #36- 38   3d    azimuth of intermediate principal error (deg E of N)
            #39- 40   2d    dip of intermediate principal error (deg)
            #41- 50   10.4f magnitude of intermediate principal error (km)
            #51- 53   3d    azimuth of largest principal error (deg E of N)
            #54- 55   2d    dip of largest principal error (deg)
            #56- 65   10.4f magnitude of largest principal error (km)
            #66- 75   10.4f error in latitude    (km)
            #76- 85   10.4f error in longitude   (km)
            #86- 97   12d   local event id # (right justified)
            #98-109  *12d   data center id # (right justified)
            #------- -----   -------

            # if line length > 172 chars, look for additional location information
            if len( line.strip() ) > 172:

                # copy remaining $addloc block into new zero-offset string
                # $addloc block can go from column 174 to column 282 (109 columns, list indices 173-281)
                # NOTE: $addloc block can contain fewer characters, has no required fields
                addloc_info = line[173:]
                
                try:
                    if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                        ori.quality = OriginQuality()
                    ori.quality.associatedPhaseCount = int( addloc_info[8:12].strip() )
                except: 
                    pass
                    
                try:
                    first_motion_cnt = int( addloc_info[16:20].strip() )
                    fm = FocalMechanism( ''.join( ( 'smi:', auth_id, '/focalmechanism/', curr_id ) ) )
                    fm.stationPolarityCount = first_motion_cnt
                    fm.add( ev, 'focalMechanism' )
                    ev.preferredFocalMechanismID = fm.publicID
                except:
                    pass
                
                ## if lat/lon uncertainties in km are given, convert to degrees and add to quantity
                try:
                    ori.latitude.uncertainty = float( addloc_info[65:75].strip() ) / 111.0
                except:
                    pass
                
                try:
                    # if latitude is +/- 90 degrees, set longitude error to 0.0
                    # avoid division by zero
                    if abs( ori.latitude.value ) == 90.0:
                        ori.longitude.uncertainty = 0.0
                    else:
                        ori.longitude.uncertainty = float( addloc_info[75:85].strip() ) / ( 111.0 * math.cos( ori.latitude.value ) )
                except:
                    pass

    # -------------------------------------------------------------------------
    
    def importPDECompressed( self, input, **kwargs ):
        """
        import data from USGS/NEIC (PDE) catalog in "compressed" format, one event per line
        
        get data from web form:
            http://neic.usgs.gov/neis/epic/epic_global.html

        data example from web site (note: ruler shown below is not part of data)

           01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345

            GS    1987 0130222942.09  -60.063 -26.916 47  D171.416.2237.0Z266.90MsBRK  6.80MsPAS  153172...FG...........       
            GS    1987 0208183358.39   -6.088 147.689 54     1.31     7.4Z217.60MsBRK  7.00MsPAS  2072907C.FG..........M       
            GS    1987 0305091705.28  -24.388 -70.161 62  G  1.216.5667.3Z167.20MsBRK  7.00MsPAS  1223036C.FG.....T.....       
            GS    1987 0306041041.96    0.151 -77.821 10  G  1.226.5516.9Z217.00MsBRK  6.70MsPAS  106344.C.FG..........S       
            GS    1987 0706024942.78  -14.074 167.828 47     0.985.9526.6Z237.10MsBRK  6.50MsPAS  186304.F..G.....T.....       
            GS    1987 0903064013.91  -58.893 158.513 33  N  1.035.9437.3Z207.70MsBRK  6.90MsPAS  167292....G...........       
            GS    1987 1006041906.08  -17.940-172.225 16  G  1.036.7367.3Z267.30MsBRK  7.20MsPAS  174455.F.FG.....T.....       
            GS    1987 1016204801.64   -6.266 149.060 47  G  1.245.9267.4Z227.70MsBRK  7.00MsPAS  1923458D.FG.....T.....       
            GS    1987 1025165405.69   -2.323 138.364 33  N  1.126.2437.0Z196.70MsBRK  6.70MsPAS  201290.F.F............       
            GS    1987 1117084653.32AS 58.586-143.270 10  G      6.6766.9Z197.00MsBRK  7.00MLPMR  0155815FUFG.....T.....       
            GS    1987 1130192319.59   58.679-142.786 10  G      6.7577.6Z117.70MsBRK  7.10MLPMR  0155846DUFG.X...T.....       

            NOTE: column numbers are zero-offset (C style)

            0      blank
            1-5    Catalog Source (a5)
            6-11   Year (a6)
            12-13   Month (i2)
            14-15   Day (i2)
            16-24   Origin Time (f9.2)
            25-26   Coordinate/OT Auth.(a2)
            27-33   Latitude (f7.3) [-=South] 
            34-41   Longitude (f8.3) [-=West]
            42-44   Depth (i3)
            47      Depth Control Designator (a1)
            48-49   pP Phases( i2)
            50-53   Std. Dev.(4.2)
            54-56   mb magnitude (f3.1)
            57-58   mb obs (i2)
            59-61   Ms magnitude (f3.1)
            62      Z/H Component (a1)
            63-64   Ms obs. (i2)
            65-68   Magnitude1 (f4.2)
            69-70   Mag1. Scale (a2)
            71-75   Mag1 Donor (a5)
            76-79   Magnitude2 (f4.2)
            80-81   Mag2. Scale (a2)
            82-86   Mag2 Donor (a5)
            87-89   Region Number (i3)
            90-92   Sta. No./Qual. (a3)
            93      Io value (a1)
            94      Cultural Effect (a1)
            95      Isoseismal Map (a1)
            96      Fault Plane Sol. (a1)
            97      Moment Flag (a1)
            98      ISC Depth Flag (a1)
            99      IDE Flag (a1)
            100     Preferred Flag (a1)
            101     Flag (a1)
            102-108 Phenomena Codes (7a1)
            109-115 Radial Distance (a7)

        correct illegal time components: NO
        
        kwargs:
            authorityID - authority id used in public ids, default: 'PDE'
        """

        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'authorityID' in kwargs and isinstance( kwargs['authorityID'], basestring ):
            auth_id = kwargs['authorityID']
        else:
            auth_id = 'PDE'
            
        line_ctr = 0
        for line in istream:

            line_ctr += 1
            
            # in order to accept a location, we require focal time, lat, lon to be present
            # depth can be missing
            # ignore lines shorter than 42 characters
            # use rstrip() in order to keep leading space char
            if len( line.rstrip() ) < 42:
                continue

            try:
                # first get required fields
                # NOTE: minutes, seconds, and depth can be missing!
                curr_source = line[1:6].strip()
                
                curr_year   = int( line[6:12].strip() )
                curr_month  = int( line[12:14] )
                curr_day    = int( line[14:16] )
                
                curr_hour   = int( line[16:18] )
                curr_lat    = float( line[27:34].strip() )
                curr_lon    = float( line[34:42].strip() )
            except:
                print " importPDECompressed - error in line %s: %s" % ( line_ctr, line )
                continue

            # if minutes part is not given, set to 0
            try:
                curr_minute = int( line[18:20].strip() )
            except:
                curr_minute = 0
                
            # if seconds part is not given, set to 0.0
            try:
                curr_second = float( line[20:25].strip() )
            except:
                curr_second = 0.0

            ## define event id
            # - use year + datetime block 
            curr_id = ''.join( ( str(curr_year), line[12:25].strip() ) )
            
            # create event
            ev = Event( ''.join( ( 'smi:', auth_id, '/event/', curr_id ) ) )
            ev.add( self.eventParameters, 'event' )
            # self.eventParameters.event.append( ev )

            # create origin
            ori = Origin( ''.join( ( 'smi:', auth_id, '/origin/', curr_id ) ) )
            ori.add( ev, 'origin' )
            
            ori.latitude  = RealQuantity( curr_lat )
            ori.longitude = RealQuantity( curr_lon )
            ori.time      = TimeQuantity( ( curr_year, curr_month, curr_day,
                                            curr_hour, curr_minute, curr_second ) )
            
            # set preferred origin
            ev.preferredOriginID = ori.publicID

            ## get optional fields

            # depth
            try:
                curr_depth = float( line[42:45].strip() )
                ori.depth  = RealQuantity( curr_depth )
            except:
                pass
                    
            # location authority
            # map this to creationInfo.agencyID
            try:
                curr_locsource = line[25:27].strip()
                
                if len( curr_locsource ) > 0:
                    if ( not hasattr( ori, 'creationInfo' ) ) or ( not isinstance( ori.creationInfo, CreationInfo ) ):
                        ori.creationInfo = CreationInfo()
                    ori.creationInfo.agencyID = curr_locsource
            except:
                pass

            # depth control designator
            try:
                curr_depthcontrol = line[47:48].strip()

                if len( curr_depthcontrol ) > 0:
                    oc = Comment( ''.join( ( 'PDE:depth_control_designator=', curr_depthcontrol ) ) )
                    ori.comment.append( oc )
                        
                if curr_depthcontrol in ( 'A', 'a' ):
                    ori.depthType = 'operator assigned'
                elif curr_depthcontrol in ( 'D', 'd' ):
                    ori.depthType = 'constrained by depth phases'
                elif curr_depthcontrol in ( 'N', 'n' ):
                    ori.depthType = 'other'
                elif curr_depthcontrol in ( 'G', 'g' ):
                    ori.depthType = 'other'
                elif curr_depthcontrol in ( 'S', 's' ):
                    ori.depthType = 'other'
            except: 
                  pass

            # number of pP phases for depth determination
            try:    
                curr_pP_cnt = int( line[48:50].strip() )
                if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                    ori.quality   = OriginQuality()
                ori.quality.depthPhaseCount = curr_pP_cnt
            except: 
                pass
                
            # standard deviation of arrival time residuals
            try:    
                curr_std_dev = float( line[50:54].strip() )
                if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                    ori.quality   = OriginQuality()
                ori.quality.standardError = curr_std_dev
            except: 
                pass

            # Flinn-Engdahl region
            try:
                curr_fe_region = line[87:90].strip()
                if len( curr_fe_region ) > 0:
                    desc = EventDescription( curr_fe_region )
                    desc.type = 'Flinn-Engdahl region'
                    
                    ev.description.append( desc )
            except: 
                  pass
                
            # up to 4 magnitudes
            mag_indices = [ { 'from': 54, 'to': 57, 'obs_from': 57, 'obs_to': 59, 'mag_type': 'mb' },
                            { 'from': 59, 'to': 62, 'obs_from': 63, 'obs_to': 65, 'mag_type': 'Ms' },
                            { 'from': 65, 'to': 69, 'magtype_from': 69, 'magtype_to': 71, 'magsrc_from': 71, 'magsrc_to': 76 },
                            { 'from': 76, 'to': 80, 'magtype_from': 80, 'magtype_to': 82, 'magsrc_from': 82, 'magsrc_to': 87 } ]
                            
            for mag_ctr in xrange( 0, 4 ):
                
                validMagBlock = True
                try:
                    curr_mag = float( line[mag_indices[mag_ctr]['from']:mag_indices[mag_ctr]['to']].strip() )
                except:
                    validMagBlock = False

                if validMagBlock is True:

                    mag = Magnitude( ''.join( ( 'smi:', auth_id, '/magnitude/', curr_id, '/', str(mag_ctr+1) ) ) )
                    mag.add( ev, 'magnitude' )
                    mag.mag = RealQuantity( curr_mag )
                    mag.setOriginAssociation( ori.publicID )

                    # magnitude type
                    try:
                        if mag_ctr in ( 0, 1 ):
                            curr_magtype = mag_indices[mag_ctr]['mag_type']
                        else:    
                            curr_magtype = line[mag_indices[mag_ctr]['magtype_from']:mag_indices[mag_ctr]['magtype_to']].strip()

                        if len( curr_magtype ) > 0:
                            mag.type = curr_magtype
                    except:
                        pass

                    # source of magnitude determination (for mb and Ms it is 'NEIS')
                    try:
                        if mag_ctr in ( 0, 1 ):
                            curr_magsource = 'NEIS'
                        else:
                            curr_magsource = line[mag_indices[mag_ctr]['magsrc_from']:mag_indices[mag_ctr]['magsrc_to']].strip()
                            
                        if len( curr_magsource ) > 0:
                            if ( not hasattr( mag, 'creationInfo' ) ) or ( not isinstance( mag.creationInfo, CreationInfo ) ):
                                mag.creationInfo = CreationInfo()
                            mag.creationInfo.agencyID = curr_magsource
                    except:
                        pass
                                
                    # number of observations (only for mb and Ms)
                    if mag_ctr in ( 0, 1 ):
                        try:
                            curr_obs = int( line[mag_indices[mag_ctr]['obs_from']:mag_indices[mag_ctr]['obs_to']].strip() )
                            mag.stationCount = curr_obs
                        except:
                            pass

            ## set preferred magnitude
            # - take first given magnitude from the 4 positions as preferred
            # - usually this is mb
            # - if neither mb nor Ms is given, use first of additional magnitudes
            # - this is usually ML or Mw
             
            if len( ev.magnitude ) > 0:
                ev.preferredMagnitudeID = ev.magnitude[0].publicID

    # -------------------------------------------------------------------------
    
    def importJMADeck( self, input, **kwargs ):
        """
        import data from Japanese JMA catalog in "deck" format
        
        data example (note: ruler shown below is not part of data)

            123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456

            J2007040100030476 011 370929 018 1363520 053  76511322V   511   4167OFF NOTO PENINSULA       29K
            _N.TGIH2511 h 1P   00030729S   030908 6899 01  2 5235 01  2 3867 00  01                 7 4MM  1
            _N.SHKH2507 h 1P   00030929S   031241 2047 01  3 3035 01  3  953 01  31                 7 4MM  1
            _HAKUI  575 N 1P   00031046S   031464                       1985 01  51                 7 4MM  1
            _E.WAJ  811 % 1P   00031064S   031496                        405 01  54                 7 4MM  1
            _N.AMZH2508 h 1P   00031106S   031527 4930 01  5 4550 01  4 3348 00  4J                 7 4MM  1
            E
            J2007040110182065 033 304356 099 1423717 561 47     41V   571   8329FAR E OFF IZU ISLANDS    24S
            U2007040110181893     303624     1422862     10     36B         8   SOUTH OF HONSHU, JAPAN      
            _CHIJI3 589 N 1EP  10191592                                 1335 03 47J                 7 4N   1
            _KUWANO 692 N 1EP  10191638                                                             7 4N   1
            _BS1OBS 159 Q 1EP  10192333ES  201007                                                   7 4NN  1
            _BS2OBS 160 Q 1ES  10201622                                                             7 4N   1
            _BS3OBS 161 Q 1EP  10192722ES  201757                                                   7 4NN  1
            _N.ST2H 983 r 1EP  10192957ES  201998 4129 01 57 6330 01 57 1880 02 69J                 7 4NN  1
            _BS4OBS 162 Q 1EP  10193229ES  202262                        244 01  3J                 7 4NN  1
            _N.KIBH2770 s 1EP  10195294ES  210031                                                   7 4NN  1
            _MATSUS  67 G 1EP  10200041ES  211542                        106 03  2J                 7 4NN  1
            _MATSUS  67 H01P   1020010 S   21169                                                    7 4   R0
            E
            U2007040114175318    -175652    -1783258    608     49B         9   FIJI REGION
            J2007040114274820                                          8    8400FAR FIELD                  F
            _RYOKAM 561 N 1P   14274820                                                             7 4M   1
            _MATSUS  67 G 1P   14275230                                                             7 4M   1
            _WACHI  592 N 1P   14275740                                                             7 4M   1
            _AIDA   608 N 1P   14280122                                                             7 4M   1
            _SAIJYO 605 N 1P   14280515                                                             7 4M   1
            _KAMIAS 505 N 1P   14281401                                                             7 4M   1
            _SHIMAM 511 N 1P   14281389                                                             7 4M   1
            _MATSUS  67 H01P   1427522                                                              7 4   R0
            E

            123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456

            NOTE: this is an illegal line, 2nd magnitude does not have correct format
                  we interpret this as a magnitude 8.0
                                                                   ||
            U198703061839542     -24122     - 70044      45     57B8 S      9   NEAR COAST OF N-CHILE  5

            123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456

            NOTE: there can be a block like the following after the hypocenter line
                  we ignore all lines that are not hypocenter, phase, comment, or end
        
            J199201012203556  01  34084  06  135355  05  681 19 39D40V1112  5189NE WAKAYAMA PREF           K
            MI JMAM       1N AM  91 67 206 10 300 20 273 39 -123 133 58  -66 92  52  91  35        92 1LOW
            NUNZEND  93 S 1            # 93  32 4390  130 1620                             U       92 1JMA
            NSHIMO2 130 S 1            #130  34  750  131  650                             U       92 1JMA
            NKIRISH 155 S 1            #155  31 5380  130 5240                             U       92 1JMA
            NOOSAKI 364 S 1            #364  35   28  139  607                             U       92 1JMA
            NASAHI  408 S 1            #408  36  733  137 5117                             U       92 1JMA
            NHOKIGI 412 S 1            #412  34 5098  139  238                             U       92 1JMA
            NOKABE  431 S 1            #431  34 5700  138 1523                             U       92 1JMA
            NNAKAIZ 432 S 1            #432  34 5477  138 5980                             U       92 1JMA
            NTENRYU 433 S 1            #433  34 5447  137 5312                             U       92 1JMA
            NKOZU   440 S 1            #440  34 1177  139  835                             U       92 1JMA
            NINUYAM 451 S 1            #451  35 2098  137  172                             U       92 1JMA
            NMIKAWA 452 S 1            #452  34 4575  137 2820                             U       92 1JMA
            NTSUKEC 454 S 1            #454  35 3920  137 2797                             U       92 1JMA
            NOSHIKA 455 S 1            #455  35 3478  138  293                             U       92 1JMA
            NTOYONE 456 S 1            #456  35  813  137 4462                             D       92 1JMA

            123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
            
            NOTE: there are invalid event blocks like the following (focal time information is not complete)
                  such blocks are discarded

            J198302010039                                              81       INSUFFICIENT N OF DATA
            _MATSUS  670U 1IP  0039408 S   5000                                            U       83 2
            E

            NOTE: there can be phase lines like the following, which are discard
            _JNEMUR4320 j28M   2322               5829 30   18780 33    5094 42   2                 0 1

            
        all lines have 96 regular chars, with trailing hex '0A'
        separator lines start with 'E' (hex '45'), then 95x space char (hex '20')

        NOTE: there can be blank lines, e.g., line 158 in catalog file 198301.deck.Z
              we just skip them

        first there are hypocenter lines starting with:
            J: JMA location
            U: USGS location
            I: ISC location
        then there are optional comment lines starting with 'C'
        then there are phase lines starting with '_' which can contain one OR two phases
        NOTE: there can be other line types which will be ignored
        finally there is a separator line starting with 'E'

        phases belong to JMA location
        JMA location is not always preferred (sometimes it is incomplete)

        preferred origin:
         - JMA if complete
         - other (USGS) if complete
         - if no one is complete, use JMA

        correct illegal time components: YES (use QPUtils.fixTimeComponents)
        
        kwargs: nopicks        = True - do not read in pick lines
                jmaonly        = True - import only JMA origins
                minimumDataset = True - only read basic information (save memory)
                authorityID           - authority id used in public ids, default: 'JMA'
        """
        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        if 'minimumDataset' in kwargs and kwargs['minimumDataset'] is True:
            minConf = True
        else:
            minConf = False

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'authorityID' in kwargs and isinstance( kwargs['authorityID'], basestring ):
            auth_id = kwargs['authorityID']
        else:
            auth_id = 'JMA'
            
        # time shift in hours of Japan Standard Time
        time_delta_jst = 9.0

        mag_indices = [ { 'from': 52, 'to': 54, 'magtype_from': 54, 'magtype_to': 55 },
                        { 'from': 55, 'to': 57, 'magtype_from': 57, 'magtype_to': 58 } ]
                                
        line_ctr = 0

        phasesMode     = None
        hypocenterMode = None
        newEventMode   = None
        commentMode    = None

        comment_str_xml = ''
        
        for line in istream:

            line_ctr += 1

            # if line is blank, skip
            if len( line.strip() ) == 0:
                continue
            
            # check if it is hypocenter, comment, pick, or separator line
            line_type = line[0]

            if line_type in ( 'J', 'U', 'I' ):

                # previous mode can be: None, Hypocenter
                if hypocenterMode is None:
                    newEventMode    = True
                    validHypocenter = False
                    jma_origin      = None
                else:
                    newEventMode = False
                    
                hypocenterMode = True
                phasesMode     = False
                
            elif line_type == '_':

                # previous mode can be: Hypocenter, Phases, Comment
                hypocenterMode = False
                phasesMode     = True

                # if nopicks set or no valid hypocenter exists, skip line
                if (    ( 'nopicks' in kwargs and kwargs['nopicks'] is True )
                     or validHypocenter is False ):
                    continue

            elif line_type == 'C':

                # previous mode can be: Hypocenter, Comment

                hypocenterMode = False
                commentMode    = True
                
            elif line_type == 'E':
                phasesMode     = None
                hypocenterMode = None
                newEventMode   = None
                commentMode    = None

                # if no valid hypocenter has been found (no event object exists), go to next line
                if validHypocenter is True:
                
                    # add comment to event
                    if len( comment_str_xml.strip() ) > 0:

                        if minConf is not True:
                            ct = Comment( comment_str_xml.strip() )
                            ev.comment.append( ct )

                        comment_str_xml = ''

                    # assign preferred origin and preferred magnitude
                    if len( ev.origin ) == 1:
                        ev.preferredOriginID = ev.origin[0].publicID

                        # get first magnitude for that origin, if there are any
                        if len( ev.magnitude ) > 0:
                            ev.preferredMagnitudeID = ev.getMagnitudes( ev.origin[0] )[0].publicID
                        
                    else:
                        # there are more than one origins
                        # look for JMA origin and see if it is complete

                        # TODO: at first, use JMA origin even if it is incomplete
                        for my_ori in ev.origin:
                            if my_ori.creationInfo.agencyID == 'JMA':
                                ev.preferredOriginID = my_ori.publicID

                                # get first magnitude for that origin
                                ori_mag = ev.getMagnitudes( my_ori )
                                if len( ori_mag ) > 0:
                                    ev.preferredMagnitudeID = ori_mag[0].publicID

                # proceed to next input line
                continue
                
            else:
                # print " importJMADeck: ignoring unknown line type in line %s: %s" % ( line_ctr, line )
                continue

            if hypocenterMode is True:

                # on entering this block for the first hypocenter line, validHypocenter is False
                # set this to True if a good hypocenter (with full time information) has been found
                try:
                    # first get required fields
                    # NOTE: coordinates and depth need not be there
                    # NOTE: sometimes seconds are not there, this is not a valid hypocenter
                    curr_source = line[0]
                    
                    curr_year   = int( line[1:5] )
                    curr_month  = int( line[5:7] )
                    curr_day    = int( line[7:9] )
                    
                    curr_hour   = int( line[9:11] )
                    curr_minute = int( line[11:13] )

                    # seconds are given w/o decimal point
                    curr_second_int  = int( line[13:15] )
                    curr_second_frac = line[15:17].strip()

                except:
                    # print " importJMADeck: time information in hypocenter line not complete, line %s: %s" % ( line_ctr, line )
                    continue

                # if only JMA hypocenters are requested, discard all others
                if (     'jmaonly' in kwargs and kwargs['jmaonly'] is True
                     and curr_source != 'J' ):
                    # print " importJMADeck: discarding non-JMA hypocenter, line %s: %s" % ( line_ctr, line )
                    continue
                    
                ## define event id
                # - use agency + datetime block
                pick_ctr = 0
                curr_id  = line[0:17].strip()
                
                # create event
                ev = Event()
                ev.add( self.eventParameters, 'event' )
                
                if minConf is not True:
                    ev.publicID = ''.join( ( 'smi:', auth_id, '/event/', curr_id ) )

                # create origin
                ori = Origin()
                ori.add( ev, 'origin' )
                
                if minConf is not True:
                    ori.publicID = ''.join( ( 'smi:', auth_id, '/origin/', curr_id ) )
                    
                ## check if time components are well-behaved and fix if necessary
                # seconds and minutes are sometimes set to 60
                # possibly hours can be set to 24
                timeCorrection = fixTimeComponents( curr_hour, curr_minute, curr_second_int )

                try:
                    curr_second = float( '.'.join( ( str(int(timeCorrection['component'][2])), curr_second_frac ) ) )
                except:
                    print " importJMADeck: illegal time format (seconds) in hypocenter line %s: %s" % ( line_ctr, line )
                    continue

                # if we have arrived here, we call this a valid hypocenter (valid time information is there)
                validHypocenter = True
                    
                focal_time_utc = DateTime( curr_year,
                                           curr_month,
                                           curr_day,
                                           timeCorrection['component'][0],
                                           timeCorrection['component'][1],
                                           curr_second ) - TimeDelta( time_delta_jst )
                                           
                                           
                focal_time_utc = adjustDateTime( timeCorrection['increaseFlag'], focal_time_utc )
                ori.time = TimeQuantity( focal_time_utc )
                
                # set agency
                ori.creationInfo = CreationInfo()

                if curr_source == 'J':
                    ori.creationInfo.agencyID = 'JMA'
                    jma_origin                = ori

                    # TODO: set preferred origin
                    ev.preferredOriginID = ori.publicID
                
                elif curr_source == 'U':
                    ori.creationInfo.agencyID = 'USGS'
                elif curr_source == 'I':
                    ori.creationInfo.agencyID = 'ISC'

                ## get optional fields

                # latitude is given as degrees + decimal minutes (w/o decimal point)
                try:
                    curr_lat_deg      = line[22:24].strip()
                    curr_lat_min_int  = line[24:26].strip()
                    curr_lat_min_frac = line[26:28].strip()

                    curr_lat_min      = float( '.'.join( ( curr_lat_min_int, curr_lat_min_frac ) ) )
                    curr_lat          = float( curr_lat_deg ) + ( curr_lat_min / 60.0 )

                    ori.latitude = RealQuantity( curr_lat )

                except:
                    pass

                # longitude is given as degrees + decimal minutes (w/o decimal point)
                try:
                    curr_lon_deg      = line[33:36].strip()
                    curr_lon_min_int  = line[36:38].strip()
                    curr_lon_min_frac = line[38:40].strip()
                    curr_lon_min      = float( '.'.join( ( curr_lon_min_int, curr_lon_min_frac ) ) )
                    curr_lon          = float( curr_lon_deg ) + ( curr_lon_min / 60.0 )

                    ori.longitude = RealQuantity( curr_lon )

                except:
                    pass

                # depth
                try:
                    curr_depth_int  = line[44:47].strip()
                    curr_depth_frac = line[47:49].strip()

                    # if depth was determined using 'depth slice method', no fraction is there
                    if len( curr_depth_frac ) > 0:
                        curr_depth = float( '.'.join( ( curr_depth_int, curr_depth_frac ) ) )
                    else:
                        curr_depth = float( curr_depth_int )

                    ori.depth = RealQuantity( curr_depth )

                except:
                    pass

                ## do not read uncertainties if minimum configuration is selected
                if minConf is not True:

                    # focal time error (seconds)
                    try:
                        curr_second_err_int  = line[17:19].strip()
                        curr_second_err_frac = line[19:21].strip()
                        curr_second_err      = float( '.'.join( ( curr_second_err_int,
                                                                    curr_second_err_frac ) ) )

                        ori.time.uncertainty = curr_second_err
                    except:
                        pass
                    
                    # latitude minutes error
                    try:
                        curr_lat_err_min_int  = line[28:30].strip()
                        curr_lat_err_min_frac = line[30:32].strip()
                        curr_lat_err          = float( '.'.join( ( curr_lat_err_min_int,
                                                                   curr_lat_err_min_frac ) ) ) / 60.0
                        ori.latitude.uncertainty = curr_lat_err
                    except:
                        pass
                    
                    # longitude minutes error
                    try:
                        curr_lon_err_min_int  = line[40:42].strip()
                        curr_lon_err_min_frac = line[42:44].strip()
                        curr_lon_err          = float( '.'.join( ( curr_lon_err_min_int,
                                                                   curr_lon_err_min_frac ) ) ) / 60.0
                        ori.longitude.uncertainty = curr_lon_err
                    except:
                        pass
                    
                    # depth error
                    try:
                        # NOTE: format %3.2f, so 9.99 km is maximum error!
                        curr_depth_err_int  = line[49:50].strip()
                        curr_depth_err_frac = line[50:52].strip()
                        curr_depth_err = float( '.'.join( ( curr_depth_err_int, curr_depth_err_frac ) ) )

                        ori.depth.uncertainty = curr_depth_err
                    except:
                        pass
                  
                ## magnitudes
                for mag_ctr in xrange( 2 ):

                    curr_mag_code = line[mag_indices[mag_ctr]['from']:mag_indices[mag_ctr]['to']]

                    # something there?
                    if len( curr_mag_code.strip() ) > 0:

                        # magnitude entry should have 2 characters, but sometimes a charcater is missing
                        # last character missing (this occurs, e.g., in the 1987/03 USGS entry U198703061839542):
                        #  - assume missing character as '0'
                        # first character missing (don't know if there is a case):
                        # - illegal magnitude, don't register

                        # check if one component is whitespace
                        if len( curr_mag_code ) > len( curr_mag_code.strip() ):

                            if len( curr_mag_code[0].strip() ) == 0:

                                # not fixable, discard entry
                                print " importJMADeck: illegal one-character magnitude format %s in hypocenter line %s: %s" % (
                                                                                     curr_mag_code, line_ctr, line )
                                continue

                            elif len( curr_mag_code[1].strip() ) == 0:

                                # fix magnitude code
                                curr_mag_code = ''.join( ( curr_mag_code[0], '0' ) )
                              
                        try:
                            # is magnitude code numeric? (mag. >= -0.9)
                            curr_mag_code_int = int( curr_mag_code )

                            if curr_mag_code_int >= 0:

                                # mag >= 0 (F2.1 w/o decimal point)
                                try:
                                    curr_mag = float( '.'.join( ( curr_mag_code[0], curr_mag_code[1] ) ) )
                                except:
                                    print " importJMADeck: illegal positive numeric magnitude format %s in hypocenter line %s: %s" % ( curr_mag_code, line_ctr, line )
                                    continue
                            else:

                                # -1, -2, ..., -9 (-0.9 <= mag <= 0.1)
                                try:
                                    curr_mag = float( '.'.join( ( '-0', curr_mag_code[1] ) ) )
                                except:
                                    print " importJMADeck: illegal negative numeric magnitude format %s in hypocenter line %s: %s" % ( curr_mag_code, line_ctr, line )
                                    continue
                        except:

                            # mag. <= -1.0, code not numeric, use letter code
                            # first char has to be letter A, B, or C
                            # second char has to be integer number

                            if (     curr_mag_code[0] in ( 'A', 'B', 'C' )
                                 and curr_mag_code[1].isdigit() ):
                            
                                letter_code = { 'A': '-1', 'B': '-2', 'C': '-3' }
                                curr_mag = float( '.'.join( ( letter_code[curr_mag_code[0]], curr_mag_code[1] ) ) )

                            else:
                                print " importJMADeck: illegal magnitude format %s in hypocenter line %s: %s" % (
                                                                                        curr_mag_code, line_ctr, line )
                                continue

                        mag = Magnitude()
                        if minConf is not True:
                            mag.publicID = ''.join( ( 'smi:', auth_id, '/magnitude/', curr_id, '/', str(mag_ctr+1) ) )
                    
                        mag.add( ev, 'magnitude' )
                        mag.mag = RealQuantity( curr_mag )
                        mag.setOriginAssociation( ori.publicID )

                        # magnitude type
                        try:
                            curr_magtype = line[mag_indices[mag_ctr]['magtype_from']:mag_indices[mag_ctr]['magtype_to']].strip()
                        except:
                            curr_magtype = ''
                            
                        if curr_magtype == 'J':
                            mag.type = 'MJ'
                        elif curr_magtype in ( 'D', 'd' ):
                            mag.type = 'MD'
                        elif curr_magtype in ( 'V', 'v' ):
                            mag.type = 'MV'
                        elif curr_magtype == 'B':
                            mag.type = 'mb'
                        elif curr_magtype == 'S':
                            mag.type = 'MS'
                        else:
                            mag.type = 'unknown'

                        # source of magnitude determination
                        # mb: USGS, MS: ?, other: JMA
                        if mag.type not in ( 'MS', 'unknown' ):

                            if mag.type == 'mb':
                                curr_magsource = 'USGS'
                            else:
                                curr_magsource = 'JMA'

                            if minConf is not True:
                                mag.creationInfo = CreationInfo()
                                mag.creationInfo.agencyID = curr_magsource


                ## set first as preferred magnitude
                if len( ev.magnitude ) > 0:
                    ev.preferredMagnitudeID = ev.magnitude[0].publicID

                ## get "subsidiary information": event type
                try:
                    subsidiary = line[60].strip()

                    # valid classifications are
                    # '1': Natural earthquake
                    # '2': Insufficient number of JMA station
                    # '3': Artificial event
                    # '4': Noise
                    # '5': Low frequency earthquake

                    # we consider '1', '2', and '5' as earthquakes

                    if len( subsidiary ) > 0:
                        
                        if subsidiary in ( '1', '2', '5' ):
                            ev.type = 'earthquake'

                        # comment: 'JMA:subsidiary=<subsidiary>'
                        oc = Comment( ''.join( ( 'JMA:subsidiary=', subsidiary ) ) )
                        ev.comment.append( oc )
                        
                except:
                    pass
                        
                ## geographical region
                try:
                    curr_region_str = line[68:92].strip()

                    if len( curr_region_str ) > 0:
                        curr_region_str_xml = saxutils.escape( curr_region_str )

                        if minConf is not True:
                            descr = EventDescription( curr_region_str_xml, 'region name' )
                            ev.description.append( descr )
                except:
                    pass
                    
                ## origin quality
                try:
                    curr_stations_used = int( line[92:95].strip() )
                    if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                        ori.quality   = OriginQuality()
                    if minConf is not True:
                        ori.quality.usedStationCount = curr_stations_used
                except:
                    pass

            elif phasesMode is True:

                # check if there is a JMA origin for picks
                if jma_origin is None:
                    error_str = "importJMADeck: no JMA origin to associate phases to for event %s" % ev.publicID
                    raise ValueError, error_str

                # first phase type
                curr_phase_first = line[15:19].strip()

                # ignore line if empty
                if len( curr_phase_first ) == 0:
                    continue
                
                ## required fields
                try:
                    curr_sta_code             = line[1:7].strip()
                    curr_date_day             = int( line[13:15] )
                    curr_time_first_hours     = int( line[19:21] )
                    curr_time_first_mins      = int( line[21:23] )
                    curr_time_first_secs_int  = int( line[23:25] )
                    curr_time_first_secs_frac = line[25:27]
                    
                except:
                    print " importJMADeck: format error in phase line %s: %s" % ( line_ctr, line )
                    continue

                # create pick
                pick_ctr = pick_ctr + 1
                pick_id  = line[19:27].strip()

                pick = Pick()
                if minConf is not True:
                    pick.publicID = ''.join( ( 'smi:', auth_id, '/event/', curr_id, '/pick/', pick_id ) )
                            
                # set pick time, use year and month from event

                # check if date/day given in phase line is smaller than date/day of focal time:
                # - either pick is in following month (day of focal time is at end of month)
                #   -> add one month to date/month, difference between days >~ 25
                # - or pick time is earlier than focal time AND focal time is shortly after midnight
                #   -> no correction, difference between days = 1
                #
                # it can happen that pick time is earlier than focal time
                # strange, but we will ignore this

                timeCorrection = fixTimeComponents( curr_time_first_hours,
                                                    curr_time_first_mins,
                                                    curr_time_first_secs_int )

                try:
                    curr_time_first_secs = float( '.'.join( ( str(int(timeCorrection['component'][2])),
                                                              curr_time_first_secs_frac ) ) )
                except:
                    print " importJMADeck: illegal time format in phase line %s: %s" % ( line_ctr, line )
                    continue

                pick_time_utc = DateTime( curr_year,
                                          curr_month,
                                          curr_date_day,
                                          timeCorrection['component'][0],
                                          timeCorrection['component'][1],
                                          curr_time_first_secs ) - TimeDelta( time_delta_jst )
                                           
                pick_time_utc = adjustDateTime( timeCorrection['increaseFlag'], pick_time_utc )

                if ( curr_day - curr_date_day ) > 1:

                    # add one month
                    pick_time_utc += RelativeDateTime( months=+1 )
                        
                pick.time = TimeQuantity( pick_time_utc )

                # set waveform id
                pick.waveformID = WaveformStreamID()
                pick.waveformID.stationCode = curr_sta_code

                # preliminary: network code 'JMA'
                pick.waveformID.networkCode = 'JMA'
                
                pick.add( ev, 'pick' )
                
                # create arrival: phase 
                arrv = Arrival()
                arrv.pickID = pick.publicID
                arrv.phase  = Phase( curr_phase_first )
                arrv.add( jma_origin, 'arrival' )

                ## look is second phase is given
                
                # second phase type
                curr_phase_second = line[27:31].strip()

                # ignore line if empty
                if len( curr_phase_second ) == 0:
                    continue
                
                ## required fields
                try:
                    curr_time_second_mins      = int( line[31:33] )
                    curr_time_second_secs_int  = int( line[33:35] )
                    curr_time_second_secs_frac = line[35:37]
                except:
                    print " importJMADeck: time format (seconds) error in phase line %s: %s" % ( line_ctr, line )
                    continue

                # create pick
                pick_ctr = pick_ctr + 1
                pick_id  = ''.join( ( str(curr_time_first_hours), line[31:37].strip() ) )

                pick = Pick()
                if minConf is not True:
                    pick.publicID = ''.join( ( 'smi:', auth_id, '/event/', curr_id, '/pick/', pick_id ) )
                
                # set pick time, use date from event and hours from first pick of line

                # check if pick time is earlier than time of first pick if hours from first pick is used
                # -> then second pick is in following hour
                # first use date from first pick, if required, add one hour
                # get components from DateTime object of first pick
                # no correction for time zone required here!

                timeCorrection = fixTimeComponents( pick_time_utc.hour,
                                                    curr_time_second_mins,
                                                    curr_time_second_secs_int )

                try:
                    curr_time_second_secs = float( '.'.join( ( str(int(timeCorrection['component'][2])),
                                                               curr_time_second_secs_frac ) ) )
                except:
                    print " importJMADeck: illegal time format in phase line %s: %s" % ( line_ctr, line )
                    continue

                pick_time_second_utc = DateTime( pick_time_utc.year,
                                                 pick_time_utc.month,
                                                 pick_time_utc.day,
                                                 timeCorrection['component'][0],
                                                 timeCorrection['component'][1],
                                                 curr_time_second_secs )
                                                 
                pick_time_second_utc = adjustDateTime( timeCorrection['increaseFlag'], pick_time_second_utc )

                if pick_time_second_utc < pick_time_utc:

                    # add one hour
                    pick_time_second_utc += TimeDelta( 1.0 )
                        
                pick.time = TimeQuantity( pick_time_second_utc )

                # set waveform id
                pick.waveformID = WaveformStreamID()
                pick.waveformID.stationCode = curr_sta_code

                # preliminary: network code 'JMA'
                pick.waveformID.networkCode = 'JMA'
                
                pick.add( ev, 'pick' )
                
                # create arrival: phase 
                arrv = Arrival()
                arrv.pickID = pick.publicID
                arrv.phase  = Phase( curr_phase_second )
                arrv.add( jma_origin, 'arrival' )

            elif commentMode is True:

                # append to comment text string
                # comment object is created and added to Event when data block ends ('E')
                comment_str_xml = ' '.join( ( comment_str_xml, saxutils.escape( line[2:].strip() ) ) )

    # -------------------------------------------------------------------------
    
    def importGSE2_0Bulletin( self, input, **kwargs ):
        """
        import earthquake catalog data in GSE2.0 Bulletin format (as used, e.g., for the INGV catalog)
        
        data example from INGV (note: rulers shown at beginning and end are not part of data)

            123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012

            BEGIN GSE2.0
            MSG_TYPE DATA
            MSG_ID 2008-03-21_16:15:35 ITA_NDC
            E-MAIL info@example.com
            DATA_TYPE BULLETIN GSE2.0

            EVENT 00011371
            Date       Time            Lat       Lon    Depth    Ndef Nsta Gap    Mag1  N    Mag2  N    Mag3  N  Author          ID
                rms   O.T. Error    Smajor Sminor Az        Err  mdist   Mdist     Err        Err        Err     Quality

            2008/02/16 01:22:14.9      43.935    10.289      5.1      11    6 170  Md 2.2  6  Ml 1.7  4             ITA_NDC   00011371
                0.32    +-  0.18       2.3    1.3   91    +-  2.7    0.15   0.72   +-0.3      +-0.1      +-        m i ke

            ITALY (Alpi Apuane)
            Sta    Dist   EvAz     Phase         Date       Time  TRes  Azim  AzRes  Slow  SRes Def  SNR        Amp   Per   Mag1   Mag2 Arr ID
            MAIM    0.15    98 m   Pg      2008/02/16 01:22:18.8   0.4                          T                                       00149600
            MAIM    0.15    98 m   Sg      2008/02/16 01:22:21.1   0.1                          T               273   .21 Md 1.9 Ml 1.9 00149601
            VLC     0.23    17 m   Pg      2008/02/16 01:22:20.4   0.1                          T                                       00149602
            VLC     0.23    17 m   Sg      2008/02/16 01:22:24.1  -0.1                          T                96   .19 Md 2.1 Ml 1.6 00149603
            BDI     0.26    60 m   Pg      2008/02/16 01:22:20.8   0.1                          T                                       00149604
            BDI     0.26    60 m   Sg      2008/02/16 01:22:24.5  -0.5                          T               124   .70 Md 2.0 Ml 1.7 00149605
            PII     0.27   141 m   Pg      2008/02/16 01:22:21.1  -0.0                          T                                       00149606
            PII     0.27   141 m   Sg      2008/02/16 01:22:25.2  -0.5                          T                         Md 2.0        00149607
            ERBM    0.49    10 m   Pg      2008/02/16 01:22:25.8   0.2                          T                                       00149608
            ERBM    0.49    10 m   Sg      2008/02/16 01:22:33.4   0.0                          T               120   1.3 Md 2.5 Ml 1.7 00149609
            SC2M    0.72   311 m   Pg      2008/02/16 01:22:29.2  -0.2                          T                         Md 2.4        00149610

            .

            EVENT 00011372
            Date       Time            Lat       Lon    Depth    Ndef Nsta Gap    Mag1  N    Mag2  N    Mag3  N  Author          ID
                rms   O.T. Error    Smajor Sminor Az        Err  mdist   Mdist     Err        Err        Err     Quality

            2008/02/16 02:49:56.0      43.696    10.660      9.7      12    7 169  Md 2.1  7  Ml 1.6  4             ITA_NDC   00011372
                0.52    +-  0.19       3.3    1.7   21    +-  1.8    0.10   0.69   +-0.1      +-0.2      +-        m i ke

            ITALY (Valdarno inferiore)
            Sta    Dist   EvAz     Phase         Date       Time  TRes  Azim  AzRes  Slow  SRes Def  SNR        Amp   Per   Mag1   Mag2 Arr ID
            PII     0.10   285 m   Pg      2008/02/16 02:49:59.4   0.3                          T                                       00149612
            PII     0.10   285 m   Sg      2008/02/16 02:50:00.8  -0.4                          T                         Md 2.1        00149613
            CRMI    0.24    67 m   Pg      2008/02/16 02:50:01.8  -0.0                          T                                       00149614
            CRMI    0.24    67 m   Sg      2008/02/16 02:50:05.8  -0.2                          T                45   .20 Md 1.9 Ml 1.4 00149615
            BDI     0.37   353 m   Pg      2008/02/16 02:50:03.9  -0.1                          T                                       00149616
            BDI     0.37   353 m   Sg      2008/02/16 02:50:08.9  -0.8                          T                80   .29 Md 1.9 Ml 1.8 00149617
            CSNT    0.51   116 m   Pg      2008/02/16 02:50:06.3  -0.1                          T                                       00149618
            CSNT    0.51   116 m   Sg      2008/02/16 02:50:13.3  -0.6                          T                48   .36 Md 2.2 Ml 1.7 00149619
            VLC     0.50   337 m   Pg      2008/02/16 02:50:06.4   0.1                          T                                       00149620
            VLC     0.50   337 m   Sg      2008/02/16 02:50:13.8   0.0                          T                40   .86 Md 2.0 Ml 1.5 00149621
            SEI     0.62    54 m   Pg      2008/02/16 02:50:08.8   0.5                          T                83   .32 Md 2.2        00149622
            VMG     0.69    67 m   Pg      2008/02/16 02:50:10.1   0.7                          T                12   .16 Md 2.3        00149623

            STOP

            123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
            
            NOTE: the point separating two event blocks seems to be a non-standard INGV addition
                  according to the standard, event blocks are separated by two blank lines

            NOTE: some INGV bulletin files have some deviations from the standard:
                  (1) blank line between the two origin lines
                  (2) magnitude-related fields are in the wrong columns
                  (3) no blank line between last header line (DATA_TYPE BULLETIN GSE2.0) 
                      and first 'EVENT NNNNNNNN' line

            NOTE: if kwarg 'checkMessageHeader' is set to True, the header lines are not checked
                  for correctness. However, there has to be at least *one* non-blank header line.

            NOTE: the 'STOP' line at the end finishes processing (entries after it are ignored),
                  but is not required for correct processing
            
        correct illegal time components: NO
        
        kwargs: fieldIndices              - dictionary of positions of input fields, used if
                                            input file is non-standard, like older INGV bulletins 
                                            (if not set, use default)
                nopicks            = True - do not read in phase pick lines
                                            (default: False)
                minimumDataset     = True - only read basic information (save memory)
                                            (default: False)
                checkMessageHeader = True - check if correct header lines are at beginning of stream
                                            (default: False)
                authorityID               - authority ID used in publicID attributes
                                            (default: 'local')
                networkCode               - network code used for phase lines
                                            (default: 'XX')
        """
        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        if 'minimumDataset' in kwargs and kwargs['minimumDataset'] is True:
            minConf = True
        else:
            minConf = False

        if 'checkMessageHeader' in kwargs and kwargs['checkMessageHeader'] is True:
            checkMessageHeader = True
        else:
            checkMessageHeader = False

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'authorityID' in kwargs and isinstance( kwargs['authorityID'], basestring ):
            auth_id = kwargs['authorityID']
        else:
            auth_id = 'local'

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'networkCode' in kwargs and isinstance( kwargs['networkCode'], basestring ):
            network_code = kwargs['networkCode']
        else:
            network_code = 'XX'

        # if input field positions are given, replace default
        if 'fieldIndices' in kwargs and isinstance( kwargs['fieldIndices'], dict ):
            fieldIndices = kwargs['fieldIndices']
        else:

            # standard-conforming default field positions
            # NOTE: these are Python-style, i.e. first position is zero-offset, and last position is zero-offset plus one
            fieldIndices = { 'author_from': 104,
                             'author_to': 112,
                             'id_from': 114,
                             'id_to': 122,
                             'mag': ( { 'magtype_from': 71, 'magtype_to': 73, 'mag_from': 73, 'mag_to': 77,
                                        'st_cnt_from': 78, 'st_cnt_to': 80, 'magerr_from': 74, 'magerr_to': 77 },
                                      { 'magtype_from': 82, 'magtype_to': 84, 'mag_from': 84, 'mag_to': 88,
                                        'st_cnt_from': 89, 'st_cnt_to': 91, 'magerr_from': 85, 'magerr_to': 88 },
                                      { 'magtype_from': 93, 'magtype_to': 95, 'mag_from': 95, 'mag_to': 99,
                                        'st_cnt_from': 100, 'st_cnt_to': 102, 'magerr_from': 96, 'magerr_to': 99 } )
                           }
                                
        line_ctr = 0

        startMode      = True
        headerMode     = None
        newEventMode   = None
        originMode     = None
        regionMode     = None
        phaseMode      = None
        eventEndMode   = None
        
        commentMode     = None
        comment_str_xml = ''
        
        for line in istream:

            line_ctr += 1

            # if line is blank, change mode and skip

            # startMode
            if startMode is True:
                if len( line.strip() ) == 0:
                    continue
                else:
                    startMode   = False
                    headerMode  = True
                    header_line = 0
                    
            # header mode
            if headerMode is True:
                
                # header block starts with BEGIN line
                # first blank line is interpreted as end of header block
                # if no blank line after header block, detect 'EVENT' line
                if len( line.strip() ) == 0:
                    headerMode   = False
                    newEventMode = True
                    continue
                elif line.strip().startswith('EVENT'):
                    # first event block has started without separating blank line
                    # don't go to next line, fall through to new event mode block
                    headerMode   = False
                    newEventMode = True
                else:
                    header_line += 1
                    if checkMessageHeader is True:

                        # check header lines 1 (BEGIN) and 5 (DATA_TYPE) for correctness, ignore others
                        if (    ( ( header_line == 1 ) and ( upper( line.strip() ) != 'BEGIN GSE2.0' ) )
                             or ( ( header_line == 5 ) and ( upper( line.strip() ) != 'DATA_TYPE BULLETIN GSE2.0' ) ) ):
                                error_str = " importGSE2_0Bulletin - illegal GSE2.0 header line %s: %s" % ( line_ctr, line )
                                raise ValueError, error_str
                        else:
                            continue
                    else:
                        continue

            # event end mode
            if eventEndMode is True:

                if string.upper( line.strip() ) == 'STOP':

                    # end processing loop
                    break
                
                elif len( line.strip() ) == 0 or line.strip() == '.':
                    continue
                else:

                    # new regular line, new data block starts
                    eventEndMode = False
                    newEventMode = True
                    
            if newEventMode is True:
                
                # event block starts with EVENT line
                # then the origin header lines follow: ignore them
                # first blank line is interpreted as begin of first origin block
                if len( line.strip() ) == 0:
                    newEventMode = False
                    originMode   = True
                    origin_line  = 0
                    continue

                elif line.strip().startswith('EVENT'):

                    # create event
                    ev = Event()
                    ev.add( self.eventParameters, 'event' )
                    # self.eventParameters.event.append( ev )

                    # get event id 
                    try:
                        ev_id = line[6:].strip()
                    except:
                        error_str = " importGSE2_0Bulletin - event id missing in EVENT line %s: %s" % ( line_ctr, line )
                        raise ValueError, error_str

                    if minConf is not True:
                        ev.publicID = ''.join( ( 'smi:', auth_id, '/event/', ev_id ) )

                    continue

                else:
                    # ignore non-blank and non-EVENT line
                    continue
                    
            if originMode is True:

                # origin block consists of 2-line origin info pairs, separated by blank lines
                # origin block is terminated by region info line

                # NOTE: in older INGV bulletin data sets, origin line pairs can be separated with blank line

                # blank line: just ignore
                if len( line.strip() ) == 0:
                    continue

                elif len( line[0].strip() ) > 0:
                    
                    # non-blank line with non-empty first column
                    # this is either 1st event line or region line which indicates the next input block

                    # check if the first entry in line is a valid date
                    match_str = r'\d{4}/\d{2}/\d{2}\s+'
                
                    if re.match( match_str, line.strip() ):

                        origin_line = 1
                    
                    else:
                        # this line is already the region line, set region info within origin block
                        origin_line = 0
                        originMode  = False
                        regionMode  = True

                else:
                    
                    # non-blank line with empty first column: 2nd event line
                    # if current origin line is not 1, we have an error
                    # 2nd origin line cannot stand alone

                    if origin_line == 1:
                        origin_line = 2
                    else:
                        error_str = " importGSE2_0Bulletin - illegal format in origin line %s: %s" % ( line_ctr, line )
                        raise ValueError, error_str
                
                if origin_line == 1:
                    try:
                        # require time and lat/lon
                        
                        curr_year   = int( line[0:4] )
                        curr_month  = int( line[5:7] )
                        curr_day    = int( line[8:10] )
                        
                        curr_hour   = int( line[11:13] )
                        curr_minute = int( line[14:16] )
                        curr_second = float( line[17:21].strip() )

                        # ignored in QuakeML <= 1.0.2: fixed time flag

                        curr_lat    = float( line[25:33].strip() )
                        curr_lon    = float( line[34:43].strip() )

                        # ignored in QuakeML <= 1.0.2: fixed epicenter flag

                        # create origin, set publicID later
                        ori = Origin()

                        ori.time      = TimeQuantity( ( curr_year, curr_month, curr_day,
                                                        curr_hour, curr_minute, curr_second ) )
                        ori.latitude  = RealQuantity( curr_lat )
                        ori.longitude = RealQuantity( curr_lon )

                        ori.add( ev, 'origin' )

                    except:
                        print " importGSE2_0Bulletin - error in origin line %s: %s" % ( line_ctr, line )
                        continue

                    ## optional fields

                    # origin ID
                    try:    ori_id = line[fieldIndices['id_from']:fieldIndices['id_to']].strip()
                    except: ori_id = ev_id

                    if minConf is not True:
                        ori.publicID = ''.join( ( 'smi:', auth_id, '/origin/', ori_id ) )
                    
                    # depth
                    try:
                        curr_depth = float( line[47:52].strip() )
                        ori.depth  = RealQuantity( curr_depth )
                    except:
                        pass

                    # ignored in QuakeML <= 1.0.2: fixed depth flag

                    if minConf is not True:

                        # number of used phases
                        try:    
                            used_phase_cnt = int( line[56:60].strip() )
                            if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                                ori.quality   = OriginQuality()
                            ori.quality.usedPhaseCount = used_phase_cnt
                        except: 
                            pass
                            
                        # number of used stations
                        try:    
                            used_sta_cnt = int( line[61:65].strip() )
                            if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                                ori.quality   = OriginQuality()
                            ori.quality.usedStationCount = used_sta_cnt
                        except: 
                            pass
                        
                        # azimuthal gap
                        try:    
                            azimuthal_gap = float( line[66:69].strip() )
                            if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                                ori.quality   = OriginQuality()
                            ori.quality.azimuthalGap = azimuthal_gap
                        except: 
                            pass

                    ## magnitudes
                    mag_arr = []
                    for mag_ctr in xrange( 3 ):

                        try:
                            curr_mag_str = line[fieldIndices['mag'][mag_ctr]['mag_from']:fieldIndices['mag'][mag_ctr]['mag_to']].strip()
                        except:
                            curr_mag_str = ''

                        # something there?
                        if len( curr_mag_str ) > 0:
                            
                            mag = Magnitude()
                            mag.add( ev, 'magnitude' )
                            mag.setOriginAssociation( ori.publicID )
                                                        
                            try:
                                mag.mag = RealQuantity( float( curr_mag_str ) )
                            except:
                                error_str = " importGSE2_0Bulletin - illegal magnitude format in EVENT line %s: %s" % ( line_ctr, line )
                                raise ValueError, error_str

                            # magnitude type
                            try:
                                mag.type = line[fieldIndices['mag'][mag_ctr]['magtype_from']:fieldIndices['mag'][mag_ctr]['magtype_to']].strip()
                            except:
                                mag.type = 'unknown'

                            if minConf is not True:

                                # public ID
                                mag.publicID = ''.join( ( 'smi:', auth_id, '/magnitude/', ori_id, '/', str(mag_ctr+1) ) )
                                
                                # station count
                                try:
                                    mag.stationCount = \
                                        int( line[fieldIndices['mag'][mag_ctr]['st_cnt_from']:fieldIndices['mag'][mag_ctr]['st_cnt_to']].strip() )
                                except:
                                    pass

                            # fill mag_arr for later identification of magnitude instance
                            mag_arr.append( mag  )

                        else:
                            mag_arr.append( None )

                    # Author: map to creationInfo.agencyID
                    try:
                        curr_locsource = line[fieldIndices['author_from']:fieldIndices['author_to']].strip()
                        
                        if len( curr_locsource ) > 0:
                            if ( not hasattr( ori, 'creationInfo' ) ) or ( not isinstance( ori.creationInfo, CreationInfo ) ):
                                ori.creationInfo = CreationInfo()
                            ori.creationInfo.agencyID = curr_locsource
                    except:
                        pass

                    # set preferred origin
                    # TODO: set first origin as preferred
                    ev.preferredOriginID = ori.publicID

                    ## set first magnitude as preferred one
                    if len( ev.magnitude ) > 0:
                        ev.preferredMagnitudeID = ev.magnitude[0].publicID

                # everything in origin line 2 is ignored in minimum configuration    
                elif ( origin_line == 2 ) and ( minConf is not True ):

                    # rms - standard deviation of arrival time residuals
                    try:
                        curr_std_dev = float( line[5:10].strip() )
                        if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                            ori.quality   = OriginQuality()
                        ori.quality.standardError = curr_std_dev
                    except:
                        pass

                    # focal time error (seconds)
                    try:    ori.time.uncertainty = float( line[15:21].strip() )
                    except: pass

                    ## horizontal error ellipse

                    # add only if complete description is there
                    try:
                        ou = OriginUncertainty()
                        ori.originUncertainty.append( ou )

                        ou.minHorizontalUncertainty        = float( line[25:31].strip() )
                        ou.maxHorizontalUncertainty        = float( line[32:38].strip() )
                        ou.azimuthMaxHorizontalUncertainty = float( line[40:43].strip() )
                        ou.preferredDescription            = 'uncertainty ellipse'

                    except:
                        pass

                    # depth Error
                    try:    ori.depth.uncertainty = float( line[49:54].strip() )
                    except: pass

                    # minumum distance to station (degrees)
                    # NOTE: this is given in degrees in GSE2.0
                    # for QuakeML version <= 1.0.1 we have to transform to km
                    try:
                        min_dist = float( line[56:62].strip() )
                        if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                            ori.quality   = OriginQuality()
                        ori.quality.minimumDistance = min_dist * 111.0
                    except:
                        pass

                    # maximum distance to station (degrees)
                    # NOTE: this is given in degrees in GSE2.0
                    # for QuakeML version <= 1.0.1 we have to transform to km
                    try:
                        max_dist = float( line[63:69].strip() )
                        if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                            ori.quality = OriginQuality()
                        ori.quality.maximumDistance = max_dist * 111.0
                    except:
                        pass

                    # magnitude errors, mags from 1st origin line are saved in mag_arr
                    for mag_ctr, curr_mag in enumerate( mag_arr ):

                        if curr_mag is None:
                            continue
                        else:
                            try:
                                curr_mag.mag.uncertainty = \
                                    float( line[fieldIndices['mag'][mag_ctr]['magerr_from']:fieldIndices['mag'][mag_ctr]['magerr_to']].strip() )
                            except:
                                pass

                    # antype -> Origin.evaluationMode
                    # TODO: no match for 'g' (guess) in QuakeML
                    # we map 'g' to 'manual' and create a comment for origin
                    try:
                        curr_ori_mode = line[104:105].strip()
                        if curr_ori_mode in ( 'm', 'M' ):
                            ori.evaluationMode = 'manual'
                        elif curr_ori_mode in ( 'a', 'A' ):
                            ori.evaluationMode = 'automatic'
                        elif curr_ori_mode in ( 'g', 'G' ):
                            ori.evaluationMode = 'manual'

                            # comment: 'GSE2.0:antype=g'
                            oc = Comment( ''.join( ( 'GSE2.0:antype=', curr_ori_mode ) ) )
                            ori.comment.append( oc )

                    except:
                        pass

                    # TODO: loctype -> no good match in QuakeML
                    # loctype is ignored

                    # evtype -> Event.type
                    # TODO: no match in QuakeML for 'suspected' classification
                    # currently, we do not discriminate 'known' and 'suspected'
                    # if classification is 'unknown', we do not set the Event.type
                    # 'rockburst' -> 'other', 'mine explosion' -> 'quarry blast'
                    # 'experimental explosion' -> 'explosion'
                    # add comment to event with GSE2.0 classification
                    try:
                        curr_ev_remarks = line[108:110].strip()
                        if curr_ev_remarks in ( 'ke', 'se', 'KE', 'SE' ):
                            ev.type = 'earthquake'
                        elif curr_ev_remarks in ( 'kr', 'sr', 'KR', 'SR' ):
                            ev.type = 'other'
                        elif curr_ev_remarks in ( 'ki', 'si', 'KI', 'SI' ):
                            ev.type = 'induced earthquake'
                        elif curr_ev_remarks in ( 'km', 'sm', 'KM', 'SM' ):
                            ev.type = 'quarry blast'
                        elif curr_ev_remarks in ( 'kx', 'sx', 'KX', 'SX' ):
                            ev.type = 'explosion'
                        elif curr_ev_remarks in ( 'kn', 'sn', 'KN', 'SN' ):
                            ev.type = 'nuclear explosion'
                        elif curr_ev_remarks in ( 'ls', 'LS' ):
                            ev.type = 'landslide'

                        # comment: 'GSE2.0:evtype=<evtype>'
                        if len( curr_ev_remarks ) > 0:
                            oc = Comment( ''.join( ( 'GSE2.0:evtype=', curr_ev_remarks ) ) )
                            ev.comment.append( oc )

                    except:
                        pass

                elif regionMode is True:

                    # change mode to phases mode
                    regionMode     = False
                    phasesMode     = True
                    phase_line_ctr = 0
                    
                    if minConf is not True:

                        ## geographical region
                        curr_region_str = line.strip()

                        if len( curr_region_str ) > 0:
                            curr_region_str_xml = saxutils.escape( curr_region_str )

                            descr = EventDescription( curr_region_str_xml, 'region name' )
                            ev.description.append( descr )



            elif phasesMode is True:

            # Sta    Dist   EvAz     Phase         Date       Time  TRes  Azim  AzRes  Slow  SRes Def  SNR        Amp   Per   Mag1   Mag2 Arr ID
            # CRMI    0.24    67 m   Sg      2008/02/16 02:50:05.8  -0.2                          T                45   .20 Md 1.9 Ml 1.4 00149615

                # phase block starts with 1 line of header information
                # the following lines are phase lines
                # phase block and whole event section is terminated by blank line(s)

                phase_line_ctr += 1

                # if first line of block (header), skip
                if phase_line_ctr == 1:
                    continue

                else:

                    # if blank line, phase block has ended
                    if len( line.strip() ) == 0:
                        phasesMode   = False
                        eventEndMode = True
                        continue

                    elif 'nopicks' in kwargs and kwargs['nopicks'] is True:

                        # if keyword argument 'nopicks' set, skip line
                        continue

                    else:

                        ##  read regular phase line
                        try:

                            # required fields: Sta Phase DateTime
                            curr_sta_code   = line[0:5].strip()
                            curr_phase_code = line[23:30].strip()

                            curr_year       = int( line[31:35] )
                            curr_month      = int( line[36:38] )
                            curr_day        = int( line[39:41] )
                        
                            curr_hour       = int( line[42:44] )
                            curr_minute     = int( line[45:47] )
                            curr_second     = float( line[48:52].strip() )
                            
                        except:
                            print " importGSE2_0Bulletin - error in phase line %s: %s" % ( line_ctr, line )
                            continue

                        # if pick id is not provided, use phase line number
                        try:
                            pick_id = line[124:132].strip()
                            if len( pick_id ) == 0:
                                pick_id = str( phase_line_ctr )
                        except:
                            pick_id = str( phase_line_ctr )

                        # create pick object and add to event
                        pick = Pick()
                        if minConf is not True:
                            pick.publicID = ''.join( ( 'smi:', auth_id, '/event/', ori_id, '/pick/', pick_id ) )

                        base_time_utc = DateTime( curr_year,
                                                  curr_month,
                                                  curr_day,
                                                  curr_hour,
                                                  curr_minute,
                                                  0.0 )
                                                  
                        pick_time_utc = base_time_utc + DateTimeDeltaFromSeconds( curr_second )
                        pick.time = TimeQuantity( pick_time_utc )

                        # set waveform id
                        pick.waveformID = WaveformStreamID()
                        pick.waveformID.stationCode = curr_sta_code
                        pick.waveformID.networkCode = network_code
                        pick.add( ev, 'pick' )

                        # create arrival object and add to origin
                        arrv = Arrival()
                        arrv.pickID = pick.publicID
                        arrv.phase  = Phase( curr_phase_code )
                        arrv.add( ori, 'arrival' )
                        
                        # optional fields: Azim Slow Dist EvAz TRes AzRes SRes Def Arr
                        # AzRes, SRes, and Def are not present in QuakeML version <= 1.0.1
                        # Arr (arrival id) is also not present in QuakeML
                        if minConf is not True:

                            pick.phaseHint = Phase( curr_phase_code )
                            
                            try:    pick.azimuth = RealQuantity( float( line[59:64].strip() ) )
                            except: pass

                            try:    pick.slowness = RealQuantity( float( line[72:77].strip() ) )
                            except: pass

                            # distance is in degrees, needs to be converted to km
                            try:    arrv.distance = float( line[6:12].strip() ) * 111.0
                            except: pass

                            try:    arrv.azimuth = float( line[13:18].strip() )
                            except: pass

                            try:    arrv.residual = float( line[53:58].strip() )
                            except: pass
                            
                            ## Amplitude

                            # optional fields: SNR Amp Per
                            # create amplitude object if amplitude value is given
                            try:
                                curr_amp = Amplitude()
                                curr_amp.publicID = ''.join( ( 'smi:', auth_id, '/event/', ori_id, '/amplitude/', pick_id ) )

                                amp_value    = float( line[94:103].strip() )
                                curr_amp.amp = RealQuantity( amp_value )

                                curr_amp.type       = 'A'
                                curr_amp.pickID     = pick.publicID
                                curr_amp.waveformID = pick.waveformID
                                
                                curr_amp.add( ev, 'amplitude' )
                            except:
                                curr_amp = None

                            # NOTE: Amplitude.snr will be included in QuakeML version > 1.0.1
                            #try:    curr_amp.snr = float( line[88:93].strip() )
                            #except: pass

                            try:    curr_amp.period = RealQuantity( float( line[104:109].strip() ) )
                            except: pass

                            ## StationMagnitude(s), optional fields: Mag1 Mag2
                            sta_mag_indices = [ { 'mag_from': 112, 'mag_to': 116, 'magtype_from': 110, 'magtype_to': 112 },
                                                { 'mag_from': 119, 'mag_to': 123, 'magtype_from': 117, 'magtype_to': 119 } ]

                            # 1st potential station magnitude
                            for sta_mag_idx in xrange( 2 ):
                                try:
                                    sta_mag      = StationMagnitude()

                                    mag_value    = float(
                               line[sta_mag_indices[sta_mag_idx]['mag_from']:sta_mag_indices[sta_mag_idx]['mag_to']].strip() )
                                    sta_mag.mag  = RealQuantity( mag_value )
                                    sta_mag.type = \
                               line[sta_mag_indices[sta_mag_idx]['magtype_from']:sta_mag_indices[sta_mag_idx]['magtype_to']].strip()

                                    sta_mag.add( ori, 'stationMagnitude' )

                                    sta_mag.waveformID = pick.waveformID

                                    if curr_amp is not None:
                                        sta_mag.stationAmplitudeID = curr_amp.publicID
                                    
                                except:
                                    pass
                                    
                    continue

    def importOGS_HPL( self, input, **kwargs ):
        """
        import earthquake catalog data in HPL format as used by OGS
        
        data example from OGS (note: ruler shown below is not part of data)

        123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890

             9 770508 0619                                   0.49                                                                                                      2
             9 BUA                eP   0619 24.90                                                   FL   80 0.5  iS   28.80                        BUA     gg
             9 COLI               eP   0619 25.90                                                   FL   70 0.5  iS   30.00                        COLI    gg
        ^MOGGIO UDINESE (FRIULI)

            10 770508 1659  8.76 46-21.69  13- 9.51   1.62   0.78  5 16 317 1 0.12  0.3  2.1 C B/D 0.56 10  6 0.00 0.10  0  0.0  0.0  2  0.8  0.2 9  0.3  0.1          3
            10 BAD   15.6 155  96    4 1659  0.00******  2.68  0.00        0.00   0  0  0.00 0      FLD          iS   13.30  4.54 -0.23   0.8      BAD      g
            10 BUA   16.3 189  96 eP   1659 11.60  2.84  2.80  0.00  0.04  1.25   0  0  0.00 0      FLD 100 0.7  iS   13.70  4.94 -0.05   1.2      BUA     gg
            10 COLI  30.5 147  93 eP   1659 14.10  5.34  5.23  0.00  0.11  0.83   0  0  0.00 0      FLD 100 0.8  iS   18.20  9.44  0.13   0.8      COLI    gg
        ^ZAGARIE (SLOVENIA)

            11 770508 22 9 27.33 46- 5.83  15-19.11  18.07   2.26  5 61 336 1  .14  3.4  1.6 D C/D 7.15 10  6  .00  .13  0   .0   .0  1  2.3   .0 4  2.6  2.1          3
            11 LJU   61.5 264 106 eP   22 9 38.20 10.87 10.96   .00  -.09  1.35   0  0   .00 0      SLD          eS   46.90 19.57   .07   1.3      LJU     gg
            11 CEY   79.8 240 103 eP   22 9 41.50 14.17 13.99   .00   .18   .84   0  0   .00 0      SLD          eS   52.00 24.67  -.24    .8      CEY     gg
            11 TRI  128.1 250  59    4 2210   .00 32.67 21.10   .00         .00   0  0   .00 0      SLD 650 2.3  eS    5.00 37.67   .11    .7      TRI      g                                                   

        Note:
            (1) lines are very long in this format (> 200 chars)
            (2) there are localized and non-localized EQs in the catalog. blocks of localized events
                have a starting line that begins with '^' in the first column and gives the region.
                a blank line (whitespace only) follows immediately
                after that, a hypocenter line and phase lines follow 
                non-localized event blocks have the same sequence of hypocenter and phase lines,
                but have no starting ^-line with region
            (3) all lines which are not ^-lines with region info start with a sequence number that takes the
                first 6 columns, right-adjusted (!)
            (4) there can be comment lines after last phase line that start with '*'

        correct illegal time components: NO
        
        kwargs: nopicks            = True - do not read in phase pick lines
                minimumDataset     = True - only read basic information (save memory)
                authorityID               - authority id used in public ids, default: 'local'
                networkCode               - network code used for phase lines, default: 'XX'
        """
        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        if 'minimumDataset' in kwargs and kwargs['minimumDataset'] is True:
            minConf = True
        else:
            minConf = False

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'authorityID' in kwargs and isinstance( kwargs['authorityID'], basestring ):
            auth_id = kwargs['authorityID']
        else:
            auth_id = 'OGS'

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'networkCode' in kwargs and isinstance( kwargs['networkCode'], basestring ):
            network_code = kwargs['networkCode']
        else:
            network_code = 'XX'
                                
        line_ctr = 0

        newEventMode   = None
        originMode     = None
        phaseMode      = None
        skipMode       = True
        
        comment_str_xml = ''
        
        for line in istream:

            line_ctr += 1

            # check if it's a new localized event
            if line[0] == '^':

                newEventMode = True
                skipMode     = False

                # get region string
                region_str_xml = saxutils.escape( line[1:].strip() )
                continue

            elif line[0] == '*':

                # comment line, skip
                continue

            else:

                if skipMode is True :
                    continue

                elif newEventMode is True:

                    # skip one blank line after region info
                    originMode   = True
                    newEventMode = False
                    continue
                    
                elif originMode is True:

                    # get hypocenter data

                    # create event
                    ev = Event()
                    ev.add( self.eventParameters, 'event' )
                    # self.eventParameters.event.append( ev )

                    # get event sequence number of data set
                    try:
                        ev_seq = int( line[0:6].strip() )
                    except:
                        error_str = " importOGS_HPL - no valid sequence number in origin line %s: %s" % ( line_ctr, line )
                        raise ValueError, error_str
                    
                    try:
                        
                        # require time, lat/lon
                        
                        curr_year_str = line[7:9]
                        curr_month    = int( line[9:11] )
                        curr_day      = int( line[11:13] )

                        # catalog starts in 1977
                        # make years >= 77 -> 19XX, < 77 -> 20XX
                        if int( curr_year_str ) >= 77:
                            curr_year = int( ''.join( ( '19', curr_year_str ) ) )
                        else:
                            curr_year = int( ''.join( ( '20', curr_year_str ) ) )
                            
                        curr_hour   = int( line[14:16] )
                        curr_minute = int( line[16:18] )
                        curr_second = float( line[19:24].strip() )

                        # latitude and longitude are given as degrees and decimal minutes
                        curr_lat_deg = float( line[25:27].strip() )
                        curr_lat_min = float( line[28:33].strip() )
                        curr_lat     = curr_lat_deg + ( curr_lat_min / 60.0 )
                        
                        curr_lon_deg = float( line[35:37].strip() )
                        curr_lon_min = float( line[38:43].strip() )
                        curr_lon     = curr_lon_deg + ( curr_lon_min / 60.0 )
                        
                        # create origin, set publicID later
                        ori = Origin()

                        ori.time      = TimeQuantity( ( curr_year, curr_month, curr_day,
                                                        curr_hour, curr_minute, curr_second ) )
                        ori.latitude  = RealQuantity( curr_lat )
                        ori.longitude = RealQuantity( curr_lon )
                        
                        ori.add( ev, 'origin' )

                    except:
                        print " importOGS_HPL - error in origin line %s: %s" % ( line_ctr, line )
                        skipMode = True
                        continue

                    ## optional fields
                    
                    # depth
                    try:
                        curr_depth = float( line[45:50].strip() )
                        ori.depth  = RealQuantity( curr_depth )
                    except:
                        pass

                    ## origin/event id is not provided in HPL format
                    # get origin id from date and time, show as ISO datetime with 2 decimal places for seconds
                    ori_id = mxDateTime2ISO( ori.time.value.datetime, secondsdigits = 2 )
                    
                    if minConf is not True:
                        ori.publicID = ''.join( ( 'smi:', auth_id, '/origin/', ori_id ) )

                    ## set preferred origin
                    ev.preferredOriginID = ori.publicID
                        
                    ## magnitude
                    try:
                        curr_mag = float( line[52:57].strip() )
                        
                        mag = Magnitude()
                        mag.add( ev, 'magnitude' )
                        mag.setOriginAssociation( ori.publicID )

                        # duration magnitude (local?)
                        mag.mag  = RealQuantity( curr_mag )
                        mag.type = 'Md'

                        if minConf is not True:

                            mag.publicID = ''.join( ( 'smi:', auth_id, '/magnitude/', ori_id ) )

                            # station count
                            try:
                                mag.stationCount = int( line[125:127].strip() )
                            except:
                                pass

                    except:
                        pass

                    ## set first magnitude as preferred one
                    if len( ev.magnitude ) > 0:
                        ev.preferredMagnitudeID = ev.magnitude[0].publicID
                            
                    if minConf is not True:
                        ev.publicID = ''.join( ( 'smi:', auth_id, '/event/', ori_id ) )

                        ## geographical region
                        descr = EventDescription( region_str_xml, 'region name' )
                        ev.description.append( descr )
                        
                        # number of used phases
                        try:    
                            used_phase_cnt = int( line[58:60].strip() )
                            if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                                ori.quality   = OriginQuality()
                            ori.quality.usedPhaseCount = used_phase_cnt
                        except: 
                            pass

                        # azimuthal gap
                        try:    
                            azimuthal_gap = float( line[64:67].strip() )
                            if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                                ori.quality   = OriginQuality()
                            ori.quality.azimuthalGap = azimuthal_gap
                        except: 
                            pass

                        # standardError (rms)
                        try:
                            curr_rms_residual_phases = float( line[70:74].strip() )
                            if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                                ori.quality   = OriginQuality()
                            ori.quality.standardError = curr_rms_residual_phases
                        except: 
                            pass
                            
                        # horizontal error (km)
                        try:
                            che = float( line[74:79].strip() )
                            ou = OriginUncertainty()
                            ou.horizontalUncertainty = che
                            ou.add( ori, 'originUncertainty' )
                        except:
                            pass

                        # depth error (km)
                        try:    ori.depth.uncertainty = float( line[79:84].strip() )
                        except: pass
                        
                        # number of associated stations
                        try:    
                            ass_sta_cnt = int( line[99:101].strip() )
                            if ( not hasattr( ori, 'quality' ) ) or ( not isinstance( ori.quality, OriginQuality ) ):
                                ori.quality   = OriginQuality()
                            ori.quality.associatedStationCount = ass_sta_cnt
                        except: 
                            pass

                    phaseMode      = True
                    phase_line_ctr = 0
                    originMode     = False
                    continue

                elif phaseMode is True:

                    ## read phase line

                    # if keyword argument 'nopicks' set, skip line
                    if 'nopicks' in kwargs and kwargs['nopicks'] is True:
                        continue

                    phase_line_ctr += 1
                    
                    # get sequence number
                    try:
                        curr_ev_seq = int( line[0:6].strip() )
                    except:
                        error_str = " importOGS_HPL - no valid sequence number in phase line %s: %s" % ( line_ctr, line )
                        raise ValueError, error_str
                    
                    # check if there's a change in sequential number
                    # this would indicate a location line of a new non-localized event
                    # -> skip subsequent lines until next '^' localized event line 
                    if ev_seq != curr_ev_seq:
                        skipMode = True
                        continue

                    ## required: station code, phase block, arrival hour, minute, second
                    try:
                        curr_sta_code    = line[7:11].strip()
                        phase_block      = line[26:29].strip()

                        curr_hour        = int( line[31:33].strip() )
                        curr_minute      = int( line[33:35].strip() )
                        curr_second      = float( line[36:41].strip() )

                    except:
                        print " importOGS_HPL - error in phase line %s: %s" % ( line_ctr, line )
                        continue

                    # no pick id provided: use phase line number
                    pick_id = str( phase_line_ctr )

                    # set base time (hour/minute) from which pick times are measured
                    # curr_hour and curr_minute are trusted to be in valid range
                    base_time_utc = DateTime( curr_year,
                                              curr_month,
                                              curr_day,
                                              curr_hour,
                                              curr_minute,
                                              0.0 )

                    # sometimes first (P) pick seems to be invalid, then seconds are set to 0.0
                    # and a sequence of '*' chars follow
                    if not ( curr_second == 0.0 and line[41] == '*' ):
                        
                        # create pick object and add to event
                        pick = Pick()
                        if minConf is not True:
                            pick.publicID = ''.join( ( 'smi:', auth_id, '/event/', ori_id, '/pick/', pick_id ) )

                        pick_time_utc = base_time_utc + DateTimeDeltaFromSeconds( curr_second )
                        pick.time     = TimeQuantity( pick_time_utc )

                        # set waveform id
                        pick.waveformID = WaveformStreamID()
                        pick.waveformID.stationCode = curr_sta_code
                        pick.waveformID.networkCode = network_code
                        pick.add( ev, 'pick' )

                        # create arrival object and add to origin
                        arrv           = Arrival()
                        arrv.pickID    = pick.publicID
                        arrv.phase     = Phase( 'P' )
                        arrv.add( ori, 'arrival' )

                        ## optional fields
                        
                        if minConf is not True:

                            pick.phaseHint = Phase( 'P' )

                            try:    phase_onset = phase_block[0]
                            except: phase_onset = ''

                            try:    phase_polarity = phase_block[2]
                            except: phase_polarity = ''
                        
                            # epicentral distance (km)
                            try:    arrv.distance = float( line[12:17].strip() )
                            except: pass

                            # azimuth
                            try:    arrv.azimuth = float( line[18:21].strip() )
                            except: pass

                            # angle of incidence (= slowness azimuth)
                            try:    pick.azimuth = RealQuantity( float( line[22:25].strip() ) )
                            except: pass

                            # P phase onset (emergent/impulsive)
                            if phase_onset in ( 'e', 'E' ):
                                pick.onset = 'emergent'
                            elif phase_onset in ( 'i', 'I' ):
                                pick.onset = 'impulsive'

                            # P polarity
                            # NOTE: change for QuakeML version > 1.0.1
                            if phase_polarity == '+':
                                pick.polarity = 'up'
                            elif phase_polarity == '-':
                                pick.polarity = 'down'

                            # P-arrival residual (sec)
                            try:    arrv.residual = float( line[60:65].strip() )
                            except: pass

                            # weight in hypocentral solution for this arrival
                            try:    arrv.weight = float( line[67:71].strip() )
                            except: pass

                            ## duration/amplitude for magnitude
                            # QuakeML TimeWindow requires an absolute point in time for duration
                            # we use pick time of P phase
                            # measured duration goes into 'end' component of TimeWindow, 'begin' is set to zero
                            amp = None
                            try:
                                amp            = Amplitude()

                                tw             = TimeWindow( 0.0, float( line[96:99].strip() ) )
                                tw.reference   = pick.time.value
                                amp.timeWindow = tw
                                        
                                amp.add( ev, 'amplitude' )
                                        
                            except:
                                pass

                            if amp is not None:
                                
                                amp.publicID   = ''.join( ( 'smi:', auth_id, '/event/', ori_id, '/amplitude/', pick_id ) )
                                
                                amp.type       = 'END'
                                amp.pickID     = pick.publicID
                                amp.waveformID = pick.waveformID

                            ## duration magnitude at station
                            sta_mag = None
                            try:
                                sta_mag     = StationMagnitude()
                                sta_mag.mag = RealQuantity( float( line[100:103].strip() ) )

                                sta_mag.add( ori, 'stationMagnitude' )
                                    
                            except:
                                pass

                            if sta_mag is not None:

                                sta_mag.type        = 'Md'
                                sta_mag.waveformID  = pick.waveformID

                                try:    sta_mag.stationAmplitudeID = amp.publicID
                                except: pass

                    ## S arrival
                    try:
                        phase_block   = line[105:107].strip()
                        curr_second_s = float( line[110:115].strip() )
                    except:
                        continue

                    # no pick id provided: use phase line number
                    phase_line_ctr += 1
                    pick_id = str( phase_line_ctr )

                    # create pick object and add to event
                    pick = Pick()
                    if minConf is not True:
                        pick.publicID = ''.join( ( 'smi:', auth_id, '/event/', ori_id, '/pick/', pick_id ) )

                    pick_time_utc = base_time_utc + DateTimeDeltaFromSeconds( curr_second_s )
                    pick.time     = TimeQuantity( pick_time_utc )

                    # set waveform id
                    pick.waveformID = WaveformStreamID()
                    pick.waveformID.stationCode = curr_sta_code
                    pick.waveformID.networkCode = network_code
                    pick.add( ev, 'pick' )

                    # create arrival object and add to origin
                    arrv = Arrival()
                    arrv.pickID = pick.publicID
                    arrv.phase  = Phase( 'S' )
                    arrv.add( ori, 'arrival' )
                    
                    if minConf is not True:

                        pick.phaseHint = Phase( 'S' )

                        try:    phase_onset = phase_block[0]
                        except: phase_onset = ''
                        
                        # S phase onset (emergent/impulsive)
                        if phase_onset in ( 'e', 'E' ):
                            pick.onset = 'emergent'
                        elif phase_onset in ( 'i', 'I' ):
                            pick.onset = 'impulsive'

                        # residual of S arrival (sec)
                        try:    arrv.residual = float( line[122:127].strip() )
                        except: pass

                        # weight in hypocentral solution for S arrival
                        try:    arrv.weight = float( line[130:133].strip() )
                        except: pass
                        
                    continue

    def importShadow2000( self, input, **kwargs ):
        """

        NOTE: WORK IN PROGRESS

        import earthquake catalog data in Shadow 2000 (aka Hypoinverse-2000 archive) format
        this is used by the Northern California Seismic Network (NCSN)
        see documentation in http://www.ncedc.org/ftp/pub/doc/ncsn/shadow2000.pdf

        data example from NCSN (note: ruler shown below is not part of data)

                 10        20        30        40        50        60        70        80        90        100       110       120       130       140       150       160                          
        12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890

        201001010241047539 2332121 1489 1826    10198  2   7 8249  8422935  59149     33    2  54  64  9      80     9       D 13 D149 80         71328725D149  80        2FNC01
        ORV  BK  HHN    4201001010241             1375IS 1 -54           0       0 284120          310             0J  --       
        ORV  BK  HHZ IPU1201001010241  989 -22 64        0                   0     284120 7       5310 91     65    JD --       
        ABJ  NC  EHZ IPD0201001010241  948  -3138 1316IS 1  -6          69 -16 -28 252123 0       7169160    714 271JD --       
        AOH  NC  EHZ IPD0201001010241  759   3138  977IS 1   2          69  -9 -16  16174 0       9210149    406 516JD --       
        OBHB NC  EHZ IPD1201001010241 1093  -6 69        0                   7     345114 2       6327149    126    JD --       
        OCR  NC  EHZ EPD2201001010241 1669  44  0        0                  32     686 98 2      11321153      0    JD --       
        OGO  NC  EHZ EP 3201001010241 1192 -27  3        0                   6     430108 2       7313122      0    JD --       
        OHC  NC  EHZ IPU0201001010241  903   7138        0                 -22     212129 0       8253157    346    JD --       
        OSU  NC  EHZ IPU0201001010241 1434   1138        0                  62     537103 2       9255120    581    JD --       
        OWY  NC  EHZ IPD0201001010241  901  -9138        0                 -15     218128 0       6289133    273    JD --       
        MGL  WR  EHZ IPD0201001010241 1377   5138        0                  -2     540103          330       699    J  --       
                                                                        71328725

        Notes:
            (1) the first (summary) line of an event block is very long (168 chars)
            (2) station/phase lines have 120 chars
            (3) after the phase lines, a terminator line with the EQ ID follows. The ID is at
                position 63-72 (10 chars)
            (4) there is no separator line between event blocks, but we will ignore whitespace-only
                lines between blocks
            (5) shadow lines (starting with '$') are ignored

        kwargs: 
            nopicks            = True - do not read in phase pick lines
            minimumDataset     = True - only read basic information (save memory)
            authorityID               - authority id used in public ids, default: 'local'
        """
        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        if 'minimumDataset' in kwargs and kwargs['minimumDataset'] is True:
            minConf = True
        else:
            minConf = False

        # TODO: needs to be encoded for XML, check if valid format for XML attribute
        if 'authorityID' in kwargs and isinstance( kwargs['authorityID'], basestring ):
            auth_id = kwargs['authorityID']
        else:
            auth_id = 'NCSN'
                                
        line_ctr = 0

        newEventMode   = None
#        originMode     = None
        phaseMode      = None
        terminatorMode = None
#        skipMode       = True
        
        comment_str_xml = ''
        
        for line_ctr, line in enumerate( istream ):

            line_length = len( line )

            # check if line is whitespace only, skip
            if len( line.strip() ) == 0:
                continue

            elif len( line ) == 168:

                # event summary line
                newEventMode = True

                # get hypocenter data

                # create event
                ev = Event()
                ev.add( self.eventParameters, 'event' )
                # self.eventParameters.event.append( ev )

            elif len( line ) == 72:

                # terminator line
                phaseMode = False
                terminatorMode = True

                # get ID string
                region_str_xml = saxutils.escape( line[1:].strip() )

            elif len( line ) == 120:

                # phase line
                newEventMode = False
                phaseMode = True

            else:
                # no valid line
                pass

    def importSHAREPotsdamCSV( self, input, **kwargs ):
        """
        Import EQ catalog in format provided by Potsdam group for SHARE, 
        (SHEEC, version 20110808)
        This is based on an Microsoft Excel spreadsheet that has been
        exported to CSV, with "|" as separator character, and double quotes
        around text fields.

        We will use moment magnitudes in field "Mw".
        Depth is optional.
        Skip events with no magnitude.
        Use only mainshocks.
        
        There are additional fields not present in QuakeML data model:
        
        ASA: ID of area zone event falls into
        ASR: ID of area zone event is attributed to for activity computation
        Comp: can be 'c', or empty. event magnitude below or above completeness level
        Main: 1: mainshock, 0: triggered event
        
        Example:
        
        "En"|"Year"|"Mo"|"Da"|"Ho"|"Mi"|"Se"|"Ax"|"Reg"|"Lat"|"Lon"|"LatUnc"|"LonUnc"|"TEpi"|"H"|"HUnc"|"TH"|"Io"|"TIo"|"Mw"|"MwUnc"|"TMw"|"Root"|"Nmdp"|"Ix"|"CCode"|"CM"|"CTM"|"ASA"|"ASR"|"TA
S"|"CSZ"|"CSZname"|"Comp"|"Main"|"event_id"
        210|1000|"03"|"29"||||"Saint-Amand"|"SCR"|50.183|4.237|50.0|50.0|"bx"||||"4-5"|"bx"|3.73|0.30|"bx"|"Alexandre, 1990"|4|"NC"||||124|124|"A"|"E"|"FRAB"||1|"210"
        500|1005||||||"AREZZO"|"APD"|43.463|11.879|30.0|30.0|"bx"||||"7"|"bx"|5.20|0.33|"wm"|"Castelli et al., 1996"|4|"7-8"|"CPTI004"|||307|307|"A"|"K"|"CITA"||1|"500"
        600|1005||||||"CASSINO"|"APD"|41.488|13.831|30.0|30.0|"bx"||||"7"|"bx"|5.20|0.33|"wm"|"Figliuolo & Mar., 2002"|5|"7"|"CPTI004"|||308|308|"A"|"K"|"CITA"||1|"600"
        800|1009||||||"Offshore Portugal"|"TSZ"|36.000|-10.700|||"cat"||||||||"nd"|"LNEC, 1986"|||"MAMV001"|||268|268|"A"|"B"|"OFFP"||0|"800"

        OBSOLETE:
        "year","mo","d","h","mi","lat","lon","depth","Int","M","type","Mw","Ref"
        1900,1,8,14,45,43.7,16.7,,7,,,4.9,"HHM"
        1900,1,14,9,53,43.5,27,10,7,5.9,"w",5.9,"Onc"
        1900,1,21,0,39,46.17,14.5,7,4.5,,,3.5,"ZivS"

        correct illegal time components: YES (was: NO)
        """

        CSV_DELIMITER = '|'
        CSV_QUOTECHAR = '"'
        
        YEAR_IDX = 1
        MONTH_IDX = 2
        DAY_IDX = 3
        HOUR_IDX = 4
        MIN_IDX = 5
        SEC_IDX = 6
        LAT_IDX = 9
        LON_IDX = 10
        LAT_ERR_IDX = 11
        LON_ERR_IDX = 12
        DEPTH_IDX = 14
        DEPTH_ERR_IDX = 15
        MAG_IDX = 19
        MAG_ERR_IDX = 20
        ASA_IDX = 28
        ASR_IDX = 29
        COMP_IDX = 33
        MAIN_IDX = 34
        
        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        reader = csv.reader( istream, delimiter=CSV_DELIMITER, 
            quotechar=CSV_QUOTECHAR )

        # over lines (= events) in input stream
        for line_ctr, line in enumerate( reader ):

            # skip first header line
            if line_ctr == 0:
                continue

            # skip if event is not a mainshock
            if int(line[MAIN_IDX]) != 1:
                continue
            
            # create event
            ev = Event()
              
            # create origin
            ori = Origin()
            ori.add( ev, 'origin' )
            
            ori.longitude = RealQuantity( float( line[LON_IDX] ) )
            ori.latitude  = RealQuantity( float( line[LAT_IDX] ) )

            # optional
            try:
                ori.depth = RealQuantity( float( line[DEPTH_IDX] ) )
            except ValueError:
                ori.depth = RealQuantity( numpy.nan )
            
            # month, day, hour and minute are not always defined, set to 0 or 1 
            # if not available
            try:
                month = int( line[MONTH_IDX] )
            except ValueError:
                month = 1
                
            try:
                day = int( line[DAY_IDX] )
            except ValueError:
                day = 1

            try:
                hour = int( line[HOUR_IDX] )
            except ValueError:
                hour = 0

            try:
                minute = int( line[MIN_IDX] )
            except ValueError:
                minute = 0

            try:
                sec = float( line[SEC_IDX] )
            except ValueError:
                sec = 0.0
                
            # get time components
            ori.time = TimeQuantity( ( int( line[YEAR_IDX] ), month, day, 
                hour, minute, sec ) )
            
            # create magnitude
            mag = Magnitude()
            mag.add( ev, 'magnitude' )
            
            try:
                mag.mag = RealQuantity( float( line[MAG_IDX] ) )
            except ValueError:
                continue
            
            mag.setOriginAssociation( ori.publicID )
            
            # set preferred origin and magnitude
            ev.preferredOriginID    = ori.publicID
            ev.preferredMagnitudeID = mag.publicID
            
            # add event if everything went well
            ev.add( self.eventParameters, 'event' )

    def importSHEEC( self, input, input_format='csv', **kwargs ):
        """
        Import EQ catalog in SHEEC CSV or TXT format as provided by 
        Jochen Wössner for SHARE.
        
        input_format can be 'csv' (default) or 'txt'
        
        CSV as in catalog of 2011-05-03:
        
        lon,lat,depth / km,year,month,day,minute,hour,second,Mw,Event_ID
        4.237,50.183,NaN,1000,3,29,0,0,0,3.73,210
        11.879,43.463,NaN,1005,1,1,0,0,0,5.19,500
        13.831,41.488,NaN,1005,1,1,0,0,0,5.19,600
        
        Separator char: ','
        Missing values can be 'NaN'
        First line is column headers.

        TXT (ASCII rows/columns) as in catalog of 2011-05-30:
        
        Long      Lat       Depth Year   Mo  Day    Ho   Mi  Sec    Mw         Event_ID  
          4.2370  50.1830     NaN 1000    3   29     0    0  0.0000 3.73            210
        -10.4400  39.2700   57.00 2006   12   31    14   16 43.0000 4.40 20061231214500
         -4.0900  35.0200    7.00 2006   12   31    14   55 52.0000 4.20 20061231214501

        
        Missing values are 'NaN'
        First line is column headers.
        
        We will interpret magnitudes as moment magnitudes (Mw).
        Depth is optional.
        Skip events with no magnitude.

        correct illegal time components: YES
        """

        CSV_DELIMITER = ','
        
        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        if input_format.lower() == 'csv':
            reader = csv.reader( istream, delimiter=CSV_DELIMITER )
        elif input_format.lower() == 'txt':
            reader = istream
        else:
            raise ValueError, \
                "importSHEEC: unsupported input file format %" % input_format

        # over lines (= events) in ZMAP input stream
        for line_ctr, line in enumerate( reader ):

            # skip first header line
            if line_ctr == 0:
                continue

            if input_format.lower() == 'csv':
                parameters = line
            else:
                parameters = line.strip().split()

            # create event
            ev = Event()
              
            # create origin
            ori = Origin()
            ori.add( ev, 'origin' )
            
            ori.longitude = RealQuantity( float( parameters[0].strip() ) )
            ori.latitude  = RealQuantity( float( parameters[1].strip() ) )

            # depth: optional
            if parameters[2].strip().lower() != 'nan':
                try:
                    ori.depth = RealQuantity( float( parameters[2].strip() ) )
                except Exception:
                    pass
            
            # if day, hour, minute, and second are not defined, set to 0 or 1
            # NOTE: assumes that these components are not NaN
            try:
                month = int( parameters[4].strip() )
            except Exception:
                month = 1

            try:
                day = int( parameters[5].strip() )
            except Exception:
                day = 1

            try:
                hour = int( parameters[6].strip() )
            except Exception:
                hour = 0

            try:
                minute = int( parameters[7].strip() )
            except Exception:
                minute = 0

            try:
                second = float( parameters[8].strip() )
            except Exception:
                second = 0.0

            # dirty fix for special case of 31th June
            if month == 6 and day == 31:
                month = 7
                day = 1

            # get time components
            timeCorrection = fixTimeComponents( hour, minute, second )

            try:
                focal_time = DateTime( int( parameters[3].strip() ),
                                       month,
                                       day,
                                       timeCorrection['component'][0],
                                       timeCorrection['component'][1],
                                       timeCorrection['component'][2] )
            except Exception, e:
                error_str = "e: %s, %s" % ( e, ( int( parameters[3].strip() ), 
                    int( parameters[4].strip() ), day, hour, minute, second ) )
                print error_str
                continue

            focal_time = adjustDateTime( timeCorrection['increaseFlag'], 
                focal_time )
            ori.time = TimeQuantity( focal_time )

            # magnitude: required, skip event if not given
            mag = Magnitude()
            mag.add( ev, 'magnitude' )

            if parameters[9].strip().lower() == 'nan':
                continue

            try:
                mag.mag = RealQuantity( float( parameters[9].strip() ) )
            except Exception:
                continue

            mag.setOriginAssociation( ori.publicID )
            
            # set preferred origin and magnitude
            ev.preferredOriginID    = ori.publicID
            ev.preferredMagnitudeID = mag.publicID

            # add event if everything went well
            ev.add( self.eventParameters, 'event' )
 
    def exportAtticIvy( self, output, **kwargs ):
        """ 
        Export catalog to format required by Roger Musson's AtticIvy code.

        One line per event (in this implementation).

        data example (note: ruler shown below is not part of data)

            0         1         2         3         4         5
            123456789012345678901234567890123456789012345678901234567890

            YYYY MM DD HH IIAABB PPPPPP LLLLLLL  EE RRR FFF WWWWKKKSSSSS
            1690  8 27 20  0 1 1  51.83   -4.38  00 4.3 0.0 1.00  0
            1727  7 19  4  0 1 1  51.57   -3.76  00 4.8 0.0 1.00 25
            1734 10 25  3 50 1 1  50.20   -0.70  00 4.1 0.0 1.00 14

        Fields: YYYY, MM, DD, HH, II    date/time components, no seconds used
                AA, BB                  ?
                PPPPPP, LLLLLLL         latitude, longitude
                EE                      horizontal error of epicenter in km, 00 if not known
                RRR                     magnitude
                FFF                     error of magnitude, 0.0 if not known
                WWWW                    weight of this location 
                                        (can have multiple locations per event)
                KKK                     hypocenter depth in km, 0 if not known
                SSSSS                   ?
        """
        
        header_line = \
            "YYYY MM DD HH IIAABB PPPPPP LLLLLLL  EE RRR FFF WWWWKKKSSSSS"

        if isinstance( output, basestring ):
            ostream = writeQPData( output, **kwargs )
        else:
            ostream = output

        ostream.write("%s\n" % header_line)

        for ev in self.eventParameters.event:

            # check if event has preferred origin and coordinates, 
            # otherwise skip
            try:
                ori = ev.getPreferredOrigin()
                curr_lon = ori.longitude.value
                curr_lat = ori.latitude.value
            except:
                continue

            try:
                mag = ev.getPreferredMagnitude()
                mag_value = mag.mag.value
            except:
                continue

            # if origin has no depth, set depth column to zero
            if hasattr( ori, 'depth' ) and ori.depth is not None:
                depth_value = ori.depth.value
            else:
                depth_value = 0.0

            try:
                depth_value = int(ori.depth.value)
            except Exception:
                depth_value = 0

            hzErrorFound = False
            
            # look if explicit horizontal error is given in OriginUncertainty object
            # this overrides possible separate lat/lon errors
            if len( ori.originUncertainty ) > 0:
                ou = ori.originUncertainty[0]

                if hasattr( ou, 'horizontalUncertainty' ):
                    try:
                        horizontal_uncertainty_value = ou.horizontalUncertainty
                        hzErrorFound = True
                    except:
                        pass

            # if no explicit horizontal error is given, compute horizontal error from lat/lon errors
            if hzErrorFound is False:

                if ( hasattr( ori.longitude, 'uncertainty' ) and hasattr( ori.latitude, 'uncertainty' ) ):

                    try:
                        curr_lon_err = ori.longitude.uncertainty
                        curr_lat_err = ori.latitude.uncertainty
                        horizontal_uncertainty_value = math.sqrt( 
                            math.pow( curr_lat_err * 111.0, 2 ) +
                            math.pow( curr_lon_err * math.cos(
                                curr_lat * math.pi/180.0) * 111.0, 2 ) )
                        hzErrorFound = True
                    except:
                        pass


            # set horizontal error to zero, if not given
            if hzErrorFound is False:
                horizontal_uncertainty_value = 0.0

            # set magnitude error to zero, if not given
            if ( hasattr( mag.mag, 'uncertainty' ) and ( 
                mag.mag.uncertainty is not None ) ):
                magnitude_uncertainty_value = mag.mag.uncertainty
            else:
                magnitude_uncertainty_value = 0.0

            line_arr = ( '%4i' % ori.time.value.datetime.year,
                         '%3i' % ori.time.value.datetime.month,
                         '%3i' % ori.time.value.datetime.day,
                         '%3i' % ori.time.value.datetime.hour,
                         '%3i' % ori.time.value.datetime.minute,
                         ' 1 1',
                         '%7.2f' % ori.latitude.value,
                         '%8.2f' % ori.longitude.value,
                         '  %02i' % horizontal_uncertainty_value,
                         ' %3.1f' % mag_value,
                         ' %3.1f' % magnitude_uncertainty_value,
                         ' 1.00',
                         '%3i' % depth_value,
                         '     ' )

            ostream.write( "%s\n" % ''.join( line_arr ) )
        
        ostream.close()



    def importHypoInverse( self, input, **kwargs ):
        """
        import NCSN event/phase data as obtained via the website:
        http://www.ncedc.org/ncedc/catalog-search.html
        originally written by Michael Lewis, USC, 2009(?)

            sample of three events below (additional '            ' at the beginning of every line added)

            201001010109535538 4531122 4303  247     9 67  2   8 7074 124317 6  50 20     34    0  50 119  9      40    17       D 11 D 20 40         71328690D 20  40        2FNC01
            DES  BG  DPZ IPD0201001010109 5407  -7196        0                 -10      21135 5       6 47 81    686    JD --       
            EPR  BG  DPE    4201001010109             5545ES 2  62           0       0  23132          113             0J  --       
            EPR  BG  DPZ EPU2201001010109 5420  -7 39        0                   0      23132 5       3113 23     82    JD --       
            FNF  BG  DPZ IPD1201001010109 5446 -11 98        0                  -8      46108 5       3292 30    476    JD --       
            MNS  BG  DPZ EPU2201001010109 5457   0 39        0                   0      41111 6      49  0278     91    JD --       
            PFR  BG  DPZ IPU1201001010109 5415 -12 98        0                  -2      24130 5       3263 22    453    JD --       
            SSR  BG  DPZ IPU1201001010109 5413   3 98        0                  -9      17143 5       4162 40    553    JD --       
            GAXB NC  EHZ IPD1201001010109 5473  10 98        0                 -30      60 97 2       2215 11    720    JD --       
            GBG  NC  EHZ EPU2201001010109 5528   6 39        0                   6      72 94 2       2 24 20    323    JD --       
            GCR  NC  EHZ IPU0201001010109 5406   8196        0                 -25      20137 0       3359 60    613    JD 02       
            GSG  NC  EHZ EP 3201001010109 5574 -62  0        0                  24     124 72 2       2  3 -7      0    JD 02       
                                                                            71328690
            201001010128404138 4945122 4888  219    15 90  0   422770  35337 6  26 29     21    0  25  33 12      40    27       D 15 D 29 40         71328695D 29  40        2FNC01
            AL3  BG  DPZ EP 3201001010128 4133   0 25        0                   0      37111 7       3275 23     19    JD --       
            AL4  BG  DPZ IPD1201001010128 4116   4125        0                   0      24127 5       3311 26    422    JD --       
            AL5  BG  DPZ EP 3201001010128 4155  -4 25        0                   0      50101 5       3290 35     26    JD --       
            CLV  BG  DPZ IPU1201001010128 4113   0125        0                  -2      26124 7       3 53 23    283    JD --       
            DRK  BG  DPZ EPD2201001010128 4119 -15 50        0                  -6      41108 7       4165 44     48    JD --       
            FUM  BG  DPZ EP 3201001010128 4127 -13 25        0                  -2      42108 7       5146 67     10    JD --       
            SB4  BG  DPZ IPD0201001010128 4095  -2250        0                  -9      21132 7       4216 46    743    JD --       
            SQK  BG  DPZ IPU0201001010128 4083  -3250        0                  -4       4165 7       5102 62    673    JD --       
            GAC  NC  EHZ EPU2201001010128 4198  13 50        0                  -9      68 87 2       2322  6    127    JD --       
            GAXB NC  EHZ EPU2201001010128 4283  -8 50        0                 -30     136 68          158       126    J  --       
            GBG  NC  EHZ EPU2201001010128 4285  -5 50        0                   6     115 70           95        89    J  --       
            GCR  NC  EHZ EPU2201001010128 4226  -8 50        0                 -25     102 74          123        89    J  02       
            GDXB NC  HHZ IPD0201001010128 4101   5250        0                 -18      25126 2       4136 53    615    JD --       
            GPM  NC  EHZ EPD2201001010128 4289   8 50        0                  -6     116 70 2       3281 60    170    JD 02       
            GSG  NC  EHZ IPU1201001010128 4281  -5125        0                  24     103 74 2       3 62 -3    553    JD 02       
                                                                            71328695
            201001010128445938 4953122 4892  257     7106  1   135255 21516834  56 65     37    0 122 178  4      10     0       D  7 D 65 10         71033519D 65  10        2FNC01
            AL4  BG  DPZ EP 3201001010128 4535   3 20        0                   0      23133 5       3309 25    303    JD --       
            AL6  BG  DPZ EPU2201001010128 4575  -1 40        0                   0      49105 7       1234-31    550    JD --       
            CLV  BG  DPZ EP 3201001010128 4536   1 20        0                  -2      26129 5      15 56166    239    JD --       
            DRK  BG  DPZ EP 3201001010128 4550  -8 20        0                  -6      43111 5       5165 71     24    JD --       
            SB4  BG  DPZ IPD0201001010128 4522   1200        0                  -9      22134 5       5213 65    950    JD --       
            SQK  BG  DPZ IPU0201001010128 4511   0200        0                  -4       5166 5       7117 93    954    JD --       
            GDXB NC  HHZ IPU0201001010128 4521   0200        0                 -18      27129 0       5137 65    977    JD --       
                                                                            71033519
            """

        if isinstance( input, basestring ):
            istream = getQPDataSource( input, **kwargs )
        else:
            istream = input

        event_line  = False
        event_ctr   = 0

        for line in istream:

            ncsn_pars = line.split()
            # check type of line
            if len( line ) == 169:

                ## event line
                event_line = True
                pick_ctr   = 0

                # create event
                curr_id = ncsn_pars[0]
                ev = Event( ''.join( ('smi:NCSN/event/', curr_id) ) )
                self.eventParameters.event.append( ev )
                
                # create origin
                ori = Origin( ''.join( ('smi:NCSN/origin/', curr_id) ) )
                ori.add( ev, 'origin' )
                ori.latitude  = RealQuantity( float( float(line[16:18])+(float(line[19:23])/6000)  ) )
                ori.longitude = RealQuantity( float( (float(line[23:26])+(float(line[27:31])/6000))*-1  ) )
                ori.depth     = RealQuantity( float( float(line[31:36])/100 ) )
                 
                # get time components
                datetime_str = ncsn_pars[0]
                year   = datetime_str[0:4]
                month  = datetime_str[4:6]
                day    = datetime_str[6:8]
                hour   = datetime_str[8:10]
                mins   = datetime_str[10:12]
                sec    = (float ( datetime_str[12:16]) )/100
                if float(sec) >= 0:
                   ori.time = TimeQuantity( ( int(year), int(month), int(day),
                                           int(hour), int(mins), float(sec) ) )
                else:
                   ori.time = TimeQuantity( ( int(year), int(month), int(day),
                                           int(hour), int(mins), float('0.00') ) )

                # create magnitude
                mag = Magnitude( ''.join( ('smi:NCSN/magnitude/', curr_id) ) )
                # tmag=line[36:39].strip()
                if len(line[147:150].strip()) == 0:
                    raise 'no_magnitude_found'
                else:
                    tmag=float(line[147:150].strip())/100
                mag.mag = RealQuantity( float(tmag) )
                # TODO: correct magnitude type
                mag.type = 'mb'

                mag.add( ev, 'magnitude' )
                mag.setOriginAssociation( ori.publicID )
                ev.preferredMagnitudeID = mag.publicID
                if len(line[70:73].strip()) > 0:
                   mag = Magnitude()
                   mag.mag = RealQuantity( float(line[70:74].strip())/100 )
                   mag.type = 'Md'  #Duration magnitude
                   mag.add( ev, 'magnitude' )
                   mag.setOriginAssociation( ori.publicID )
                if len(line[123:126].strip()) > 0:
                   mag = Magnitude()
                   mag.mag = RealQuantity( float(line[123:126].strip())/100 )
                   mag.type = 'M'+line[122].strip()  #external magnitude computed by UCB
                   mag.add( ev, 'magnitude' )
                   mag.setOriginAssociation( ori.publicID )
                if len(line[130:133].strip()) > 0: 
                   mag = Magnitude()
                   mag.mag = RealQuantity( float(line[130:133].strip())/100 )
                   mag.type = 'M'+line[129].strip()  #Alternate amplitude magnitude
                   mag.add( ev, 'magnitude' )
                   mag.setOriginAssociation( ori.publicID )
                if len(line[147:150].strip()) > 0:
                   mag = Magnitude()
                   mag.mag = RealQuantity( float(line[147:150].strip())/100 )
                   mag.type = 'M'+line[146].strip()  #Prefered magnitude
                   mag.add( ev, 'magnitude' )
                   mag.setOriginAssociation( ori.publicID )
                if len(line[155:158].strip()) > 0:
                   mag = Magnitude()
                   mag.mag = RealQuantity( float(line[155:158].strip())/100 )
                   mag.type = 'M'+line[154].strip()  #Alternate coda mag
                   mag.add( ev, 'magnitude' )
                   mag.setOriginAssociation( ori.publicID )

                # set preferred origin
                ev.preferredOriginID    = ori.publicID
                
                event_ctr += 1

            elif len( line ) > 81 and len( line ) < 150:
                # if no event line defined: error
                if not event_line:
                    raise ValueError, "QPCatalog::importHypoInverse - phase line without event"
                ## phase line
                pick_ctr    = pick_ctr + 1

                # create pick
                curr_pickid = ''.join( ('smi:NCSN/event/', curr_id, '/pick/', str(pick_ctr)) )
                pick = Pick( curr_pickid )
                
                # get pick time
                year   = line[17:21]
                monnth  = line[21:22]
                day    = line[23:25]
                hour   = line[25:27]
                mins   = line[27:29]
                if len(line[29:34].strip()) < 1:
                    pick.time = TimeQuantity( ( int(year), int(month), int(day),
                                           int(hour), int(mins), float('0.00') ) )
                elif int( (float (line[29:34]))/100 ) == 60:
                    if int(mins) == 59 and int(hour) < 23:
                       pick.time = TimeQuantity(( int(year), int(month), int(day),
                                          int(int(hour)+1), int('0'), float('0.00') ) )
                    elif int(mins) == 59 and hour == 23:
                       pick.time = TimeQuantity(( int(year), int(month), int(day),
                                           int(hour), int(mins), float('59.9') ) )
                    else:
                       pick.time = TimeQuantity(( int(year), int(month), int(day),
                                           int(hour), int(int(mins)+1), float('0.00') ) )
                else:
                    sec    = (float (line[29:34]))/100
                    pick.time = TimeQuantity(( int(year), int(month), int(day),
                                           int(hour), int(mins), float(sec) ) )

                # get waveform id
                pick.waveformID = WaveformStreamID( ncsn_pars[1], ncsn_pars[0], ncsn_pars[2] )
                pick.add( ev, 'pick' )

                # create arrival
                arrv = Arrival()
                arrv.pickID   = curr_pickid
                arrv.phase    = Phase( line[14:15] )
                if len( line[74:79].strip() ) < 1:
                    arrv.distance = float (-10)
                else:
                    arrv.distance = (float( line[74:79] ))/10
                arrv.add( ori, 'arrival' )

            elif len( line ) < 79:
                # skip empty line
                continue
            else:
                print line
                print len( line )
                raise ValueError, "QPCatalog::importHypoInverse - format error in input stream"




    # -------------------------------------------------------------------------
                
    def cut( self, polygon=None, grid=None, geometry=None, **kwargs ):
        """
        polygon: tuple, list, numpy.array OR QPPolygon object
        grid:    QPGrid object
        geometry: Shapely aggregate geometry object
        
        kwargs:  min* and max* for lat, lon, depth, time, magnitude
                 default: cutting limits are included in result
                 if limit should be excluded, set
                 kwarg min*_excl=True (example: minmag_excl=True)

                 removeNaN=True: remove values that are set to 'NaN' when cutting (default: do not remove)
            
        strategy for cutting: event is deleted if one origin/magnitude meets criterion
            loop over event array (from end)
                (1) loop over origins:
                    - check if lat/lon falls into polygon: if not, break origin loop, delete event, continue event loop
                      keep event if it falls on polygon boundary
                      (Polygon: isInside(), QPPolygon: isInsideOrOnBoundary())
                    - if cut criterion met for lat/lon/depth/time, break origin loop, delete event, continue event loop
                (2) loop over magnitudes: if cut criterion met, break mag loop, delete event, continue event loop
        """
        cut_params = kwargs.keys()
        
        # check if polygon is given
        if polygon is not None:
            if isinstance( polygon, QPPolygon.QPPolygon ):
                poly_area = polygon
            else:
                poly_area = QPPolygon.QPPolygon( polygon )
        else:
            poly_area = None
            
        for curr_ev_idx in reversed( xrange( len( self.eventParameters.event ) ) ):
            
            cut_ev = False
            
            # need to go into origins?
            if    poly_area is not None \
               or grid is not None \
               or geometry is not None \
               or 'minlat' in cut_params or 'maxlat' in cut_params \
               or 'minlon' in cut_params or 'maxlon' in cut_params \
               or 'mindepth' in cut_params or 'maxdepth' in cut_params \
               or 'mintime' in cut_params or 'maxtime' in cut_params:
                
                for curr_ori in self.eventParameters.event[curr_ev_idx].origin:
                    
                    # check for polygon
                    if poly_area is not None:
                        if not poly_area.isInsideOrOnBoundary( float( curr_ori.longitude.value ),
                                                               float( curr_ori.latitude.value ) ):
                            cut_ev = True
                            break
                        
                    # check for grid
                    if grid is not None:
                        
                        # check if grid is correct object
                        try:
                            isinstance( grid, QPGrid.QPGrid )
                        except:
                            raise TypeError, 'QPCatalog::cut - grid has wrong type'
                        
                        if not grid.inGrid( float(curr_ori.latitude.value), 
                                            float(curr_ori.longitude.value),
                                            float(curr_ori.depth.value) ):
                            cut_ev = True
                            break    
                    
                    # check for geometry
                    if geometry is not None:
                        ev_point = shapely.geometry.Point( 
                            float( curr_ori.longitude.value ),
                            float( curr_ori.latitude.value ) ) 
                        if not ( geometry.contains( ev_point ) or 
                                 geometry.touches( ev_point ) ):
                            cut_ev = True
                            break

                    # cut latitude
                    if (    'removeNaN' in cut_params
                         and kwargs['removeNaN'] is True
                         and (    curr_ori.latitude.value is None
                               or numpy.isnan( curr_ori.latitude.value ) ) ):
                        cut_ev = True
                        break
                    else:
                        if 'minlat' in cut_params:
                            if 'minlat_excl' in cut_params and kwargs['minlat_excl'] is True:
                                if float(curr_ori.latitude.value) <= kwargs['minlat']:
                                    cut_ev = True
                                    break
                            else:
                                if float(curr_ori.latitude.value) < kwargs['minlat']:
                                    cut_ev = True
                                    break
                                    
                        if 'maxlat' in cut_params:
                            if 'maxlat_excl' in cut_params and kwargs['maxlat_excl'] is True:
                                if float(curr_ori.latitude.value) >= kwargs['maxlat']:
                                    cut_ev = True
                                    break
                            else:
                                if float(curr_ori.latitude.value) > kwargs['maxlat']:
                                    cut_ev = True
                                    break
                
                    # cut longitude
                    if (    'removeNaN' in cut_params
                         and kwargs['removeNaN'] is True
                         and (    curr_ori.longitude.value is None
                               or numpy.isnan( curr_ori.longitude.value ) ) ):
                        cut_ev = True
                        break
                    else:
                        if 'minlon' in cut_params:
                            if 'minlon_excl' in cut_params and kwargs['minlon_excl'] is True:
                                if float(curr_ori.longitude.value) <= kwargs['minlon']:
                                    cut_ev = True
                                    break
                            else:
                                if float(curr_ori.longitude.value) < kwargs['minlon']:
                                    cut_ev = True
                                    break
                                    
                        if 'maxlon' in cut_params:
                            if 'maxlon_excl' in cut_params and kwargs['maxlon_excl'] is True:
                                if float(curr_ori.longitude.value) >= kwargs['maxlon']:
                                    cut_ev = True
                                    break
                            else:
                                if float(curr_ori.longitude.value) > kwargs['maxlon']:
                                    cut_ev = True
                                    break
                    
                    # cut depth
                    if (    'removeNaN' in cut_params
                         and kwargs['removeNaN'] is True
                         and (    curr_ori.depth.value is None
                               or numpy.isnan( curr_ori.depth.value ) ) ):
                        cut_ev = True
                        break
                    else:
                        if 'mindepth' in cut_params:
                            if 'mindepth_excl' in cut_params and kwargs['mindepth_excl'] is True:
                                if float(curr_ori.depth.value) <= kwargs['mindepth']:
                                    cut_ev = True
                                    break
                            else:
                                if float(curr_ori.depth.value) < kwargs['mindepth']:
                                    cut_ev = True
                                    break
                                
                        if 'maxdepth' in cut_params:
                            if 'maxdepth_excl' in cut_params and kwargs['maxdepth_excl'] is True:
                                if float(curr_ori.depth.value) >= kwargs['maxdepth']:
                                    cut_ev = True
                                    break
                            else:
                                if float(curr_ori.depth.value) > kwargs['maxdepth']:
                                    cut_ev = True
                                    break
                    
                    # cut time
                    if 'mintime' in cut_params:
                        mintime = QPDateTime( kwargs['mintime'] )
                        if 'mintime_excl' in cut_params and kwargs['mintime_excl'] is True:
                            if curr_ori.time.value <= mintime:
                                cut_ev = True
                                break
                        else:
                            if curr_ori.time.value < mintime:
                                cut_ev = True
                                break

                    if 'maxtime' in cut_params:
                        maxtime = QPDateTime( kwargs['maxtime'] )
                        if 'maxtime_excl' in cut_params and kwargs['maxtime_excl'] is True:
                            if curr_ori.time.value >= maxtime:
                                cut_ev = True
                                break
                        else:
                            if curr_ori.time.value > maxtime:
                                cut_ev = True
                                break
                
                # origin loop finished    
                if cut_ev is True:
                    self.eventParameters.event.pop( curr_ev_idx )
                    continue     
                   
            # need to go into magnitudes?
            if 'minmag' in cut_params or 'maxmag' in cut_params:

                for curr_mag in self.eventParameters.event[curr_ev_idx].magnitude:

                    if (    'removeNaN' in cut_params
                         and kwargs['removeNaN'] is True
                         and (    curr_mag.mag.value is None
                               or numpy.isnan( curr_mag.mag.value ) ) ):
                        cut_ev = True
                        break
                    else:
                        if 'minmag' in cut_params:
                            if 'minmag_excl' in cut_params and kwargs['minmag_excl'] is True:
                                if float(curr_mag.mag.value) <= kwargs['minmag']:
                                    cut_ev = True
                                    break
                            else:
                                if float(curr_mag.mag.value) < kwargs['minmag']:
                                    cut_ev = True
                                    break

                        if 'maxmag' in cut_params:
                            if 'maxmag_excl' in cut_params and kwargs['maxmag_excl'] is True:
                                if float(curr_mag.mag.value) >= kwargs['maxmag']:
                                    cut_ev = True
                                    break
                            else:
                                if float(curr_mag.mag.value) > kwargs['maxmag']:
                                    cut_ev = True
                                    break    
                
                # mag loop finished
                if cut_ev is True:
                    self.eventParameters.event.pop( curr_ev_idx )
                    continue     

    def rebin( self, binsize=0.1, allorigins=True ):
        
        # loop over events
        for curr_ev in self.eventParameters.event:
            
            # if allorgins = False, rebin only preferred origins, otherwise all
            if allorigins is False:
                curr_ori = curr_ev.getPreferredOrigin()
                
                # rebin all magnitudes
                for curr_mag_idx in curr_ev.getMagnitudesIdx( curr_ori ):
                    curr_mag_value = curr_ev.magnitude[curr_mag_idx].mag.value
                    curr_ev.magnitude[curr_mag_idx].mag.value = self.mag_rebin( curr_mag_value, binsize )
                
            else:
                for curr_ori in curr_ev.origin:
                    
                    # rebin all magnitudes
                    for curr_mag_idx in curr_ev.getMagnitudesIdx( curr_ori ):
                        curr_mag_value = curr_ev.magnitude[curr_mag_idx].mag.value
                        curr_ev.magnitude[curr_mag_idx].mag.value = self.mag_rebin( curr_mag_value, binsize )

    def size( self ):
        return len( self.eventParameters.event )
    
    def timeSpan( self ):
        """
        Compute time span of events (use preferred origins).
        
        Returns a triple of time difference (in years), start time, end time
        """
        time_start = None
        time_end = None
        
        for curr_ev in self.eventParameters.event:
            curr_ori = curr_ev.getPreferredOrigin()
            curr_time = curr_ori.time.value
            
            if time_start is None:
                time_start = curr_time
                time_end = curr_time
            else:
                if curr_time < time_start:
                    time_start = curr_time
                elif curr_time > time_end:
                    time_end = curr_time
                    
        time_diff = diffQPDateTime( time_end, time_start )
        return ( time_diff.days / 365.25, time_start.datetime, time_end.datetime )
        
    def mag_rebin( self, curr_mag, binsize=0.1 ):
        """
        mag_rebin
        """
        return binsize * round( float(curr_mag) / binsize ) 
        
    def getFmd( self, allorigins=False, **kwargs ):
        """
        getFmd
        """
        # self.rebin( binsize )
        self.frequencyMagnitudeDistribution = \
            qpfmd.FrequencyMagnitudeDistribution( 
                self.eventParameters, **kwargs )
        return self.frequencyMagnitudeDistribution
    
    def getCumulativeDistribution( self ):
        """
        getCumulativeDistribution
        """
        self.cumulativeDistribution = cumuldist.CumulativeDistribution(
            self.eventParameters )
        return self.cumulativeDistribution 

    def toCompact( self ):

        compact = QPCatalogCompact.QPCatalogCompact()
        compact.update( self.catalog )

        return compact

    def fromCompact( self, compact ):
        pass

def main():
    pass
    
if __name__ == '__main__':
    main()
  

