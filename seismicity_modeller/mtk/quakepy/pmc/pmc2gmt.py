#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# quakepy/pmc/pmc2gmt.py
# $Id: pmc2gmt.py 242 2009-10-20 04:26:16Z fab $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2007-2009 by Fabian Euchner and Danijel Schorlemmer     #
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

pmc2gmt.py [options]

extract PMC results from QPGrid XML and write to GMT column format

reads XML data from stdin (default), or from file (switch -f)
uses lxml.etree default parser which creates a tree in memory
if switch -g is specified, use lxml.etree target parser (needs less memory, but slower)

commandline options:

    -f FILE    filename of XML file to process. can be g'zipped, with mandatory extension .gz
    -p VALUE   extract magnitude of completeness for given probability
    -m VALUE   extract probability for given magnitude of completeness
    -d VALUE   extract data for given depth layer, in km (default: first)
    -s         extract station information
    -g         use lxml.etree target parser (needs less memory, but slower)
    -z         input XML is gnuzipped (only required if read from stdin) 
    -n         input XML has no namespace prefix

input:  QPGrid XML, from file or from stdin
output: GMT format, to standard output

requires package lxml (version >= 2.0)
"""

__version__  = '$Id: pmc2gmt.py 242 2009-10-20 04:26:16Z fab $'
__revision__ = '$Revision: 242 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import os
import sys
import getopt
import gzip

from lxml import etree

MM_MAGNITUDE, MM_PROBABILITY, MM_STATION = range(1, 4)

# namespace of QPGrid elements
QP_NS = 'http://quakepy.org/xmlns/quakepy/1.0'

# format string for longitude and latitude
LONLAT_FMT = '%8.4f'

# format string for magnitudes
MAG_FMT = '%5.2f'

# format string for probabilities
PROB_FMT = '%e'

# format string for station line
STATION_FMT = '%14.8f%14.8f% 8s% 8s% 8s'

class CellTarget( object ):
    """
    define target parser for <cell> elements of QPGrid XML
    """

    def __init__( self, **kwargs ):

        self.output  = []
        
        self.ns      = kwargs['ns']
        self.mapMode = int( kwargs['mapMode'] )
        
        self.targetMagnitude   = float( kwargs['targetMagnitude'] )
        self.targetProbability = float( kwargs['targetProbability'] )
        
    def start( self, tag, attrib ):

        self.is_cell        = False
         
        self.is_mp          = False
        self.is_probability = False

        if self.mapMode != MM_STATION:
            
            if tag == '{%s}cell' % self.ns:
                
                self.is_cell = True

                self.curr_mag  = None
                self.curr_prob = None
            
                self.output = [ LONLAT_FMT % float( attrib['lon'] ),
                                LONLAT_FMT % float( attrib['lat'] ) ]
                
            elif tag == '{%s}mp' % self.ns:

                self.is_mp = True
                self.pmcdata_mp = {}

                self.curr_prob = attrib['probability']
                
            elif tag == '{%s}probability' % self.ns:

                self.is_probability = True
                self.pmcdata_prob = {}

                self.curr_mag  = attrib['magnitude']

        elif tag == '{%s}PMCStation' % self.ns:

            # if location code is not present in XML, output will be 'None' in last column
            print STATION_FMT % ( float( attrib['longitude'] ),
                                  float( attrib['latitude'] ),
                                  attrib['stationCode'],
                                  attrib['networkCode'],
                                  attrib['locationCode'] )

            
    def end( self, tag ):

        if self.mapMode != MM_STATION:
            
            if tag == '{%s}cell' % self.ns:

                self.is_cell = False

                ## magnitude map, target probability given
                if self.mapMode == MM_MAGNITUDE:

                    # check if it is correct probability
                    mp_found = False

                    for prob, mp in self.pmcdata_mp.iteritems():

                        if floatEqual( float( prob ), self.targetProbability ):

                            # found desired value in cell, output mp value
                            self.output.append( MAG_FMT % float( mp ) )
                            
                            mp_found = True
                            break

                    if mp_found is not True:

                        # Correct MP not found, compute from probability list

                        # Search for target probability
                        probList = []
                        for mp, prob in self.pmcdata_prob.iteritems():
                            
                            if float( prob ) > self.targetProbability:
                                probList.append( float( mp ) )

                        if len( probList ) > 0:
                            self.output.append( MAG_FMT % min( probList ) )
                        else:
                            self.output.append( 'nan' )
                            
                ## probability map, target magnitude given
                elif self.mapMode == MM_PROBABILITY:

                    # check if it is correct magnitude
                    prob_found = False

                    for mag, prob in self.pmcdata_prob.iteritems():

                        if floatEqual( float( mag ), self.targetMagnitude ):

                            # found desired value in cell, output prob value
                            self.output.append( PROB_FMT % float( prob ) )
                            
                            prob_found = True
                            break

                    if prob_found is not True:
                        self.output.append( 'nan' )
                    
                # output lat, lon and value for current cell
                print ' '.join( self.output )

            elif tag == '{%s}mp' % self.ns:

                self.is_mp = False
                
            elif tag == '{%s}probability' % self.ns:

                self.is_probability = False


    def data( self, data ):

        # fire only on <mp> or <probability> element contents
        if self.mapMode != MM_STATION:

            # build up dictionary of values
            if self.is_mp is True:
                self.pmcdata_mp[self.curr_prob] = data

            elif self.is_probability is True:
                self.pmcdata_prob[self.curr_mag] = data

    def close( self ):
        return "parser closed"


def parse_tree( input, **kwargs ):

    qpns    = kwargs['ns']
    mapMode = int( kwargs['mapMode'] )

    targetMagnitude   = float( kwargs['targetMagnitude'] )
    targetProbability = float( kwargs['targetProbability'] )
    targetDepthLayer  = kwargs['targetDepthLayer']

    tree = etree.parse( input )

    # Extract station list
    if mapMode == MM_STATION:

        if qpns == '':
            search_pattern = 'PMCStations/PMCStation'
        else:
            search_pattern = '{%s}PMCStations/{%s}PMCStation' % ( qpns, qpns )
            
        stations = tree.findall( search_pattern )

        for curr_station in stations:

            lat = curr_station.get( 'latitude' )
            lon = curr_station.get( 'longitude' )
            stationCode  = curr_station.get( 'stationCode' )
            networkCode  = curr_station.get( 'networkCode' )
            locationCode = curr_station.get( 'locationCode' )

            # if location code is not present in XML, output will be 'None' in last column
            print STATION_FMT % ( float( lon ), float( lat ), stationCode, networkCode, locationCode )

    # Extract either detection probability or completeness magnitude
    else:

        # find right depthLayer
        if qpns == '':
            search_pattern = 'grid/depthLayer'
        else:
            search_pattern = '{%s}grid/{%s}depthLayer' % ( qpns, qpns )

        depth_layers = tree.findall( search_pattern )

        # loop over depthLayers
        depthLayerFound = False
        for dl in depth_layers:

            if depthLayerFound is True:
                sys.exit()

            # get current depth layer 'at' value
            at = dl.get( 'at' )

            if targetDepthLayer is None:
                depthLayerFound = True
            elif floatEqual( float( at ), float( targetDepthLayer ) ):
                depthLayerFound = True
            else:
                continue

            # loop over cells
            if qpns == '':
                search_pattern = 'cell'
            else:
                search_pattern = '{%s}cell' % ( qpns )

            cell = dl.findall( search_pattern )

            # loop over lon/lat cells
            for c in cell:

                # get lat, lon
                lat = c.get( 'lat' )
                lon = c.get( 'lon' )

                output = [ LONLAT_FMT % float( lon ),
                           LONLAT_FMT % float( lat ) ]
                           
                ## magnitude map, target probability given
                if mapMode == MM_MAGNITUDE:

                    # get all mp bins
                    if qpns == '':
                        search_pattern = 'PMCData/mp'
                    else:
                        search_pattern = '{%s}PMCData/{%s}mp' % ( qpns, qpns )
                
                    mp = c.findall( search_pattern )

                    # check if it is correct probability
                    mp_found = False

                    for m in mp:
                        prob = m.get( 'probability' )
                        
                        if floatEqual( float( prob ), targetProbability ):

                            # found desired value in cell, output mp value
                            output.append( MAG_FMT % float( m.text ) )

                            mp_found = True
                            break

                    if mp_found is not True:
                        
                        # Correct MP not found, compute from probability list
                        if qpns == '':
                            search_pattern = 'PMCData/probability'
                        else:
                            search_pattern = '{%s}PMCData/{%s}probability' % ( qpns, qpns )
                            
                        prob = c.findall( search_pattern )

                        # Search for target probability
                        probList = []
                        for p in prob:

                            if float( p.text ) > targetProbability:
                                probList.append( float( p.get( 'magnitude' ) ) )

                        if len( probList ) > 0:
                            output.append( MAG_FMT % min( probList ) )
                        else:
                            output.append( 'nan' )

                ## probability map, target magnitude given
                else:

                    # get all probability bins
                    if qpns == '':
                        search_pattern = 'PMCData/probability'
                    else:
                        search_pattern = '{%s}PMCData/{%s}probability' % ( qpns, qpns )
                        
                    prob = c.findall( search_pattern )

                    # check if it is correct magnitude
                    prob_found = False
                    
                    for p in prob:
                        mag = p.get( 'magnitude' )

                        if floatEqual( float( mag ), targetMagnitude ):

                            # found desired value in cell, output prob value
                            output.append( PROB_FMT % float( p.text ) )

                            prob_found = True
                            break
                        
                    if prob_found is not True:
                        output.append( 'nan' )
                    
                # output lat, lon and value for current cell
                print ' '.join( output )


def main():

    scriptname = os.path.basename( os.path.abspath( sys.argv[0] ) )
    
    # default mode: output values for completeness map
    mapMode = MM_MAGNITUDE

    # default probability level: 0.999
    targetProbability = 0.999

    # default completeness magnitude: 2.4
    targetMagnitude = 2.4

    # default: no depth layer given (use first in XML document)
    targetDepthLayer = None

    # default is input from stdin
    stdinMode = True

    # compression mode, only required for input from stdin
    zipMode = False

    # target parsing mode (default is tree parsing)
    targetMode = False
    
    # inputfile 
    inputfile = None
    
    # default namespace: QuakePy
    qpns = QP_NS

    # Read commandline arguments
    cmdParams = sys.argv[1:]
    opts, args = getopt.gnu_getopt( cmdParams,
                                    'f:p:m:d:sgznh',
                                    [ 'inputfile=',
                                      'probability=',
                                      'magnitude=',
                                      'depthlayer=',
                                      'station',
                                      'target',
                                      'compressed', 
                                      'nonamespace',
                                      'help' ] )

    for option, parameter in opts:

        if option == '-f' or option == '--inputfile':
            inputfile = parameter
            stdinMode = False
            
        if option == '-p' or option == '--probability':
            targetProbability = float( parameter )
            mapMode = MM_MAGNITUDE

        if option == '-m' or option == '--magnitude':
            targetMagnitude = float( parameter )
            mapMode = MM_PROBABILITY

        if option == '-d' or option == '--depthlayer':
            targetDepthLayer = float( parameter )

        if option == '-s' or option == '--station':
            mapMode = MM_STATION

        if option == '-g' or option == '--target':
            targetMode = True
            
        if option == '-z' or option == '--compressed':
            zipMode = True
            
        if option == '-n' or option == '--nonamespace':
            qpns = ''

        if option == '-h' or option == '--help':
            PrintHelp()
            sys.exit()

    if stdinMode is True:

        # get input data from standard input
        if zipMode == True:
            input = gzip.GzipFile( fileobj = sys.stdin )
        else:
            input = sys.stdin
        
    else:

        # get input data from file, raise exception if file does not exist
        if os.path.isfile( inputfile ):
            input = inputfile
        else:
            error_str = "%s - input file %s does not exist" % ( scriptname, inputfile )
            raise IOError( error_str )

    if targetMode is True:
        
        parser_target  = CellTarget( mapMode=mapMode,
                                     targetMagnitude=targetMagnitude,
                                     targetProbability=targetProbability,
                                     ns=qpns )

        parser         = etree.XMLParser( target=parser_target )
        results        = etree.parse( input, parser )

    else:

        parse_tree( input,
                    mapMode=mapMode,
                    targetMagnitude=targetMagnitude,
                    targetProbability=targetProbability,
                    targetDepthLayer=targetDepthLayer,
                    ns=qpns )


def PrintHelp():
    print 'converts QPGrid XML file to GMT input'
    print 'Usage: pmc2gmt.py [OPTIONS]'
    print '  Options'
    print '   -f FILE, --inputfile=<filename>    Input grid file name'
    print '   -p VALUE, --probability=<value>    Target probability for completeness map (default: 0.999)'
    print '   -m VALUE, --magnitude=<value>      Target magnitude for probability map'
    print '   -d VALUE, --depthlayer=<value>     Depth layer (in km) for which data is extracted (default: first)'
    print '   -s, --station                      Station locations for station map'
    print '   -g, --target                       Use target parser (uses less memory, but slower)'
    print '   -z, --compressed                   Input XML is gnuzipped (only required if read from stdin)'
    print '   -n, --nonamespace                  XML input has no namespace prefix'
    print '   -h, --help                         print this information'

    
def floatEqual( f1, f2, epsilon = 1e-12 ):
    """
    checks if two floating point values can be considered equal
    returns True if difference of the two values is smaller or equal to epsilon
    """
    if abs( f1 - f2 ) > epsilon:
        return False
    else:
        return True

    
if __name__ == "__main__":
    main()