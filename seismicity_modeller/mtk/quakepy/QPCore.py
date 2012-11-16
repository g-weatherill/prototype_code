# -*- coding: utf-8 -*-
#
# quakepy/QPCore.py
# $Id: QPCore.py 268 2010-09-21 12:33:15Z fab $
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
"""

__version__  = '$Id: QPCore.py 268 2010-09-21 12:33:15Z fab $'
__revision__ = '$Revision: 268 $'
__author__   = "Fabian Euchner <fabian@fabian-euchner.de>, Danijel Schorlemmer <ds@usc.edu>"
__license__  = "GPL"

import datetime
import math
import numpy
import urllib
import gzip, bz2
import types

import pyRXP

from mx.DateTime import DateTime, DateTimeType, utc
from mx.DateTime import cmp as mxdatetimecmp

# internal includes

from QPElement import *
from QPUtils import *
from QPDateTime import *


POS_TAGNAME, POS_ATTRS, POS_CHILDREN = range(3)


class QPObjectList( list ):
    pass


class QPObject( object ):
    """
    QuakePy: QPObject
    """

    # this is the difference below which two floating point values
    # are considered equal
    floatCmpEpsilon = 1.0e-9
    
    # this is the difference (in seconds) below which two [mx|QP]DateTime values
    # are considered equal
    dateTimeCmpEpsilon = 1.0e-9

    # number of decimal places for seconds
    secondsDigits = 6 # 6
    
    def __init__( self, **kwargs ):
        """
        each class that is derived from QPObject can redefine/append to the elements list
        element = QPElement( varname, xmlname, xmltype, type, vartype, parentname=None, parenttype=None )
        xmltype can be from ( element, attribute, cdata )
        vartype can be from ( basic, enum, complex, multiple )

        kwargs: parentAxis
                elementName
        """

        self.elements     = QPElementList()
        self.childXMLTree = []

        self.elementAxis = ''

        if ( 'parentAxis' in kwargs.keys() ) and ( kwargs['parentAxis'] is not None ):
            if 'elementName' in kwargs.keys():
                self.setElementAxis( kwargs['parentAxis'], kwargs['elementName'] )
            else:
                self.setElementAxis( kwargs['parentAxis'] )

    # ------------------------------------------------------------------------

    def __eq__( self, T ):
        """
        compare two QPObjects for equality
        compare all elements from elements list separately
        
        compare basic types: unicode, string, int with '==' operator
        compare basic type float with 'epsilon' range
        compare type mxDateTime / QPDateTime with 'epsilon' range

        do not compare publicID attributes and elements, return 'True' without check
        """
        
        ## XML attributes
        
        # compare attributes: can only be basic or enum
        for xmlname, varname, pytype in self._getXMLAttributeNames():
            if hasattr( self, varname ) and self.__dict__[varname] is not None:

                # no comparison for 'publicID'
                if varname == 'publicID':
                    return True
                
                if self.compareEqualBasicType( T, varname, pytype ) is False:
                    return False

        ## XML elements

        # basic type
        for xmlname, varname, pytype in self._getXMLElementNames( 'basic' ):
            if hasattr( self, varname ) and self.__dict__[varname] is not None:

                # no comparison for 'publicID'
                if varname == 'publicID':
                    return True
                
                if self.compareEqualBasicType( T, varname, pytype ) is False:
                    return False

        # enums
        for xmlname, varname, pytype in self._getXMLElementNames( 'enum' ):
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
                if self.compareEqualBasicType( T, varname, pytype ) is False:
                    return False

        # complex types
        for xmlname, varname, pytype in self._getXMLElementNames( 'complex' ):
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
              
                if hasattr( T, varname ) and T.__dict__[varname] is not None:
                    if self.__dict__[varname].__eq__( T.__dict__[varname] ) is False:
                        return False
                else:
                    return False
                
        # multiple elements
        for xmlname, varname, pytype in self._getXMLElementNames( 'multiple' ):
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
              
                if hasattr( T, varname ) and T.__dict__[varname] is not None:
                    
                    varlen = len( self.__dict__[varname] )
                    if len( T.__dict__[varname] ) != varlen:
                        return False
                    
                    for idx in xrange( varlen ):
                        if self.__dict__[varname][idx].__eq__( T.__dict__[varname][idx] ) is False:
                            return False
                else:
                    return False

        ## CDATA
        varname, pytype = self._getXMLCDATAName()
        if varname is not None:
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
                if self.compareEqualBasicType( T, varname, pytype ) is False:
                    return False
                
        # if no unequal comparison so far, return True
        return True


    def __ne__( self, T ):
      
        if self.__eq__( T ) is True:
            return False
        else:
            return True


    def compareEqualBasicType( self, T, varname, pytype, verbose=False ):
        """
        compare basic attributes 'varname' of self and another instance
        """
        
        # check if attribute 'varname' is present in both compared instances
        if (     hasattr( self, varname ) and self.__dict__[varname] is not None
             and hasattr( T, varname ) and T.__dict__[varname] is not None ):
            
            cmpEqual = True
            
            if isinstance( self.__dict__[varname], float ):
                if not floatEqual( self.__dict__[varname], T.__dict__[varname], self.floatCmpEpsilon ):
                    cmpEqual = False
            
            elif isinstance( self.__dict__[varname], DateTimeType ):
                
                # compare instances of mx.DateTime.DateTimeType
                # use special 'cmp' function because we need to use epsilon
                # cmp() from mxDateTime returns 0 for equal, -1 for smaller, 1 for greater
                if mxdatetimecmp( self.__dict__[varname], T.__dict__[varname], self.dateTimeCmpEpsilon ) != 0:
                    cmpEqual = False

            elif isinstance( self.__dict__[varname], QPDateTime ):
                
                # do not compare instances of QPDateTime using the regular 'cmp' function / operator
                # since we need to use 'epsilon'
                if cmpQPDateTime( self.__dict__[varname], T.__dict__[varname], self.dateTimeCmpEpsilon ) != 0:
                    cmpEqual = False
            else:
              
                if not ( self.__dict__[varname] == T.__dict__[varname] ):
                    cmpEqual = False
                
            if cmpEqual is False:
                if verbose is True:
                    print "%s: comparison not equal for attribute %s of %s - self %s, other %s" % ( self.__class__.__name__, 
                                                                                                    varname, 
                                                                                                    pytype, 
                                                                                                    self.__dict__[varname], 
                                                                                                    T.__dict__[varname] )
                    if isinstance( self.__dict__[varname], float ):
                        print " --> difference: %s, epsilon: %s" % ( self.__dict__[varname] - T.__dict__[varname],
                                                                     self.floatCmpEpsilon )
                    elif isinstance( self.__dict__[varname], DateTimeType ):
                        TimeDiff = self.__dict__[varname] - T.__dict__[varname]
                        print " --> difference: %s, epsilon: %s" % ( TimeDiff.seconds, self.dateTimeCmpEpsilon )
                        
                    elif isinstance( self.__dict__[varname], QPDateTime ):
                        TimeDiff = diffQPDateTime( self.__dict__[varname], T.__dict__[varname] )
                        print " --> difference: %s, epsilon: %s" % ( TimeDiff.seconds, self.dateTimeCmpEpsilon )
                        
                return False
            else:
                # if comparison was equal, return True
                return True
            
        else:
            return False

    # ------------------------------------------------------------------------
    
    def _getXMLAttributeNames( self ):
        """
        get object attributes with 'attribute' xmltype from elements list 
        """
        attr = []
        for element in self.elements:
            if element.xmltype == 'attribute':
                attr.append( [ element.xmlname, element.varname, element.pytype ] )
        return attr


    def _getXMLElementNames( self, vartype ):
        """
        get object attributes with 'element' xmltype from elements list
        for specific vartype
        
        vartype can be from ( basic, enum, complex, multiple )
        """
        elem = []
        for element in self.elements:
            if ( ( element.xmltype == 'element' ) and ( element.vartype == vartype ) ):
                elem.append( [ element.xmlname, element.varname, element.pytype ] )
        return elem


    def _getXMLCDATAName( self ):
        """
        get object attribute with 'cdata' xmltype from elements list
        
        there should not be more than one CDATA element in self.elements
        function returns the first CDATA element
        """
        for element in self.elements:
            if element.xmltype == 'cdata':
                return [ element.varname, element.pytype ]
        return [ None, None ]


    def _getXMLExtensionElements( self, elementList ):
        """
        get extension elements from a given element list
        these elements must have parentaxis and parenttype attributes
        """
        elem = []
        for element in elementList:
            if ( ( element.parentaxis is not None ) and ( element.parenttype is not None ) ):
                elem.append( element )
        return elem


    def _initMultipleElements( self ):
        for xmlname, varname, pytype in self._getXMLElementNames( 'multiple' ):
            self.__dict__[varname] = QPObjectList()


    def _getXMLSerializationString( self, attribute ):
        """
        get string for XML serialization of an attribute
        for all basic types except mxDateTime and QPDateTime use built-in 'unicode' method
        for mxDateTime use mxDateTime2ISO() from QPUtils
        for QPDateTime use method toISO()
        """
        
        # print "checking attribute %s of %s" % ( attribute, type(attribute) )
        # print "return unicode %s" % unicode( attribute )
        
        if isinstance( attribute, DateTimeType ):
            return unicode( mxDateTime2ISO( attribute, secondsdigits = self.secondsDigits ) )
        elif isinstance( attribute, QPDateTime ):
            return unicode( attribute.toISO( secondsdigits = self.secondsDigits ) )
        else:
            return unicode( attribute )

    # ------------------------------------------------------------------------

    def fromXML( self, tree, additionalElements=None ):
        """
        populate class attributes as defined in 'elements' list of the object
        from pyRXP XML tuple tree

        additionalElements is another QPElementList which is checked
        for attributes that are added to standard class layout
        """

        foundElements = []

        ## check if there are additionalElements which have to be added to
        ## this object
        if additionalElements is not None:

            for check_element in self._getXMLExtensionElements( additionalElements ):
                
                # check if elementAxis and type match
                if (     ( self.elementAxis == check_element.parentaxis )
                     and isinstance( self, check_element.parenttype ) ):

                    # append element to current elements list
                    self.elements.append( check_element )

        ## XML attributes
        if tree[POS_ATTRS] is not None:
            attr_dict = tree[POS_ATTRS]

            # loop over possible attributes
            for xmlname, varname, pytype in self._getXMLAttributeNames():
                if xmlname in attr_dict.keys():
                    self.__dict__[varname] = pytype( attr_dict[xmlname] )
                
        ## XML elements
        if tree[POS_CHILDREN] is not None:
            for child in tree[POS_CHILDREN]:

                ## get CDATA
                ## if child contains only whitespace, continue immediately
                if not isinstance( child, tuple ):
                     
                    if len( child.strip() ) > 0:
                        varname, pytype = self._getXMLCDATAName()
                        self.__dict__[varname] = pytype( child )

                    continue


                elementFound = False
                
                # basic type
                for xmlname, varname, pytype in self._getXMLElementNames( 'basic' ):

                    if child[POS_TAGNAME] == xmlname and len( child[POS_CHILDREN] ) > 0:
                        self.__dict__[varname] = pytype( child[POS_CHILDREN].pop() )
                        addUnique( foundElements, [ xmlname ] )
                        elementFound = True
                        break
                    
                if elementFound is True:
                    continue
                
                # enums
                for xmlname, varname, pytype in self._getXMLElementNames( 'enum' ):

                    if child[POS_TAGNAME] == xmlname and len( child[POS_CHILDREN] ) > 0:
                        self.__dict__[varname] = pytype( child[POS_CHILDREN].pop() )
                        addUnique( foundElements, [ xmlname ] )
                        elementFound = True
                        break

                if elementFound is True:
                    continue
                
                # complex types
                for xmlname, varname, pytype in self._getXMLElementNames( 'complex' ):

                    if child[POS_TAGNAME] == xmlname:
                        self.__dict__[varname] = pytype( parentAxis = self.elementAxis,
                                                         elementName = varname )
                        self.__dict__[varname].fromXML( child, additionalElements )
                        addUnique( foundElements, [ xmlname ] )

                if elementFound is True:
                    continue
                
                # multiple elements
                for xmlname, varname, pytype in self._getXMLElementNames( 'multiple' ):

                    if child[POS_TAGNAME] == xmlname:
                        tmp = pytype( parentAxis = self.elementAxis,
                                      elementName = varname )
                        tmp.fromXML( child, additionalElements )

                        tmp.add( self, varname )
                        addUnique( foundElements, [ xmlname ] )

                if elementFound is True:
                    continue

                ## element has not been found:
                ## append subtree to childXMLTree if not already processed
                if not child[POS_TAGNAME].strip() in foundElements:
                    self.childXMLTree.append( child )

        return True


    def toXML( self, tagname, stream ):
        """
        create XML representation for elements given in 'element' tuple
        and write XML representation to stream
        """

        stream.writelines( [ '<', tagname ] )
        
        ## XML attributes
        
        # loop over possible attributes
        for xmlname, varname, pytype in self._getXMLAttributeNames():
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
                stream.writelines( [ ' ', xmlname, '="', self._getXMLSerializationString(self.__dict__[varname]), '"' ] )

        stream.write( '>' )
        
        ## XML elements

        # basic type
        for xmlname, varname, pytype in self._getXMLElementNames( 'basic' ):
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
                stream.writelines( [ '<', xmlname, '>', self._getXMLSerializationString(self.__dict__[varname]), '</', xmlname, '>' ] )
                
        # enums
        for xmlname, varname, pytype in self._getXMLElementNames( 'enum' ):
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
                stream.writelines( [ '<', xmlname, '>', self._getXMLSerializationString(self.__dict__[varname]), '</', xmlname, '>' ] )
                
        # complex types
        for xmlname, varname, pytype in self._getXMLElementNames( 'complex' ):
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
                self.__dict__[varname].toXML( xmlname, stream )
                
        # multiple elements
        for xmlname, varname, pytype in self._getXMLElementNames( 'multiple' ):
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
                for tmp in self.__dict__[varname]:
                    tmp.toXML( xmlname, stream )

        # add non-standard elements from self.childXMLTree
        for curr_extension_node in self.childXMLTree:
            pyrxpTupleTree2XML( curr_extension_node, stream )

        ## add CDATA
        ## note: there should by only one CDATA element in self.elements
        ## we use the first element that is flagged as 'cdata'
        varname, pytype = self._getXMLCDATAName()
        if varname is not None:
            if hasattr( self, varname ) and self.__dict__[varname] is not None:
                stream.write( self._getXMLSerializationString(self.__dict__[varname]) )
                
        stream.writelines( [ '</', tagname, '>' ] )
        return True

    # ------------------------------------------------------------------------

    def add( self, parentObject, parentVariableName ):
        """
        add object (self) to list 'parentVariableName' in parent class
        (make it a child element of parentObject)
        """

        # check if parentObject type has attribute parentVariableName
        if hasattr( parentObject, parentVariableName ) and parentObject.__dict__[parentVariableName] is not None:

            # check if parentObject.parentVariableName is defined as list of
            # objects of the same type as self
            for curr_element in parentObject.elements:
                
                if (     ( parentVariableName == curr_element.varname )
                     and ( isinstance( self, curr_element.pytype ) ) ):
            
                        parentObject.__dict__[parentVariableName].append( self )
                        self.setElementAxis( parentObject.elementAxis, parentVariableName )
                        return True

            # loop has ended without success
            error_str = "QPObject::add - cannot append self to attribute %s of parent class %s" % \
                        ( parentVariableName, parentObject.__class__.__name__ )
            raise TypeError, error_str
                
        else:
            error_str = "QPObject::add - parent class %s has no attribute %s" % \
                        ( parentObject.__class__.__name__, parentVariableName )
            raise TypeError, error_str


    def addObject( self, obj, obj_elements ):
        """
        obj has to have obj.elements for own contents
        obj_elements has to be QPElement with element information for added obj
        """
        self.elements.append( obj_elements )

        # add object to instance
        self.__dict__[obj_elements.varname] = obj


    def setElementAxis( self, parentAxis, elementName=None ):

        if parentAxis is not None:

            if ( elementName is not None ) and ( elementName != '' ):
                self.elementAxis = '/'.join( ( parentAxis, elementName ) )
            else:
                self.elementAxis = parentAxis

# ----------------------------------------------------------------------------

class QPPublicObject( QPObject ):
    """
    QuakePy: QPPublicObject
    """

    # public object counter, incremented at each object creation
    publicObjectCtr = 0L

    # style of publicID
    # Note: colons are not allowed in smi URIs, so change ':' to '.' in time component
    # of timestamp
    #
    #  'full':    'smi:local/generic/2008-07-26T15:00:00'
    #  'short':   'myname' or '2008-07-26T15:00:00'
    #  'numeric': '152763'

    publicIDStyle = 'full'

    def __init__( self, publicID = None, **kwargs ):
        super( QPPublicObject, self ).__init__( **kwargs )
        QPPublicObject.publicObjectCtr += 1
        
        self.publicID = publicID


    def createPublicID( self, name=None, **kwargs ):
        """
        create publicID for a catalog object

        input:
                name: optional name string for publicID

        output:
                publicID

        kwargs:
            idstyle
                'short'   set ID without smi: format, either given string or ISO datetime
                'numeric' set publicObjectCtr as ID

        setting the ID style with a kwarg overrides QPPublicObject.publicIDStyle, but does not change the
        general setting
        """
        
        # current UTC timestamp from mxDateTime
        curr_time = utc()

        if (    ( 'idstyle' in kwargs.keys() and kwargs['idstyle'] == 'numeric' )
             or ( QPPublicObject.publicIDStyle == 'numeric' ) ):

            id = str( QPPublicObject.getPublicObjectCtr() )
            
        elif (    ( 'idstyle' in kwargs.keys() and kwargs['idstyle'] == 'short' )
               or ( QPPublicObject.publicIDStyle == 'short' ) ):
            if name is not None:
                id = str( name )
            else:
                id = mxDateTime2ISO( curr_time, 
                                     secondsdigits=self.secondsDigits, 
                                     timesepreplacechar='.' )

        else:
            if not name:
                name = 'generic'

            id = ''.join( ( 'smi:local/',
                             name,
                             '/',
                             mxDateTime2ISO( curr_time, 
                                             secondsdigits=self.secondsDigits,
                                             timesepreplacechar='.' ) ) )
            
        return id


    @classmethod
    def getPublicObjectCtr( cls ):
        return cls.publicObjectCtr


    @classmethod
    def setPublicIDStyle( cls, style ):

        # set publicIDStyle only if it is a valid entry
        if style in ( 'full', 'short', 'numeric' ):
            cls.publicIDStyle = style