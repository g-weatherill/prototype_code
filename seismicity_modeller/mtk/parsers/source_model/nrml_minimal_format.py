#!/usr/bin/env/python

'''Parser from nrml 0.4 to nhlib representation of objects'''
import numpy as np
from nhlib import geo
from nhlib import mfd
from nhlib import pmf
from nhlib import scalerel
from sources.mtk_point import mtkPointSource
from sources.mtk_area import mtkAreaSource
from sources.mtk_simple_fault import mtkSimpleFaultSource
from sources.mtk_complex_fault import mtkComplexFaultSource
#from nhlib import source
from nrml import models as nrml_models
from base import BaseSourceParser
from nrml.parsers import SourceModelParser
from shapely import wkt

_SCALE_REL_MAP = {
    'PeerMSR': scalerel.PeerMSR,
    'WC1994': scalerel.WC1994,
}


class nrmlSourceModelParser(BaseSourceParser):
    '''Parser for the nrml 0.4 source model'''
    def read_source_model(self, filename, config):
        '''Read in the source model from nrml 0.4 format (using nrml tools)'''
        input_source = SourceModelParser(filename, config['schema_validation'])
        model1 = input_source.parse()
        output_model = []
        source_model = list(model1.sources)
        #self.number_sources = len(self.source_model)
        for source in source_model:
            output_model.append(self.nrml_to_nhlib(
                source,
                config['mesh_spacing']))
                #config['bin_width'],
                #config['area_discretisation']))
        return output_model
        
    def nrml_to_nhlib(self, src, mesh_spacing=1.0): #bin_width, area_src_disc):
        """Convert a seismic source object from the NRML representation to the
        NHLib representation. Inputs can be point, area, simple fault, or 
        complex fault sources.

        See :mod:`nrml.models` and :mod:`nhlib.source`.

        :param src:
        :mod:`nrml.models` seismic source instance.
        :param float mesh_spacing:
        Rupture mesh spacing, in km.
        :param float bin_width:
        Truncated Gutenberg-Richter MFD (Magnitude Frequency Distribution) bin
        width.
        :param float area_src_disc:
        Area source discretization, in km. Applies only to area sources.
        If the input source is known to be a type other than an area source,
        you can specify `area_src_disc=None`.
        :returns:
        The NHLib representation of the input source.
        """
        # The ordering of the switch here matters because:
        # - AreaSource inherits from PointSource
        # - ComplexFaultSource inherits from SimpleFaultSource
        if isinstance(src, nrml_models.AreaSource):
            return _area_to_mtk(src)
        elif isinstance(src, nrml_models.PointSource):
            return _point_to_mtk(src)
        elif isinstance(src, nrml_models.ComplexFaultSource):
            return _complex_to_mtk(src, mesh_spacing)
        elif isinstance(src, nrml_models.SimpleFaultSource):
            return _simple_to_mtk(src, mesh_spacing)


def _point_to_mtk(src):
    """Convert a NRML point source to the MTK Source Equivalent.

    See :mod:`nrml.models` and :mod:`nhlib.source`.

    :param src:
    :class:`nrml.models.PointSource` instance.
    """
    # Loads point geometry
    shapely_pt = wkt.loads(src.geometry.wkt)
    geometry_input = np.array([shapely_pt.x, shapely_pt.y])
    point = mtkPointSource(
        src.id, 
        src.name, 
        src.trt, 
        src.rupt_aspect_ratio,
        geometry_input,
        upper_depth = src.geometry.upper_seismo_depth,
        lower_depth=src.geometry.lower_seismo_depth)
    return point
    
def _area_to_mtk(src):
    """Convert a NRML area source to the NHLib equivalent.

    See :mod:`nrml.models` and :mod:`nhlib.source`.
    :param src:
    :class:`nrml.models.AreaSource` instance.
    :returns:
    The NHLib representation of the input source.
    """
    shapely_poly = wkt.loads(src.geometry.wkt)
    geometry_input = np.array(list(shapely_poly.exterior.coords))
    area = mtkAreaSource(
        src.id, 
        src.name, 
        src.trt, 
        src.rupt_aspect_ratio,
        geometry_input,
        upper_depth=src.geometry.upper_seismo_depth,
        lower_depth=src.geometry.lower_seismo_depth)   
    return area



def _simple_to_mtk(src, fault_mesh_spacing=1.0):
    """Convert a NRML simple fault source to the MTK equivalent.

    See :mod:`nrml.models` and :mod:`nhlib.source`.

    :param src:
    :class:`nrml.models.SimpleFaultSource` instance.
    """
    shapely_line = wkt.loads(src.geometry.wkt)
    input_geometry = np.array(shapely_line.coords)
    simple_fault = mtkSimpleFaultSource(
        src.id,
        src.name,
        src.trt,
        src.rupt_aspect_ratio,
        input_geometry,
        src.geometry.dip,
        src.geometry.upper_seismo_depth,
        src.geometry.lower_seismo_depth,
        mesh_spacing=fault_mesh_spacing)
    return simple_fault
    
def _complex_to_mtk(src, fault_mesh_spacing=1.0):
    """Convert a NRML complex fault source to the MTK equivalent.

    See :mod:`nrml.models` and :mod:`nhlib.source`.

    :param src:
    :class:`nrml.models.ComplexFaultSource` instance.
    :param float mesh_spacing:
    Rupture mesh spacing, in km.
    :returns:
    The NHLib representation of the input source.
    """   
    # Parse geometry - top edge
    input_geometry = [_line_wkt_to_array(src.geometry.top_edge_wkt)]
    
    if src.geometry.int_edges:
        # Source has intermediate edges
        for edge in src.geometry.int_edges:
            #int_edge = wkt.loads(edge)
            input_geometry.append(
                _line_wkt_to_array(edge))
    
    # Parse bottom edge
    input_geometry.append(
        _line_wkt_to_array(src.geometry.bottom_edge_wkt))
    #print input_geometry    
    complex_fault = mtkComplexFaultSource(
        src.id,
        src.name,
        src.trt,
        src.rupt_aspect_ratio,
        input_geometry,
        mesh_spacing = fault_mesh_spacing)
    return complex_fault

def _line_wkt_to_array(linestring):
    '''Converts a wkt string to a numpy array using regular expressions'''
    substring = linestring[
        linestring.find("(")+1:linestring.find(")")]
    data = substring.split(",")
    array_out = np.empty([len(data), 3], dtype=float)
    for iloc, datum in enumerate(data):
        values = np.array(datum.split())
        array_out[iloc, :] = values.astype(float)
    return array_out
        

def _mfd_to_nhlib(src_mfd, bin_width):
    """Convert a NRML MFD to an NHLib MFD.

    :param src_mfd:
    :class:`nrml.models.IncrementalMFD` or :class:`nrml.models.TGRMFD`
    instance.
    :param float bin_width:
    Optional. Required only for Truncated Gutenberg-Richter MFDs.
    :returns:
    The NHLib representation of the MFD. See :mod:`nhlib.mfd`.
    """
    if isinstance(src_mfd, nrml_models.TGRMFD):
        assert bin_width is not None
        return mfd.TruncatedGRMFD(
            a_val=src_mfd.a_val, b_val=src_mfd.b_val, 
            min_mag=src_mfd.min_mag, max_mag=src_mfd.max_mag, 
            bin_width=bin_width
        )
    elif isinstance(src_mfd, nrml_models.IncrementalMFD):
        return mfd.EvenlyDiscretizedMFD(
            min_mag=src_mfd.min_mag, bin_width=src_mfd.bin_width,
            occurrence_rates=src_mfd.occur_rates
        )


def _source_type(src_model):
    """Given of the source types defined in :mod:`nrml.models`, get the
    `source_type` for a :class:`~openquake.db.models.ParsedSource`.
    """
    if isinstance(src_model, nrml_models.AreaSource):
        return 'area'
    elif isinstance(src_model, nrml_models.PointSource):
        return 'point'
    elif isinstance(src_model, nrml_models.ComplexFaultSource):
        return 'complex'
    elif isinstance(src_model, nrml_models.SimpleFaultSource):
        return 'simple'
