#!/usr/bin/env/python

'''Parser from nrml 0.4 to nhlib representation of objects'''

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
    def parse_source_model(self, filename, config):
        '''Read in the source model from nrml 0.4 format (using nrml tools)'''
        input_source = SourceModelParser(filename, config['schema_validation'])
        model1 = input_source.parse()
        output_model = []
        source_model = list(model1.sources)
        #self.number_sources = len(self.source_model)
        for source in source_model:
            output_model.append(self.nrml_to_nhlib(
                source,
                config['mesh_spacing'],
                config['bin_width'],
                config['area_discretisation']))
        return output_model
        
    def nrml_to_nhlib(self, src, mesh_spacing, bin_width, area_src_disc):
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
            return _area_to_nhlib(src, mesh_spacing, bin_width, area_src_disc)
        elif isinstance(src, nrml_models.PointSource):
            return _point_to_nhlib(src, mesh_spacing, bin_width)
        elif isinstance(src, nrml_models.ComplexFaultSource):
            return _complex_to_nhlib(src, mesh_spacing, bin_width)
        elif isinstance(src, nrml_models.SimpleFaultSource):
            return _simple_to_nhlib(src, mesh_spacing, bin_width)


def _point_to_nhlib(src)
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
    
def _area_to_nhlib


#~ def _point_to_nhlib(src, mesh_spacing, bin_width):
    #~ """Convert a NRML point source to the NHLib equivalent.

    #~ See :mod:`nrml.models` and :mod:`nhlib.source`.

    #~ :param src:
    #~ :class:`nrml.models.PointSource` instance.
    #~ :param float mesh_spacing: Rupture mesh spacing, in km.
    #~ :param float bin_width: Truncated Gutenberg-Richter MFD 
                            #~ (Magnitude Frequency Distribution) binwidth.
    #~ :returns: The NHLib representation of the input source.
    #~ """
    #~ shapely_pt = wkt.loads(src.geometry.wkt)

    #~ mf_dist = _mfd_to_nhlib(src.mfd, bin_width)

    #~ # nodal plane distribution:
    #~ npd = pmf.PMF(
        #~ [(x.probability,
          #~ geo.NodalPlane(strike=x.strike, dip=x.dip, rake=x.rake))
         #~ for x in src.nodal_plane_dist]
    #~ )

    #~ # hypocentral depth distribution:
    #~ hd = pmf.PMF([(x.probability, x.depth) for x in src.hypo_depth_dist])

    #~ point = mtkPointSource(
        #~ source_id=src.id,
        #~ name=src.name,
        #~ tectonic_region_type=src.trt,
        #~ mfd=mf_dist,
        #~ rupture_mesh_spacing=mesh_spacing,
        #~ magnitude_scaling_relationship=_SCALE_REL_MAP[src.mag_scale_rel](),
        #~ rupture_aspect_ratio=src.rupt_aspect_ratio,
        #~ upper_seismogenic_depth=src.geometry.upper_seismo_depth,
        #~ lower_seismogenic_depth=src.geometry.lower_seismo_depth,
        #~ location=geo.Point(shapely_pt.x, shapely_pt.y),
        #~ nodal_plane_distribution=npd,
        #~ hypocenter_distribution=hd
    #~ )

    #~ return point


def _area_to_nhlib(src, mesh_spacing, bin_width, area_src_disc):
    """Convert a NRML area source to the NHLib equivalent.

    See :mod:`nrml.models` and :mod:`nhlib.source`.
    :param src:
    :class:`nrml.models.PointSource` instance.
    :param float mesh_spacing:
    Rupture mesh spacing, in km.
    :param float bin_width:
    Truncated Gutenberg-Richter MFD (Magnitude Frequency Distribution) bin
    width.
    :param float area_src_disc:
    Area source discretization, in km. Applies only to area sources.
    :returns:
    The NHLib representation of the input source.
    """
    shapely_polygon = wkt.loads(src.geometry.wkt)
    nhlib_polygon = geo.Polygon(
        # We ignore the last coordinate in the sequence here, since it is a
        # duplicate of the first. nhlib will close the loop for us.
        [geo.Point(*x) for x in list(shapely_polygon.exterior.coords)[:-1]]
    )

    mf_dist = _mfd_to_nhlib(src.mfd, bin_width)

    # nodal plane distribution:
    npd = pmf.PMF(
        [(x.probability,
          geo.NodalPlane(strike=x.strike, dip=x.dip, rake=x.rake))
         for x in src.nodal_plane_dist]
    )

    # hypocentral depth distribution:
    hd = pmf.PMF([(x.probability, x.depth) for x in src.hypo_depth_dist])

    area = mtkAreaSource(
        source_id=src.id,
        name=src.name,
        tectonic_region_type=src.trt,
        mfd=mf_dist,
        rupture_mesh_spacing=mesh_spacing,
        magnitude_scaling_relationship=_SCALE_REL_MAP[src.mag_scale_rel](),
        rupture_aspect_ratio=src.rupt_aspect_ratio,
        upper_seismogenic_depth=src.geometry.upper_seismo_depth,
        lower_seismogenic_depth=src.geometry.lower_seismo_depth,
        nodal_plane_distribution=npd, hypocenter_distribution=hd,
        polygon=nhlib_polygon,
        area_discretization=area_src_disc
    )

    return area


def _simple_to_nhlib(src, mesh_spacing, bin_width):
    """Convert a NRML simple fault source to the NHLib equivalent.

    See :mod:`nrml.models` and :mod:`nhlib.source`.

    :param src:
    :class:`nrml.models.PointSource` instance.
    :param float mesh_spacing:
    Rupture mesh spacing, in km.
    :param float bin_width:
    Truncated Gutenberg-Richter MFD (Magnitude Frequency Distribution) bin
    width.
    :returns:
    The NHLib representation of the input source.
    """
    shapely_line = wkt.loads(src.geometry.wkt)
    fault_trace = geo.Line([geo.Point(*x) for x in shapely_line.coords])

    mf_dist = _mfd_to_nhlib(src.mfd, bin_width)

    simple = mtkSimpleFaultSource(
        source_id=src.id,
        name=src.name,
        tectonic_region_type=src.trt,
        mfd=mf_dist,
        rupture_mesh_spacing=mesh_spacing,
        magnitude_scaling_relationship=_SCALE_REL_MAP[src.mag_scale_rel](),
        rupture_aspect_ratio=src.rupt_aspect_ratio,
        upper_seismogenic_depth=src.geometry.upper_seismo_depth,
        lower_seismogenic_depth=src.geometry.lower_seismo_depth,
        fault_trace=fault_trace,
        dip=src.geometry.dip,
        rake=src.rake
    )

    return simple


def _complex_to_nhlib(src, mesh_spacing, bin_width):
    """Convert a NRML complex fault source to the NHLib equivalent.

    See :mod:`nrml.models` and :mod:`nhlib.source`.

    :param src:
    :class:`nrml.models.PointSource` instance.
    :param float mesh_spacing:
    Rupture mesh spacing, in km.
    :param float bin_width:
    Truncated Gutenberg-Richter MFD (Magnitude Frequency Distribution) bin
    width.
    :returns:
    The NHLib representation of the input source.
    """
    edges_wkt = []
    edges_wkt.append(src.geometry.top_edge_wkt)
    edges_wkt.extend(src.geometry.int_edges)
    edges_wkt.append(src.geometry.bottom_edge_wkt)

    edges = []

    for edge in edges_wkt:
        shapely_line = wkt.loads(edge)
        line = geo.Line([geo.Point(*x) for x in shapely_line.coords])
        edges.append(line)

    mf_dist = _mfd_to_nhlib(src.mfd, bin_width)

    cmplx = mtkComplexFaultSource(
        source_id=src.id,
        name=src.name,
        tectonic_region_type=src.trt,
        mfd=mf_dist,
        rupture_mesh_spacing=mesh_spacing,
        magnitude_scaling_relationship=_SCALE_REL_MAP[src.mag_scale_rel](),
        rupture_aspect_ratio=src.rupt_aspect_ratio,
        edges=edges,
        rake=src.rake
    )

    return cmplx


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
