"""Module  :mod: scalerel.wc1994 implements class WC1994"""
from math import log10
from scalerel.base import BaseASR

class WC1994(BaseASR):
    """
    Wells and Coppersmith magnitude -- rupture area relationships,
    see 1994, Bull. Seism. Soc. Am., pages 974-2002.

    Implements both magnitude-area and area-magnitude scaling relationships.
    """
    def get_std_dev_mag(self, rake):
        """
        Standard deviation on the magnitude for the WC1994 area relation.
        """
        assert rake is None or -180 <= rake <= 180
        if rake is None:
            # their "All" case
            return 0.24
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            return 0.23
        elif rake > 0:
            # thrust/reverse
            return 0.25
        else:
            # normal
            return 0.25

    def get_median_mag(self, area, rake):
        """
        Return magnitude (Mw) given the area and rake.
        Setting the rake to ``None`` causes their "All" rupture-types
        to be applied.

        :param area:
        Area in square km.
        :param rake:
        Rake angle (the rupture propagation direction) in degrees,
        from -180 to 180.
        """
        assert rake is None or -180 <= rake <= 180
        if rake is None:
            # their "All" case
            return 4.07 + 0.98 * log10(area)
        elif (-45 <= rake <= 45) or (rake > 135) or (rake < -135):
            # strike slip
            return 3.98 + 1.02 * log10(area)
        elif rake > 0:
            # thrust/reverse
            return 4.33 + 0.90 * log10(area)
        else:
            # normal
            return 3.93 + 1.02 * log10(area)