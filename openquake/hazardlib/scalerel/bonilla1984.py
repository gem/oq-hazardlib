# The Hazard Library
# Copyright (C) 2012-2014, GEM Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Module :mod:`openquake.hazardlib.scalerel.bonilla1984` implements
    :class:`BonillaEtAl1984Aspect1p0`
    :class:`BonillaEtAl1984Aspect0p5`
    :class:`BonillaEtAl1984Aspect1p5`
    :class:`BonillaEtAl1984Aspect2p0`
    :class:`BonillaEtAl1984Aspect2p5`
    :class:`BonillaEtAl1984Aspect3p0`
    :class:`BonillaEtAl1984Aspect3p5`
    :class:`BonillaEtAl1984Aspect4p0`
    :class:`BonillaEtAl1984Aspect5p0`
    :class:`BonillaEtAl1984Aspect6p0`
    :class:`BonillaEtAl1984Aspect7p0`
    :class:`BonillaEtAl1984Aspect8p0`
"""
from math import log10, sqrt
from openquake.hazardlib.scalerel.base import BaseMSRSigma, BaseASRSigma

class BonillaEtAl1984Aspect1p0(BaseMSRSigma, BaseASRSigma):
    """
    T

    """
    ASPECT = 1.0
    
    def get_median_area(self, mag, rake):
        """
        Returns the median area using the "ordinary least squares" regression
        model between surface rupture length and magnitude. The "all.lm" option
        in Table 3.
        """
        return ((10.0 ** (-2.77 + 0.629 * mag)) ** 2.0) / self.ASPECT

    def get_std_dev_area(self, rake):
        """
        The original regression gives the standard deviation for length of
        0.286. In terms of area this becomes 0.572
        """
        return 0.572

    def get_median_mag(self, area, rake):
        """

        """
        length = sqrt(area * self.ASPECT)
        return 6.04 + 0.708 * log10(length)

    def get_std_dev_mag(self, rake):
        """

        """
        return 0.306


class BonillaEtAl1984Aspect0p5(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 0.5
    """
    ASPECT = 0.5


class BonillaEtAl1984Aspect1p5(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 1.5
    """
    ASPECT = 1.5


class BonillaEtAl1984Aspect2p0(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 2.0
    """
    ASPECT = 2.0


class BonillaEtAl1984Aspect2p5(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 2.5
    """
    ASPECT = 2.5


class BonillaEtAl1984Aspect3p0(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 3.0
    """
    ASPECT = 3.0


class BonillaEtAl1984Aspect3p5(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 3.5
    """
    ASPECT = 3.5


class BonillaEtAl1984Aspect4p0(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 4.0
    """
    ASPECT = 4.0


class BonillaEtAl1984Aspect5p0(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 5.0
    """
    ASPECT = 5.0


class BonillaEtAl1984Aspect6p0(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 6.0
    """
    ASPECT = 6.0


class BonillaEtAl1984Aspect7p0(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 7.0
    """
    ASPECT = 7.0


class BonillaEtAl1984Aspect8p0(BonillaEtAl1984Aspect1p0):
    """
    Bonilla et al. (1984) for the case when aspect ratio is 8.0
    """
    ASPECT = 8.0
