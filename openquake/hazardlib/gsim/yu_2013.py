# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2014-2016 GEM Foundation, Chung-Han Chan
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.

"""
Module exports :class:`yu2013`.
"""
from __future__ import division
from __future__ import print_function

import numpy as np
from scipy.constants import g

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA


class YuEtAl2013(GMPE):
    """
    """

    #: Supported tectonic region type is subduction interface along the
    #: Sumatra subduction zone.
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types are spectral acceleration,
    #: peak ground velocity and peak ground acceleration, see table IV
    #: pag. 837
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
        PGV
    ])

    #: Supported intensity measure component is geometric mean
    #: of two horizontal components,
    #####: PLEASE CONFIRM!!!!! 140709
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    #: Supported standard deviation types is total, see equation IV page 837.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: Required site parameter is only Vs30 (used to distinguish rock
    REQUIRES_SITES_PARAMETERS = set(())

    #: Required rupture parameters are magnitude, and focal depth
    REQUIRES_RUPTURE_PARAMETERS = set(('mag',))

    #: Required distance measure is hypocentral distance,
    REQUIRES_DISTANCES = set(('rhypo', 'azimuth',))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # Check that the requested standard deviation type is available
        assert all(stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
                   for stddev_type in stddev_types)
        # Get the coefficients
        C = self.COEFFS[imt]
        if rup.mag > 6.5:
            ca = C['ua']
            cb = C['ub']
            cc = C['uc']
            cd = C['ud']
            ce = C['ue']
        else:
            ca = C['a']
            cb = C['b']
            cc = C['c']
            cd = C['d']
            ce = C['e']
        # Compute the mean value (i.e. natural logarithm of ground motion)
        mag = rup.mag
        mean = (ca +
                cb * mag + dists.azimuth +
                cc * np.log(dists.rhypo + cd * np.exp(ce * mag)))

        # Convert acceleration from gals into fraction of g
        if isinstance(imt, (PGA, SA)):
            mean = np.log(np.exp(mean)/g)
        # Get the standard deviation
        stddevs = self._compute_std(C, stddev_types, len(dists.rhypo))
        # Return results
        return mean, stddevs

    def _compute_std(self, C, stddev_types, num_sites):
        """
        """
        return C['sigma']

    #: Coefficient table for rock sites, see table 3 page 227.
    COEFFS = CoeffsTable(sa_damping=5, table="""\
        IMT           a        b          c          d           e      ua       ub         uc         ud          ue    sigma
        PGA       2.369   2.0852   -0.23564   -0.87906   -0.001363   2.369   2.0852   -0.23564   -0.87906   -0.001363   0.3478
        PGV       2.369   2.0852   -0.23564   -0.87906   -0.001363   2.369   2.0852   -0.23564   -0.87906   -0.001363   0.3478
    """)


class YuEtAl2013StableWenchuan(YuEtAl2013):

    #: Supported tectonic region type is subduction interface along the
    #: Sumatra subduction zone.
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL

    #: Coefficient table for rock sites, see table 3 page 227.
    COEFFS = CoeffsTable(sa_damping=5, table="""\
        IMT           a        b          c          d           e      ua       ub         uc         ud          ue    sigma
        PGA       2.369   2.0852   -0.23564   -0.87906   -0.001363   2.369   2.0852   -0.23564   -0.87906   -0.001363   0.3478
        PGV       2.369   2.0852   -0.23564   -0.87906   -0.001363   2.369   2.0852   -0.23564   -0.87906   -0.001363   0.3478
    """)

