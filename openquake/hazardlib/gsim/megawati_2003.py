# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2014-2016 GEM Foundation
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
Module exports :class:`MegawatiEtAl2003`.
"""
from __future__ import division
from __future__ import print_function

import numpy as np
from scipy.constants import g

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA


class MegawatiEtAl2003(GMPE):
    """
    Implements GMPE developed by Megawati et al. published as "Response
    spectral attenuation relationships for Singapore and the Malay Peninsula
    due to distant Sumatran-fault earthquakes" (2002, Earthquake Engineering &
    Structural Dynamics Volume 32, pages 2241–2265).
    """

    #: Supported tectonic region type is subduction interface along the
    #: Sumatra subduction zone.
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types are spectral acceleration,
    #: peak ground velocity and peak ground acceleration, see table IV
    #: pag. 837
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
        PGV,
        SA
    ])

    #: Supported intensity measure component is geometric mean
    #: of two horizontal components,
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    #: Supported standard deviation types is total, see equation IV page 837.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: Required site parameter is only Vs30 (used to distinguish rock
    #: and deep soil).
    #: This GMPE is for very hard rock site condition,
    #: see the abstract page 827.
    REQUIRES_SITES_PARAMETERS = set(())

    #: Required rupture parameters are magnitude, and focal depth, see
    #: equation 10 page 226.
    REQUIRES_RUPTURE_PARAMETERS = set(('mag',))

    #: Required distance measure is hypocentral distance,
    #: see equation 1 page 834.
    REQUIRES_DISTANCES = set(('rrup', 'azimuth'))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # Check that the definitions of standard deviations is in agreement
        # with the standards used in the OQ hazard library
        assert all(stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
                   for stddev_type in stddev_types)
        # Get coefficients
        C = self.COEFFS[imt]
        # Compute the mean
        mean = (self._get_magnitude_scaling(C, rup.mag) +
                self._get_distance_scaling(C, rup.mag, dists.rhypo))
        # Convert acceleration from cm/s^2 into fractions of g
        if isinstance(imt, (PGA, SA)):
            mean = np.log(np.exp(mean) / (100.0 * g))
        # Compute the standard deviation
        stddevs = self._compute_std(C, stddev_types, len(dists.rhypo))
        # Return results
        return mean, stddevs

    def _get_magnitude_scaling(self, C, mag):
        """
        Returns the magnitude scaling term
        """
        return (C["a0"] +
                C["a1"] * mag +
                C["a2"] * mag** 2.)

    def _get_distance_scaling(self, C, mag, r):
        """
        Returns the distance scaling term
        """
        return (C["a3"] * np.log(r) +
                C["a4"] * r)

    def _get_azimuthal_correction(self, azimuth):
        """
        Returns the azimuthal correction
        """
        return np.log(max(np.abs(cos(2*azimuth)),
                          C["a5"]*2*np.abs(sin(2*azimuth)) ))

    def _get_std(self, C, stddev_types, num_sites):
        """
        Get standard deviation. This GMPE supports only total standard
        deviation.
        """
        return [np.ones(num_sites) * C['sigma']]

    #: Coefficient table for rock sites, see table 3 page 227.
    COEFFS = CoeffsTable(sa_damping=5, table="""\
        IMT       a0       a1         a2         a3          a4        a5    sigma
        PGA   -8.167   2.7779  -0.045945   -1.00000   -0.001906    0.1356   0.3511
    """)
