# coding: utf-8
# The Hazard Library
# Copyright (C) 2012 GEM Foundation
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
Module exports :class:`BindiEtAl2011`.
"""
from __future__ import division

import numpy as np

from scipy.constants import g

from openquake.hazardlib.gsim.base import CoeffsTable, GMPE
from openquake.hazardlib import const
from openquake.hazardlib.imt import MMI


class BindiEtAl2011(IPE):
    """
    Implements IPE developed by Dino Bindi et al. 2011 and published
    as "Intensity prediction equations for Central Asia" (Geo-physical journal international, 2011, 187,327-337).
for a fixed depth of 15 km and epicentral distance (equation 5 in the paper) 

    """

    #: Supported tectonic region type is active shallow crust,
    #: see end of 'Introduction', page 454.
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types are MSK-64

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        MMI
    ])


    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.HORIZONTAL


    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: Required rupture parameters are magnitude 
    REQUIRES_RUPTURE_PARAMETERS = set(('mag'))

    #: Required distance measure is joyner-boore, 
    REQUIRES_DISTANCES = set(('rjb', ))

    def get_mean_and_stddevs(self, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # extract dictionaries of coefficients specific to required
        # intensity measure type
        C = self.COEFFS[imt]

        mean = self._compute_mean(C, rup.mag, dists.rjb, imt)

        stddevs = self._get_stddevs(C, stddev_types)

        return mean, stddevs  

    def _compute_mean(self, C, mag, rjb):
        """
        Compute mean value for MSK-64.
        """
        mean = C['a1']*mag +C['a2']-C['a3'] * np.log10 (np.sqrt((rjb**2 + 15**2) / 15**2))-C['a4'] * 
               (np.sqrt(rjb**2 + 15**2)-15)

        return mean

    def _get_stddevs(self, C, stddev_types, num_sites):
        """
        Return total standard deviation.
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            stddevs.append( C['sigma'] + np.zeros(num_sites))

        return stddevs

    #: Coefficient table constructed from the electronic suplements of the
    #: original paper.
    COEFFS = CoeffsTable(table="""\
    IMT       a1     a2      a3      a4      sigma
    MMI      1.049  0.686  2.706   0.0001811   0.689
        """)
