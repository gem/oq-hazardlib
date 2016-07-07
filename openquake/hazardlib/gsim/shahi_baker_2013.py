# The Hazard Library
# Copyright (C) 2013-2014, GEM Foundation
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
Module exports :class:`ShahiBaker2013`
"""
from __future__ import division

from openquake.hazardlib.gsim.campbell_bozorgnia_2008 import CampbellBozorgnia2008
import numpy as np
from math import log, exp
from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import SA


class ShahiBaker2013(CampbellBozorgnia2008):
    """
    Implements GMPE developed by Shrey Kumar Shahi published as "A
    probabilistic framework to include the effects of near-fault
    directivity in seismic hazard assessment" (2013, Doctoral dissertation,
    Stanford University, 102-140). The GMPE is developed based on Kenneth W.
    Campbell and Yousef Bozorgnia (2008) GMPE but refit with NGA-West2 database.
    """

    #: Supported tectonic region type is active shallow crust
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types is spectral acceleration
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        SA
    ])

    #: Supported intensity measure component is orientation-independent
    #: average horizontal :attr:`~openquake.hazardlib.const.IMC.RotD50`
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.RotD50

    #: Supported standard deviation type is total
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: Required distance measures are Rrup, Rjb, Rs, Rtheta, Rd, and Rphi
    REQUIRES_DISTANCES = set(('rrup', 'rjb', 'rs', 'rtheta', 'rd', 'rphi'))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # extract dictionaries of coefficients specific to required
        # intensity measure type
        C = self.COEFFS[imt]

        # compute median pga on rock (vs30=1100), needed for site response
        # term calculation
        # For spectral accelerations at periods between 0.0 and 0.25 s, Sa (T)
        # cannot be less than PGA on soil, therefore if the IMT is in this
        # period range it is necessary to calculate PGA on soil
        if isinstance(imt, SA) and (imt.period > 0.0) and (imt.period < 0.25):
            get_pga_site = True
        else:
            get_pga_site = False
        pga1100, pga_site = self._compute_imt1100(C,
                                                  sites,
                                                  rup,
                                                  dists,
                                                  get_pga_site)

        # Get the median ground motion
        mean = (self._compute_magnitude_term(C, rup.mag) +
                self._compute_distance_term(C, rup, dists) +
                self._compute_style_of_faulting_term(C, rup) +
                self._compute_hanging_wall_term(C, rup, dists) +
                self._compute_shallow_site_response(C, sites, pga1100) +
                self._compute_basin_response_term(C, sites.z2pt5)
                )

        stddevs = self._get_stddevs(C)

        return mean, stddevs

    def _compute_imt1100(self, C, sites, rup, dists, get_pga_site=False):
        """
        Computes the PGA on reference (Vs30 = 1100 m/s) rock.
        """
        # Calculates simple site response term assuming all sites 1100 m/s
        fsite = (C['c10'] + (C['k2'] * C['n'])) * log(1100. / C['k1'])
        # Calculates the PGA on rock
        pga1100 = np.exp(self._compute_magnitude_term(C, rup.mag) +
                         self._compute_distance_term(C, rup, dists) +
                         self._compute_style_of_faulting_term(C, rup) +
                         self._compute_hanging_wall_term(C, rup, dists) +
                         self._compute_basin_response_term(C, sites.z2pt5) +
                         fsite)
        # If PGA at the site is needed then remove factor for rock and
        # re-calculate on correct site condition
        if get_pga_site:
            pga_site = np.exp(np.log(pga1100) - fsite)
            fsite = self._compute_shallow_site_response(C, sites, pga1100)
            pga_site = np.exp(np.log(pga_site) + fsite)
        else:
            pga_site = None
        return pga1100, pga_site

    def _get_stddevs(self, C):
        """
        Returns the standard deviations. Table 5.5, pages 139
        """
        stddevs = []

        stddevs.append(np.array(C['std']))
        return stddevs

    COEFFS = CoeffsTable(sa_damping=5, table="""\
imt     c0      c1      c2      c3      c4      c5      c6      c7      c8      c9      c10     c11     c12     std     k1      k2      k3      n       c
0.01    -1.729  0.51    -0.444  -0.378  -2.208  0.178   5.779   0.239   -0.426  0.675   1.017   -0.03   0.292   0.54    865     -1.186  1.839   1.18    1.88
0.02    -1.693  0.511   -0.45   -0.377  -2.222  0.18    5.758   0.23    -0.446  0.697   1.062   -0.027  0.282   0.544   865     -1.219  1.84    1.18    1.88
0.03    -1.569  0.509   -0.426  -0.396  -2.258  0.179   5.766   0.21    -0.442  0.731   1.125   -0.036  0.31    0.557   908     -1.273  1.841   1.18    1.88
0.05    -1.227  0.519   -0.412  -0.441  -2.363  0.188   5.89    0.164   -0.507  0.803   1.248   -0.032  0.342   0.6     1054    -1.346  1.843   1.18    1.88
0.075   -0.67   0.527   -0.433  -0.493  -2.546  0.203   7.116   0.129   -0.585  0.819   1.442   -0.077  0.35    0.649   1086    -1.471  1.845   1.18    1.88
0.1     -0.327  0.519   -0.439  -0.496  -2.552  0.2     8.033   0.12    -0.51   0.813   1.546   -0.054  0.318   0.66    1032    -1.624  1.847   1.18    1.88
0.15    -0.145  0.504   -0.46   -0.453  -2.487  0.188   8.8     0.146   -0.404  0.765   1.789   -0.067  0.424   0.617   878     -1.931  1.852   1.18    1.88
0.2     -0.496  0.495   -0.394  -0.454  -2.302  0.173   7.704   0.167   -0.247  0.706   2       -0.008  0.445   0.592   748     -2.188  1.856   1.18    1.88
0.25    -0.894  0.511   -0.332  -0.456  -2.219  0.167   6.695   0.23    -0.167  0.642   2.149   0.001   0.345   0.575   654     -2.381  1.861   1.18    1.88
0.3     -1.158  0.519   -0.347  -0.377  -2.117  0.162   6.183   0.22    -0.207  0.676   2.273   -0.034  0.292   0.591   587     -2.518  1.865   1.18    1.88
0.4     -1.443  0.525   -0.287  -0.361  -2.101  0.165   5.493   0.285   -0.158  0.619   2.408   0.03    0.237   0.609   503     -2.657  1.874   1.18    1.88
0.5     -2.534  0.689   -0.461  -0.236  -2.017  0.151   5.024   0.316   -0.126  0.59    2.369   0.006   0.183   0.633   457     -2.669  1.883   1.18    1.88
0.75    -4.8    1.025   -0.767  -0.136  -1.985  0.155   4.294   0.354   -0.068  0.647   1.978   0.002   0.101   0.674   410     -2.401  1.906   1.18    1.88
1       -6.388  1.237   -0.83   -0.14   -1.99   0.144   4.207   0.395   0.105   0.545   1.467   0.094   0.14    0.684   400     -1.955  1.929   1.18    1.88
1.5     -8.705  1.536   -0.953  -0.197  -1.999  0.148   4.219   0.366   0.021   0.527   0.373   0.215   0.085   0.664   400     -1.025  1.974   1.18    1.88
2       -9.742  1.617   -0.735  -0.427  -2.015  0.148   4.409   0.358   0.032   0.402   -0.471  0.212   0.034   0.645   400     -0.299  2.019   1.18    1.88
3       -10.652 1.606   -0.106  -0.937  -1.857  0.124   4.633   0.264   -0.043  0.129   -0.795  0.148   0.01    0.642   400     0       2.11    1.18    1.88
4       -11.323 1.581   0.276   -1.175  -1.821  0.119   4.647   0.152   -0.02   -0.112  -0.783  0.087   0.053   0.636   400     0       2.2     1.18    1.88
5       -11.797 1.583   0.532   -1.386  -1.719  0.11    4.723   -0.024  -0.091  -0.23   -0.756  0.01    0.038   0.644   400     0       2.291   1.18    1.88
7.5     -12.623 1.548   0.525   -1.123  -1.73   0.126   4.52    -0.094  -0.404  -0.354  -0.669  0.028   0.068   0.666   400     0       2.517   1.18    1.88
10      -13.128 1.608   0.166   -0.811  -1.752  0.125   4.165   -0.301  -0.224  -0.383  -0.476  0.195   0.08    0.635   400     0       2.744   1.18    1.88
    """)
