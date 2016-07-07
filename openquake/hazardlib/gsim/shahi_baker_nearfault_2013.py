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
Module exports :class:`ShahiBakerNearFault2013`
"""
from __future__ import division

import numpy as np
from math import log, exp
from openquake.hazardlib.gsim.shahi_baker_2013 import ShahiBaker2013
from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import SA


class ShahiBakerNearFault2013(ShahiBaker2013):
    """
    Implements GMPE developed by Shrey Kumar Shahi published as "A
    probabilistic framework to include the effects of near-fault
    directivity in seismic hazard assessment" (2013, Doctoral dissertation,
    Stanford University, 102-140). The GMPE is developed based on Kenneth W.
    Campbell and Yousef Bozorgnia (2008) GMPE but explicitly includes the
    near-fault directivity effects. This class implements the model for
    RotD50 component of the elastic spectra. This implement is the model
    with unknown pulse period and unknown pulse indicator.
    """
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
                self._compute_basin_response_term(C, sites.z2pt5) +
                self._compute_pulse_response_term(C, rup.mag, imt, rup, dists)
                )

        stddevs = self._get_stddevs(C)
        return mean, stddevs

    def _get_stddevs(self, C):
        """
        Returns the standard deviations. Table 5.5, pages 139
        """

        stddevs = []

        stddevs.append(np.array(C['std_unknown']))

        return stddevs

    def _compute_pulse_response_term(self, C, mag, imt, rup, dists):
        """
        Returns the pulse response term (equation 5.20, page 119)
        """
        mean_lntp = -6.207 + 1.075 * mag
        sigma_lntp = 0.61
        d = 2 * sigma_lntp ** 2 * C['b1']
        c = np.log(imt.period) - C['b2']
        alpha_pulse = ((d * c ** 2 - mean_lntp ** 2) / (d - 1)) - (
            (d * c - mean_lntp) / (d - 1)) ** 2
        exp_log = -1 * (1 - d) * alpha_pulse / (2. * sigma_lntp ** 2)
        pulse_amp = C['b0'] * np.exp(exp_log) / (1. - d) ** 0.5

        prob_pulse = np.zeros_like(rup.rake, dtype=float)
        frv, fnm = self._get_fault_type_dummy_variables(rup.rake)

        if frv == fnm == 0.:
            prob_pulse = 1. / (
                1 + np.exp(0.79 + 0.138 * dists.rjb - 0.353 * dists.rs ** 0.5
                           + 0.02 * dists.rtheta))
        else:
            prob_pulse = 1. / (
                1 + np.exp(1.483 + 0.124 * dists.rjb - 0.688 * dists.rd ** 0.5
                           + 0.022 * dists.rphi))
        return pulse_amp * prob_pulse

    COEFFS = CoeffsTable(sa_damping=5, table="""\
imt     c0      c1      c2      c3      c4      c5      c6      c7      c8      c9      c10     c11     c12     k1      k2      k3      b0      b1      b2      std_unknown n     c
0.01    -1.729  0.51    -0.44   -0.396  -2.278  0.185   5.779   0.227   -0.493  0.675   1.017   -0.03   0.292   865     -1.186  1.839   0.72    -1.1    -0.19   0.541       1.18  1.88
0.02    -1.693  0.519   -0.455  -0.39   -2.299  0.186   5.758   0.214   -0.51   0.697   1.062   -0.027  0.282   865     -1.219  1.84    0.72    -1.1    -0.19   0.546       1.18  1.88
0.03    -1.569  0.509   -0.418  -0.419  -2.321  0.187   5.766   0.199   -0.516  0.731   1.125   -0.036  0.31    908     -1.273  1.841   0.72    -1.1    -0.19   0.558       1.18  1.88
0.05    -1.227  0.523   -0.389  -0.491  -2.467  0.201   5.89    0.143   -0.584  0.803   1.248   -0.032  0.342   1054    -1.346  1.843   0.72    -1.1    -0.19   0.603       1.18  1.88
0.075   -0.67   0.536   -0.397  -0.571  -2.673  0.22    7.116   0.092   -0.707  0.819   1.442   -0.077  0.35    1086    -1.471  1.845   0.72    -1.1    -0.19   0.655       1.18  1.88
0.1     -0.327  0.526   -0.383  -0.611  -2.705  0.225   8.033   0.082   -0.657  0.813   1.546   -0.054  0.318   1032    -1.624  1.847   0.72    -1.1    -0.19   0.668       1.18  1.88
0.15    -0.145  0.511   -0.428  -0.531  -2.606  0.206   8.8     0.116   -0.5    0.765   1.789   -0.067  0.424   878     -1.931  1.852   0.72    -1.1    -0.19   0.62        1.18  1.88
0.2     -0.496  0.495   -0.393  -0.473  -2.352  0.182   7.704   0.162   -0.297  0.706   2       -0.008  0.445   748     -2.188  1.856   0.72    -1.1    -0.19   0.592       1.18  1.88
0.25    -0.894  0.508   -0.363  -0.474  -2.246  0.172   6.695   0.226   -0.23   0.642   2.149   0.001   0.345   654     -2.381  1.861   0.72    -1.1    -0.19   0.574       1.18  1.88
0.3     -1.158  0.525   -0.425  -0.331  -2.177  0.172   6.183   0.233   -0.259  0.676   2.273   -0.034  0.292   587     -2.518  1.865   0.72    -1.1    -0.19   0.59        1.18  1.88
0.4     -1.443  0.524   -0.34   -0.307  -2.135  0.173   5.493   0.298   -0.2    0.619   2.408   0.03    0.237   503     -2.657  1.874   0.72    -1.1    -0.19   0.61        1.18  1.88 
0.5     -2.534  0.679   -0.471  -0.207  -2.032  0.16    5.024   0.316   -0.141  0.59    2.369   0.006   0.183   457     -2.669  1.883   0.72    -1.1    -0.19   0.634       1.18  1.88
0.75    -4.8    1       -0.778  -0.13   -1.986  0.159   4.294   0.354   -0.076  0.647   1.978   0.002   0.101   410     -2.401  1.906   0.72    -1.1    -0.19   0.675       1.18  1.88
1       -6.388  1.207   -0.835  -0.126  -1.992  0.156   4.207   0.385   0.097   0.545   1.467   0.094   0.14    400     -1.955  1.929   0.72    -1.1    -0.19   0.685       1.18  1.88
1.5     -8.705  1.503   -0.956  -0.204  -1.997  0.157   4.219   0.358   0.008   0.527   0.373   0.215   0.085   400     -1.025  1.974   0.72    -1.1    -0.19   0.665       1.18  1.88
2       -9.742  1.581   -0.743  -0.444  -2.014  0.157   4.409   0.343   0.024   0.402   -0.471  0.212   0.034   400     -0.299  2.019   0.72    -1.1    -0.19   0.646       1.18  1.88
3       -10.652 1.564   -0.114  -0.966  -1.841  0.133   4.633   0.26    -0.054  0.129   -0.795  0.148   0.01    400     0.      2.11    0.72    -1.1    -0.19   0.643       1.18  1.88
4       -11.323 1.543   0.264   -1.223  -1.794  0.124   4.647   0.152   -0.032  -0.112  -0.783  0.087   0.053   400     0.      2.2     0.72    -1.1    -0.19   0.637       1.18  1.88
5       -11.797 1.535   0.547   -1.454  -1.677  0.116   4.723   -0.021  -0.122  -0.23   -0.756  0.01    0.038   400     0.      2.291   0.72    -1.1    -0.19   0.647       1.18  1.88
7.5     -12.623 1.487   0.662   -1.312  -1.632  0.126   4.52    -0.096  -0.494  -0.354  -0.669  0.028   0.068   400     0.      2.517   0.72    -1.1    -0.19   0.669       1.18  1.88
10      -13.128 1.552   0.298   -1.006  -1.733  0.128   4.165   -0.294  -0.358  -0.383  -0.476  0.195   0.08    400     0.      2.744   0.72    -1.1    -0.19   0.635       1.18  1.88
    """)
