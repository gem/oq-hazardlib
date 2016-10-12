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
from openquake.hazardlib.geo import geodetic

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
     #: and deep soil).
     #: This GMPE is for very hard rock site condition,
     #: see the abstract page 827.
     REQUIRES_SITES_PARAMETERS = set(())

     #: Required rupture parameters are magnitude, and focal depth, see
     #: equation 10 page 226.
     REQUIRES_RUPTURE_PARAMETERS = set(('mag',))

     #: Required distance measure is hypocentral distance,
     #: see equation 1 page 834.
     REQUIRES_DISTANCES = set(('repi', 'azimuth'))

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
             a1ca = C['ua']
             a1cb = C['ub']
             a1cc = C['uc']
             a1cd = C['ud']
             a1ce = C['ue']
             a2ca = C['ia']
             a2cb = C['ib']
             a2cc = C['ic']
             a2cd = C['id']
             a2ce = C['ie']
         else:
             a1ca = C['a']
             a1cb = C['b']
             a1cc = C['c']
             a1cd = C['d']
             a1ce = C['e']
             a2ca = C['ma']
             a2cb = C['mb']
             a2cc = C['mc']
             a2cd = C['md']
             a2ce = C['me']

         # Compute the mean value (i.e. natural logarithm of ground motion)
         mag = rup.mag

         mean1 = a1ca + a1cb * mag + a1cc * np.log10(dists.repi + a1cd * np.exp(a1ce * mag))
         mean2 = a2ca + a2cb * mag + a2cc * np.log10(dists.repi + a2cd * np.exp(a2ce * mag))
         x = 10**mean1 * 10**mean1 * np.sin(np.radians(dists.azimuth)) * np.sin(np.radians(dists.azimuth))
         y = 10**mean2 * 10**mean2 * np.cos(np.radians(dists.azimuth)) * np.cos(np.radians(dists.azimuth))
         mean = 10**mean1 * 10**mean2 / np.sqrt(x+y)

         # Convert acceleration from gals into fraction of g
         if isinstance(imt, (PGA, SA)):
             mean = mean/g/100

         # Get the standard deviation
         stddevs = self._compute_std(C, stddev_types, len(dists.repi))
         # Return results
         return np.log(mean), stddevs

     def _compute_std(self, C, stddev_types, num_sites):
         """
         """
         return [np.ones(num_sites)*C['sigma']]

     #: Coefficient table for rock sites, see table 3 page 227.
     COEFFS = CoeffsTable(sa_damping=5, table="""\
         IMT        a      b       c      d      e     ua     ub      uc     ud     ue     ma     mb      mc     md     me     ia     ib         ic     id     ie   sigma
         PGA    1.791  0.720  -2.389  1.772  0.424  3.403  0.472  -2.389  1.772  0.424  0.983  0.713  -2.118  0.825  0.465  2.610  0.463  -2.118  0.825  0.465   0.236
         PGV   -0.547  0.840  -2.181  1.772  0.424  1.310  0.554  -2.181  1.772  0.424 -1.351  0.843  -1.945  0.825  0.465  0.569 0.549  -1.945  0.825  0.465   0.271
     """)

class YuEtAl2013Tibetan(YuEtAl2013):
     
     DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST
     COEFFS = CoeffsTable(sa_damping=5, table="""\
         IMT        a      b       c      d      e     ua     ub      uc     ud     ue     ma     mb      mc     md     me     ia     ib         ic     id     ie   sigma
         PGA    2.387  0.645  -2.416  2.647  0.366  3.807  0.411  -2.416  2.647  0.366  1.003  0.609  -1.854  0.612  0.457  2.457  0.388  -1.854  0.612  0.457   0.236
         PGV   -0.064  0.766  -2.205  2.647  0.366  1.714  0.491  -2.205  2.647  0.366  -1.301 0.741  -1.696  0.612  0.457  0.443 0.474  -1.696  0.612  0.457   0.271
     """)

class YuEtAl2013Eastern(YuEtAl2013):
    
     DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST
     COEFFS = CoeffsTable(sa_damping=5, table="""\
         IMT        a      b       c      d      e     ua     ub      uc     ud     ue     ma     mb      mc     md     me     ia     ib         ic     id     ie   sigma
         PGA    1.979  0.671  -2.315  2.088  0.399  3.533  0.432  -2.315  2.088  0.399  1.176  0.660  -2.004  0.944  0.447  2.753 0.418  -2.004  0.944  0.447   0.236
         PGV   -0.363  0.791  -2.103  2.088  0.399  1.437  0.513  -2.103  2.088  0.399 -1.147  0.788  -1.825  0.944  0.447  0.712  0.502  -1.825  0.944  0.447   0.271
     """)

class YuEtAl2013Stable(YuEtAl2013):

     DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL
     #: Coefficient table for rock sites, see table 3 page 227.
     COEFFS = CoeffsTable(sa_damping=5, table="""\
         IMT        a      b       c      d      e     ua     ub      uc     ud     ue     ma     mb      mc     md     me     ia     ib         ic     id     ie   sigma
         PGA    2.417  0.498  -2.079  2.802  0.295  3.706  0.298  -2.079  2.802  0.295  1.715  0.471  -1.723  1.295  0.331  2.690  0.321  -1.723  1.295  0.331   0.236
         PGV    0.093  0.621  -1.889  2.802  0.295  1.640  0.382  -1.889  2.802  0.295 -0.589  0.601  -1.559  1.295  0.331  0.671  0.407  -1.559  1.295  0.331   0.271
     """)




