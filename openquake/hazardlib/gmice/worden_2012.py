# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2016 GEM Foundation
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

import numpy as np
from scipy.constants import g

from openquake.hazardlib import const
from openquake.hazardlib.imt import MMI, PGA, PGV, SA
from openquake.hazardlib.gmice.base import GMICE
from openquake.hazardlib.gsim.base import CoeffsTable
#g = 9.81

class WordenEtAl2012(GMICE):
    """
    Implements the GMICE of Worden et al. 2012:
    Worden, C. B., Gerstenberger, M. C., Rhoades, D. A. and Wald, D. J., 2012,
    Probabilistic Relationship between Ground-Motion Parameters and Modified
    Mercalli Intensity in California, Bulletin of the Seismological Society
    of America, 102(1), 204 - 221

    Two alternative formulations are provided by the source. The implementation
    here does not consider magnitude or distance within the intensity
    conversion.
    """
    #: Supported tectonic region type is stable continental crust,
    #: given that the equations have been derived for central and eastern
    #: north America
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types are spectral acceleration,
    #: and peak ground acceleration
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        MMI,
        PGA,
        PGV,
        SA,
    ])

    #: Supported intensity measure component is the geometric mean of
    #two : horizontal components
    #:attr:`~openquake.hazardlib.const.IMC.AVERAGE_HORIZONTAL`,
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    #: Supported standard deviation type is only total.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: No site parameters required
    REQUIRES_SITES_PARAMETERS = set(())

    #: No rupture parameters required
    REQUIRES_RUPTURE_PARAMETERS = set(())

    #: No distance parameters required
    REQUIRES_DISTANCES = set(())

    #: Intensity range MMI 1 to 10
    INTENSITY_RANGE = (1.0, 10.0)

    #: Source Name
    SOURCE_NAME = "Worden et al. (2012)"

    #: Scale
    SCALE = 'scale_wgrw12.ps'
    
    def get_mean_intensity_and_stddevs(self, imls, sites, rup, dists, imt,
                                      stddev_types):
        """
        Implements equation 3 of Worden et al.
        """
        C = self.COEFFS[imt]
        C2 = self.COEFFS2[imt]
        # OpenQuake returns accelerations in terms of g
        if isinstance(imt, (PGA, SA)):
            # Convert from g to cm/s^2
            mmi = self.get_mean_mmi(C, C2, np.log10(imls * (100.0 * g)))
        else:
            mmi = self.get_mean_mmi(C, C2, np.log10(imls))
        stddevs = self.get_stddev_mmi(C, stddev_types, imls.shape)
        return mmi, stddevs

    def get_mean_mmi(self, C, C2, log_imls):
        """
        Return the expected macroseismic intensity from the ground motion input
        """
        mmi = C["c1"] + C["c2"] * log_imls
        idx = log_imls < C2["t1"]
        if np.any(idx):
            mmi[idx] = C2["c1"] + C2["c2"] * log_imls[idx]
        idx = log_imls >= C["t1"]
        if np.any(idx):
            mmi[idx] = C["c3"] + C["c4"] * log_imls[idx]
        return np.clip(mmi, 1.0, 10.0)
    
    def get_stddev_mmi(self, C, stddev_types, imls_shape):
        """
        Returns the total standard deviation in natural log units
        N. B. It is assumed the standard deviation is in common log units
        """
        stddevs = []
        for stddev_type in stddev_types:
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(np.log(10.0 ** (C["sigMMI"] +
                                               np.zeros(imls_shape))))
        return stddevs
    
    def get_mean_gm_and_stddevs(self, imls, sites, rup, dists, imt,
                                stddev_types):
        """
        Convert from macroseismic intensity to engineering ground motion
        parameters
        """
        C = self.COEFFS[imt]
        C2 = self.COEFFS2[imt]
        g_m = self.get_mean_gm(C, C2, imls)
        # If IMT is PGA or SA
        if isinstance(imt, (PGA, SA)):
            # Convert units from cm/s^2 to g
            g_m = np.log(g_m / (100.0 * g))
        else:
            g_m = np.log(g_m)
        stddevs = self.get_stddev_gm(C, stddev_types, imls.shape)
        return g_m, stddevs
    
    def get_mean_gm(self, C, C2, mmi):
        """
        Returns the expected ground motion from the macroseismic intensity input
        """
        logy = (mmi - C["c1"]) / C["c2"]
        idx = mmi < C2["t2"]
        if np.any(idx):
            logy[idx] = (mmi[idx] - C2["c1"]) / C2["c2"]
        idx = mmi >= C["t2"]
        if np.any(idx):
            logy[idx] = (mmi[idx] - C['c3']) / C["c4"]
        return 10.0 ** logy

    def get_stddev_gm(self, C, stddev_types, imls_shape):
        """
        """
        stddevs = []
        for stddev_type in stddev_types:
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(np.log(10.0 ** (C["siglogy"] +
                                               np.zeros(imls_shape))))
        return stddevs

    # This will interpolate coefficients - which may not be the author's
    # original intention.
    # In the original shakemap implementation the standard deviations for the
    # version considering magnitude and distance are applied to both GMICEs
    COEFFS = CoeffsTable(sa_damping=5, table="""
    imt     c1     c2     c3    c4    t1    t2  sigMMI  siglogy
    pga   1.78   1.55  -1.60  3.70  1.57  4.22    0.66     0.35
    pgv   3.78   1.47   2.89  3.16  0.53  4.56    0.63     0.38
    0.3   1.26   1.69  -4.15  4.14  2.21  4.99    0.82     0.44
    1.0   2.50   1.51   0.20  2.90  1.65  4.98    0.75     0.47
    3.0   3.81   1.17   1.99  3.01  0.99  4.96    0.89     0.64
    """)

    COEFFS2 = CoeffsTable(sa_damping=5, table="""
    imt    c1    c2     t1   t2
    pga  1.71  2.08   0.14  2.0
    pgv  4.62  2.17  -1.21  2.0
    0.3  1.15  1.92   0.44  2.0
    1.0  2.71  2.17  -0.33  2.0
    3.0  7.35  3.45  -1.55  2.0
    """)


class WordenEtAl2012MagDist(WordenEtAl2012):
    """
    """
    #: rupture magnitude required
    REQUIRES_RUPTURE_PARAMETERS = set(("mag",))

    #: rupture distance required
    REQUIRES_DISTANCES = set(("rrup",))
    
    def get_mean_intensity_and_stddevs(self, imls, sites, rup, dists, imt,
                                       stddev_types):
        """
        Implements equation 3 of Worden et al.
        """
        C = self.COEFFS[imt]
        C2 = self.COEFFS2[imt]
        # Clip the magnitudes between 3.0 and 7.3
        mag = np.clip(rup.mag, 3.0, 7.3)
        # Clip the distances between 10 km and 300 km
        rrup = np.clip(dists.rrup, 10., 300.)
        if isinstance(imt, (PGA, SA)):
            # Convert from g to cm/s^2
            mmi = self.get_mean_mmi(C, C2, np.log10(imls * (100.0 * g)),
                                    mag, rrup)
        else:
            mmi = self.get_mean_mmi(C, C2, np.log10(imls), mag, rrup)
        stddevs = self.get_stddev_mmi(C, stddev_types, imls.shape)
        return mmi, stddevs

    def get_mean_mmi(self, C, C2, log_imls, mag, rrup):
        """
        """
        mmi = C["c5"] + C["c6"] * np.log10(rrup) + C["c7"] * mag
        idx = log_imls <= C2["t1"]
        if np.any(idx):
            mmi[idx] += (C2["c1"] + (C2["c2"] * log_imls[idx]))
        idx = np.logical_and(log_imls > C2["t1"], log_imls <= C["t1"])
        if np.any(idx):
            mmi[idx] += (C["c1"] + (C["c2"] * log_imls[idx]))
        idx = log_imls > C["t1"]
        if np.any(idx):
            mmi[idx] += (C["c3"] + C["c4"] * log_imls[idx])
        return np.clip(mmi, 1.0, 10.0)
    
    def get_mean_gm_and_stddevs(self, imls, sites, rup, dists, imt,
                                stddev_types):
        """
        Convert from macroseismic intensity to engineering ground motion
        parameters
        """
        C = self.COEFFS[imt]
        C2 = self.COEFFS2[imt]
        # Clip the magnitudes between 3.0 and 7.3
        mag = np.clip(rup.mag, 3.0, 7.3)
        # Clip the distances between 10 km and 300 km
        rrup = np.clip(dists.rrup, 10., 300.)
        g_m = self.get_mean_gm(C, C2,imls, mag, rrup)
        # If IMT is PGA or SA
        if isinstance(imt, (PGA, SA)):
            # Convert units from cm/s^2 to g
            g_m = np.log(g_m / (100.0 * g))
        else:
            g_m = np.log(g_m)
        stddevs = self.get_stddev_gm(C, stddev_types, imls.shape)
        return g_m, stddevs
    
    def get_mean_gm(self, C, C2, mmi, mag, rrup):
        """
        Returns the expected ground motion from the macroseismic intensity
        """
        logy = mmi - (C["c5"] + C["c6"] * np.log10(rrup) + C["c7"] * mag)
        
        idx = mmi < C2["t2"]
        if np.any(idx):
            logy[idx] = (logy[idx] - C2["c1"]) / C2["c2"]
        idx = np.logical_and(mmi >= C2["t2"], mmi < C["t2"])
        if np.any(idx):
            logy[idx] = (logy[idx] - C["c1"]) / C["c2"]
        idx = mmi >= C["t2"]
        if np.any(idx):
            logy[idx] = (logy[idx] - C["c3"]) / C["c4"]
        
        return 10.0 ** logy
    
    COEFFS = CoeffsTable(sa_damping=5, table="""
    imt     c1     c2     c3    c4     c5     c6     c7    t1    t2  sigMMI  siglogy
    pga   1.78   1.55  -1.60  3.70  -0.91   1.02  -0.17  1.57  4.22    0.66     0.35
    pgv   3.78   1.47   2.89  3.16   0.90   0.00  -0.18  0.53  4.56    0.63     0.38
    0.3   1.26   1.69  -4.15  4.14  -1.05   0.60   0.00  2.21  4.99    0.82     0.44
    1.0   2.50   1.51   0.20  2.90   2.27  -0.49  -0.29  1.65  4.98    0.75     0.47
    3.0   3.81   1.17   1.99  3.01   1.91  -0.57  -0.21  0.99  4.96    0.89     0.64
    """)
