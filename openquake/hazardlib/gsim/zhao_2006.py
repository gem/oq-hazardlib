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
Module exports :class:`ZhaoEtAl2006Asc`, class:`ZhaoEtAl2006SInter`, and
class:`ZhaoEtAl2006SSlab`.
"""
from __future__ import division

import numpy as np
# standard acceleration of gravity in m/s**2
from scipy.constants import g

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA


class ZhaoEtAl2006Asc(GMPE):

    """
    Implements GMPE developed by John X. Zhao et al. and published as
    "Attenuation Relations of Strong Ground Motion in Japan Using Site
    Classification Based on Predominant Period" (2006, Bulletin of the
    Seismological Society of America, Volume 96, No. 3, pages 898-913).
    This class implements the equations for 'Active Shallow Crust'
    (that's why the class name ends with 'Asc').
    """
    #: Supported tectonic region type is active shallow crust, this means
    #: that factors SI, SS and SSL are assumed 0 in equation 1, p. 901.
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types are spectral acceleration,
    #: and peak ground acceleration, see paragraph 'Development of Base Model'
    #: p. 901.
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
        SA
    ])

    #: Supported intensity measure component is geometric mean
    #: of two horizontal components :
    #: attr:`~openquake.hazardlib.const.IMC.AVERAGE_HORIZONTAL`, see paragraph
    #: 'Development of Base Model', p. 901.
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    #: Supported standard deviation types are inter-event, intra-event
    #: and total, see equation 3, p. 902.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL,
        const.StdDev.INTER_EVENT,
        const.StdDev.INTRA_EVENT
    ])

    #: Required site parameters is Vs30.
    #: See table 2, p. 901.
    REQUIRES_SITES_PARAMETERS = set(('vs30', ))

    #: Required rupture parameters are magnitude, rake, and focal depth.
    #: See paragraph 'Development of Base Model', p. 901.
    REQUIRES_RUPTURE_PARAMETERS = set(('mag', 'rake', 'hypo_depth'))

    #: Required distance measure is Rrup.
    #: See paragraph 'Development of Base Model', p. 902.
    REQUIRES_DISTANCES = set(('rrup', ))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # extracting dictionary of coefficients specific to required
        # intensity measure type.
        C = self.COEFFS_ASC[imt]

        # mean value as given by equation 1, p. 901, without considering the
        # interface and intraslab terms (that is SI, SS, SSL = 0) and the
        # inter and intra event terms, plus the magnitude-squared term
        # correction factor (equation 5 p. 909).
        mean = self._compute_magnitude_term(C, rup.mag) +\
            self._compute_distance_term(C, rup.mag, dists.rrup) +\
            self._compute_focal_depth_term(C, rup.hypo_depth) +\
            self._compute_faulting_style_term(C, rup.rake) +\
            self._compute_site_class_term(C, sites.vs30) +\
            self._compute_magnitude_squared_term(P=0.0, M=6.3, Q=C['QC'],
                                                 W=C['WC'], mag=rup.mag)

        # convert from cm/s**2 to g
        mean = np.log(np.exp(mean) * 1e-2 / g)

        stddevs = self._get_stddevs(C['sigma'], C['tauC'], stddev_types,
                                    num_sites=len(sites.vs30))

        return mean, stddevs

    def _get_stddevs(self, sigma, tau, stddev_types, num_sites):
        """
        Return standard deviations as defined in equation 3 p. 902.
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                sigma_t = np.sqrt(sigma ** 2 + tau ** 2)
                stddevs.append(sigma_t + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(sigma + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(tau + np.zeros(num_sites))
        return stddevs

    def _compute_magnitude_term(self, C, mag):
        """
        Compute first term in equation 1, p. 901.
        """
        return C['a'] * mag

    def _compute_distance_term(self, C, mag, rrup):
        """
        Compute second and third terms in equation 1, p. 901.
        """
        term1 = C['b'] * rrup
        term2 = - np.log(rrup + C['c'] * np.exp(C['d'] * mag))

        return term1 + term2

    def _compute_focal_depth_term(self, C, hypo_depth):
        """
        Compute fourth term in equation 1, p. 901.
        """
        # p. 901. "(i.e, depth is capped at 125 km)".
        focal_depth = hypo_depth
        if focal_depth > 125.0:
            focal_depth = 125.0

        # p. 902. "We used the value of 15 km for the
        # depth coefficient hc ...".
        hc = 15.0

        # p. 901. "When h is larger than hc, the depth terms takes
        # effect ...". The next sentence specifies h>=hc.
        return float(focal_depth >= hc) * C['e'] * (focal_depth - hc)

    def _compute_faulting_style_term(self, C, rake):
        """
        Compute fifth term in equation 1, p. 901.
        """
        # p. 900. "The differentiation in focal mechanism was
        # based on a rake angle criterion, with a rake of +/- 45
        # as demarcation between dip-slip and strike-slip."
        return float(rake > 45.0 and rake < 135.0) * C['FR']

    def _compute_site_class_term(self, C, vs30):
        """
        Compute nine-th term in equation 1, p. 901.
        """
        # map vs30 value to site class, see table 2, p. 901.
        site_term = np.zeros(len(vs30))

        # hard rock
        site_term[vs30 > 1100.0] = C['CH']

        # rock
        site_term[(vs30 > 600) & (vs30 <= 1100)] = C['C1']

        # hard soil
        site_term[(vs30 > 300) & (vs30 <= 600)] = C['C2']

        # medium soil
        site_term[(vs30 > 200) & (vs30 <= 300)] = C['C3']

        # soft soil
        site_term[vs30 <= 200] = C['C4']

        return site_term

    def _compute_magnitude_squared_term(self, P, M, Q, W, mag):
        """
        Compute magnitude squared term, equation 5, p. 909.
        """
        return P * (mag - M) + Q * (mag - M) ** 2 + W

    #: Coefficient table obtained by joining table 4 (except columns for
    #: SI, SS, SSL), table 5 (both at p. 903) and table 6 (only columns for
    #: QC WC TauC), p. 907.
    COEFFS_ASC = CoeffsTable(sa_damping=5, table="""\
    IMT    a     b         c       d      e        FR     CH     C1     C2     C3     C4     sigma   QC      WC      tauC
    pga    1.101 -0.00564  0.0055  1.080  0.01412  0.251  0.293  1.111  1.344  1.355  1.420  0.604   0.0     0.0     0.303
    0.05   1.076 -0.00671  0.0075  1.060  0.01463  0.251  0.939  1.684  1.793  1.747  1.814  0.640   0.0     0.0     0.326
    0.10   1.118 -0.00787  0.0090  1.083  0.01423  0.240  1.499  2.061  2.135  2.031  2.082  0.694   0.0     0.0     0.342
    0.15   1.134 -0.00722  0.0100  1.053  0.01509  0.251  1.462  1.916  2.168  2.052  2.113  0.702   0.0     0.0     0.331
    0.20   1.147 -0.00659  0.0120  1.014  0.01462  0.260  1.280  1.669  2.085  2.001  2.030  0.692   0.0     0.0     0.312
    0.25   1.149 -0.00590  0.0140  0.966  0.01459  0.269  1.121  1.468  1.942  1.941  1.937  0.682   0.0     0.0     0.298
    0.30   1.163 -0.00520  0.0150  0.934  0.01458  0.259  0.852  1.172  1.683  1.808  1.770  0.670   0.0     0.0     0.300
    0.40   1.200 -0.00422  0.0100  0.959  0.01257  0.248  0.365  0.655  1.127  1.482  1.397  0.659   0.0     0.0     0.346
    0.50   1.250 -0.00338  0.0060  1.008  0.01114  0.247 -0.207  0.071  0.515  0.934  0.955  0.653  -0.0126  0.0116  0.338
    0.60   1.293 -0.00282  0.0030  1.088  0.01019  0.233 -0.705 -0.429 -0.003  0.394  0.559  0.653  -0.0329  0.0202  0.349
    0.70   1.336 -0.00258  0.0025  1.084  0.00979  0.220 -1.144 -0.866 -0.449 -0.111  0.188  0.652  -0.0501  0.0274  0.351
    0.80   1.386 -0.00242  0.0022  1.088  0.00944  0.232 -1.609 -1.325 -0.928 -0.620 -0.246  0.647  -0.0650  0.0336  0.356
    0.90   1.433 -0.00232  0.0020  1.109  0.00972  0.220 -2.023 -1.732 -1.349 -1.066 -0.643  0.653  -0.0781  0.0391  0.348
    1.00   1.479 -0.00220  0.0020  1.115  0.01005  0.211 -2.451 -2.152 -1.776 -1.523 -1.084  0.657  -0.0899  0.0440  0.338
    1.25   1.551 -0.00207  0.0020  1.083  0.01003  0.251 -3.243 -2.923 -2.542 -2.327 -1.936  0.660  -0.1148  0.0545  0.313
    1.50   1.621 -0.00224  0.0020  1.091  0.00928  0.248 -3.888 -3.548 -3.169 -2.979 -2.661  0.664  -0.1351  0.0630  0.306
    2.00   1.694 -0.00201  0.0025  1.055  0.00833  0.263 -4.783 -4.410 -4.039 -3.871 -3.640  0.669  -0.1672  0.0764  0.283
    2.50   1.748 -0.00187  0.0028  1.052  0.00776  0.262 -5.444 -5.049 -4.698 -4.496 -4.341  0.671  -0.1921  0.0869  0.287
    3.00   1.759 -0.00147  0.0032  1.025  0.00644  0.307 -5.839 -5.431 -5.089 -4.893 -4.758  0.667  -0.2124  0.0954  0.278
    4.00   1.826 -0.00195  0.0040  1.044  0.00590  0.353 -6.598 -6.181 -5.882 -5.698 -5.588  0.647  -0.2445  0.1088  0.273
    5.00   1.825 -0.00237  0.0050  1.065  0.00510  0.248 -6.752 -6.347 -6.051 -5.873 -5.798  0.643  -0.2694  0.1193  0.275
    """)


class ZhaoEtAl2006SInter(ZhaoEtAl2006Asc):

    """
    Implements GMPE developed by John X. Zhao et al and published as
    "Attenuation Relations of Strong Ground Motion in Japan Using Site
    Classification Based on Predominant Period" (2006, Bulletin of the
    Seismological Society of America, Volume 96, No. 3, pages
    898-913). This class implements the equations for 'Subduction
    Interface' (that's why the class name ends with 'SInter'). This
    class extends the
    :class:`openquake.hazardlib.gsim.zhao_2006.ZhaoEtAl2006Asc`
    because the equation for subduction interface is obtained from the
    equation for active shallow crust, by removing the faulting style
    term and adding a subduction interface term.
    """
    #: Supported tectonic region type is subduction interface, this means
    #: that factors FR, SS and SSL are assumed 0 in equation 1, p. 901.
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # extracting dictionary of coefficients specific to required
        # intensity measure type.
        C = self.COEFFS_ASC[imt]
        C_SINTER = self.COEFFS_SINTER[imt]

        # mean value as given by equation 1, p. 901, without considering the
        # faulting style and intraslab terms (that is FR, SS, SSL = 0) and the
        # inter and intra event terms, plus the magnitude-squared term
        # correction factor (equation 5 p. 909)
        mean = self._compute_magnitude_term(C, rup.mag) +\
            self._compute_distance_term(C, rup.mag, dists.rrup) +\
            self._compute_focal_depth_term(C, rup.hypo_depth) +\
            self._compute_site_class_term(C, sites.vs30) + \
            self._compute_magnitude_squared_term(P=0.0, M=6.3,
                                                 Q=C_SINTER['QI'],
                                                 W=C_SINTER['WI'],
                                                 mag=rup.mag) +\
            C_SINTER['SI']

        # convert from cm/s**2 to g
        mean = np.log(np.exp(mean) * 1e-2 / g)

        stddevs = self._get_stddevs(C['sigma'], C_SINTER['tauI'], stddev_types,
                                    num_sites=len(sites.vs30))

        return mean, stddevs

    #: Coefficient table containing subduction interface coefficients,
    #: taken from table 4, p. 903 (only column SI), and table 6, p. 907
    #: (only columns QI, WI, TauI)
    COEFFS_SINTER = CoeffsTable(sa_damping=5, table="""\
        IMT    SI     QI      WI      tauI
        pga    0.000  0.0     0.0     0.308
        0.05   0.000  0.0     0.0     0.343
        0.10   0.000  0.0     0.0     0.403
        0.15   0.000 -0.0138  0.0286  0.367
        0.20   0.000 -0.0256  0.0352  0.328
        0.25   0.000 -0.0348  0.0403  0.289
        0.30   0.000 -0.0423  0.0445  0.280
        0.40  -0.041 -0.0541  0.0511  0.271
        0.50  -0.053 -0.0632  0.0562  0.277
        0.60  -0.103 -0.0707  0.0604  0.296
        0.70  -0.146 -0.0771  0.0639  0.313
        0.80  -0.164 -0.0825  0.0670  0.329
        0.90  -0.206 -0.0874  0.0697  0.324
        1.00  -0.239 -0.0917  0.0721  0.328
        1.25  -0.256 -0.1009  0.0772  0.339
        1.50  -0.306 -0.1083  0.0814  0.352
        2.00  -0.321 -0.1202  0.0880  0.360
        2.50  -0.337 -0.1293  0.0931  0.356
        3.00  -0.331 -0.1368  0.0972  0.338
        4.00  -0.390 -0.1486  0.1038  0.307
        5.00  -0.498 -0.1578  0.1090  0.272
        """)


class ZhaoEtAl2006SSlab(ZhaoEtAl2006Asc):

    """
    Implements GMPE developed by John X. Zhao et al and published as
    "Attenuation Relations of Strong Ground Motion in Japan Using Site
    Classification Based on Predominant Period" (2006, Bulletin of the
    Seismological Society of America, Volume 96, No. 3, pages
    898-913). This class implements the equations for 'Subduction
    Slab'. (that's why the class name ends with 'SSlab'). This class
    extends the
    :class:`openquake.hazardlib.gsim.zhao_2006.ZhaoEtAl2006Asc`
    because the equation for subduction slab is obtained from the
    equation for active shallow crust, by removing the faulting style
    term and adding subduction slab terms.
    """
    #: Supported tectonic region type is subduction interface, this means
    #: that factors FR, SS and SSL are assumed 0 in equation 1, p. 901.
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTRASLAB

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # extracting dictionary of coefficients specific to required
        # intensity measure type.
        C = self.COEFFS_ASC[imt]
        C_SSLAB = self.COEFFS_SSLAB[imt]

        # to avoid singularity at 0.0 (in the calculation of the
        # slab correction term), replace 0 values with 0.1
        d = dists.rrup
        d[d == 0.0] = 0.1

        # mean value as given by equation 1, p. 901, without considering the
        # faulting style and intraslab terms (that is FR, SS, SSL = 0) and the
        # inter and intra event terms, plus the magnitude-squared term
        # correction factor (equation 5 p. 909)
        mean = self._compute_magnitude_term(C, rup.mag) +\
            self._compute_distance_term(C, rup.mag, d) +\
            self._compute_focal_depth_term(C, rup.hypo_depth) +\
            self._compute_site_class_term(C, sites.vs30) +\
            self._compute_magnitude_squared_term(P=C_SSLAB['PS'], M=6.5,
                                                 Q=C_SSLAB['QS'],
                                                 W=C_SSLAB['WS'],
                                                 mag=rup.mag) +\
            C_SSLAB['SS'] + self._compute_slab_correction_term(C_SSLAB, d)

        # convert from cm/s**2 to g
        mean = np.log(np.exp(mean) * 1e-2 / g)

        stddevs = self._get_stddevs(C['sigma'], C_SSLAB['tauS'], stddev_types,
                                    num_sites=len(sites.vs30))

        return mean, stddevs

    def _compute_slab_correction_term(self, C, rrup):
        """
        Compute path modification term for slab events, that is
        the 8-th term in equation 1, p. 901.
        """
        slab_term = C['SSL'] * np.log(rrup)

        return slab_term

    #: Coefficient table containing subduction slab coefficients taken from
    #: table 4, p. 903 (only columns for SS and SSL), and table 6, p. 907
    #: (only columns for PS, QS, WS, TauS)
    COEFFS_SSLAB = CoeffsTable(sa_damping=5, table="""\
        IMT    SS     SSL     PS      QS       WS      tauS
        pga    2.607 -0.528   0.1392  0.1584  -0.0529  0.321
        0.05   2.764 -0.551   0.1636  0.1932  -0.0841  0.378
        0.10   2.156 -0.420   0.1690  0.2057  -0.0877  0.420
        0.15   2.161 -0.431   0.1669  0.1984  -0.0773  0.372
        0.20   1.901 -0.372   0.1631  0.1856  -0.0644  0.324
        0.25   1.814 -0.360   0.1588  0.1714  -0.0515  0.294
        0.30   2.181 -0.450   0.1544  0.1573  -0.0395  0.284
        0.40   2.432 -0.506   0.1460  0.1309  -0.0183  0.278
        0.50   2.629 -0.554   0.1381  0.1078  -0.0008  0.272
        0.60   2.702 -0.575   0.1307  0.0878   0.0136  0.285
        0.70   2.654 -0.572   0.1239  0.0705   0.0254  0.290
        0.80   2.480 -0.540   0.1176  0.0556   0.0352  0.299
        0.90   2.332 -0.522   0.1116  0.0426   0.0432  0.289
        1.00   2.233 -0.509   0.1060  0.0314   0.0498  0.286
        1.25   2.029 -0.469   0.0933  0.0093   0.0612  0.277
        1.50   1.589 -0.379   0.0821 -0.0062   0.0674  0.282
        2.00   0.966 -0.248   0.0628 -0.0235   0.0692  0.300
        2.50   0.789 -0.221   0.0465 -0.0287   0.0622  0.292
        3.00   1.037 -0.263   0.0322 -0.0261   0.0496  0.274
        4.00   0.561 -0.169   0.0083 -0.0065   0.0150  0.281
        5.00   0.225 -0.120  -0.0117  0.0246  -0.0268  0.296
        """)


class ZhaoEtAl2006AscSWISS05(ZhaoEtAl2006Asc):

    """
    --------------------------------------------------------------------
    This class implments an extension of the ZhaoEtAl (2006Asc) model,
    adjusted to be used for the new Swiss Hazard Model [2014].
    1) kappa value
       K-adjustments corresponding to model 05 (lower model) - as prepared by Ben Edwards
       K-value for PGA were not provided but infered from SA[]0.01s]
       the model considers a fixed value of vs30=1100m/s
    2) small-magnitude correction
    3) single station sigma - mean inter-event adjustment
    4) single station sigma - inter-event magnitude/distance dependent
    ------------------------------------------------------------------------
    Disclaimer: these equations are modified to be used for the
    new Swiss Seismic Hazard Model [2014].
    The use of these models is the soly responsability of the hazard modeler.
    --------------------------------------------------------------------
    Model implmented by laurentiu.danciu@sed.ethz.ch
    --------------------------------------------------------------------
    """

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):

        C_ADJ = self.COEFFS_FS_ROCK[imt]
        C1_rrup = self._compute_C1_term(C_ADJ, imt, dists)
        phi_ss = self._compute_phi_ss(C_ADJ, rup, C1_rrup, imt)

        mean, stddevs = super(ZhaoEtAl2006AscSWISS05, self).\
            get_mean_and_stddevs(sites, rup, dists, imt, stddev_types)

        #: apply k-correction corresponding to the lower model [01]
        mean_corr = np.exp(
            mean) * C_ADJ['k_adj'] * self._compute_small_mag_correction_term(C_ADJ, rup.mag, dists.rrup)
        mean = np.log(mean_corr)

        std_corr = self._get_corr_stddevs(
            self.COEFFS_ASC[imt], stddev_types, len(sites.vs30), phi_ss)
        stddevs = np.array(std_corr)

        return mean, stddevs

    def _get_corr_stddevs(self, C, stddev_types, num_sites, phi_ss):
        """
        Return standard deviations adjusted for single station sigma
        as proposed to be used in the new Swiss Hazard Model [2014].
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(
                    np.sqrt(
                        C['tauC'] ** 2 +
                        phi_ss *
                        phi_ss) +
                    np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(C['sigma'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(C['tauC'] + np.zeros(num_sites))
        return stddevs

    def _compute_small_mag_correction_term(self, C, mag, rrup):
        if mag >= 3.00 and mag < 5.5:
            return 1 / np.exp(((5.50 - mag) / C['a1']) ** C['a2'] * (
                C['b1'] + C['b2'] * np.log(np.maximum(np.minimum(rrup, C['Rm']), 10) / 20)))
        elif mag >= 5.50:
            return 1
        else:
            return 1

    def _compute_C1_term(self, C, imt, dists):
        """
        Return C1 coeffs as function of Rrup as proposed by Rodriguez-Marek et al (2013)
        The C1 coeff are used to compute the single station sigma
        """
        C1_rrup = 0.0
        if (dists.rrup < C['Rc11']).any():
            C1_rrup = C['phi_11']
        elif ((dists.rrup >= C['Rc11']).any()
              and (dists.rrup <= C['Rc21']).any()):
            C1_rrup = C['phi_11'] + (C['phi_21'] - C['phi_11']) * \
                ((dists.rrup - C['Rc11']) / (C['Rc21'] - C['Rc11']))
        elif (dists.rrup > C['Rc21']).any():
            C1_rrup = C['phi_21']
        return C1_rrup

    def _compute_phi_ss(self, C, rup, C1_rrup, imt):
        """
        Return C1 coeffs as function of Rrup as proposed by Rodriguez-Marek et al (2013)
        The C1 coeff are used to compute the single station sigma
        """
        phi_ss = 0

        if rup.mag < C['Mc1']:
            phi_ss = C1_rrup
        elif rup.mag >= C['Mc1'] and rup.mag <= C['Mc2']:
            phi_ss = C1_rrup + \
                (C['C2'] - C1_rrup) * \
                ((rup.mag - C['Mc1']) / (C['Mc2'] - C['Mc1']))
        elif rup.mag > C['Mc2']:
            phi_ss = C['C2']
        return phi_ss

    COEFFS_FS_ROCK = CoeffsTable(sa_damping=5, table="""\
       IMT    k_adj           a1              a2              b1              b2              Rm             phi_11     phi_21  C2      Mc1     Mc2     Rc11        Rc21
       pga    0.893272000     1.642085E+00    1.833001E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.58000    0.47000 0.35000 5.00000 7.00000 16.00000    36.00000
       0.05   0.960326204     1.590620E+00    1.617332E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.55204    0.44903 0.40592 5.00000 7.00000 16.00000    36.00000
       0.10   0.883646750     1.478437E+00    1.608714E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.54000    0.44000 0.43000 5.00000 7.00000 16.00000    36.00000
       0.15   0.848710124     1.488927E+00    1.732820E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.58095    0.47510 0.40075 5.00000 7.00000 16.00000    36.00000
       0.20   0.829169096     1.496370E+00    1.826639E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.61000    0.50000 0.38000 5.00000 7.00000 16.00000    36.00000
       0.25   0.819497568     1.473785E+00    1.854356E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.62651    0.50000 0.37450 5.00000 7.00000 16.00000    36.00000
       0.30   0.814346585     1.455332E+00    1.877313E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.64000    0.50000 0.37000 5.00000 7.00000 16.00000    36.00000
       0.40   0.814243546     1.428722E+00    1.886767E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.61747    0.48874 0.37000 5.00000 7.00000 16.00000    36.00000
       0.50   0.812705538     1.405794E+00    1.841480E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.60000    0.48000 0.37000 5.00000 7.00000 16.00000    36.00000
       0.60   0.815801037     1.387061E+00    1.805286E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.58422    0.47211 0.37789 5.00000 7.00000 16.00000    36.00000
       0.70   0.820928683     1.371222E+00    1.775240E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.57087    0.46544 0.38456 5.00000 7.00000 16.00000    36.00000
       0.80   0.827376206     1.357502E+00    1.749618E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.55932    0.45966 0.39034 5.00000 7.00000 16.00000    36.00000
       0.90   0.831408127     1.345400E+00    1.727324E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.54912    0.45456 0.39544 5.00000 7.00000 16.00000    36.00000
       1.00   0.840799961     1.334574E+00    1.707623E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.54000    0.45000 0.40000 5.00000 7.00000 16.00000    36.00000
       1.25   0.855226821     1.168218E+00    1.437403E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.53797    0.43984 0.40000 5.00000 7.00000 16.00000    36.00000
       1.50   0.871104233     1.032296E+00    1.248681E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.53631    0.43155 0.40000 5.00000 7.00000 16.00000    36.00000
       2.00   0.891458727     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.53369    0.41845 0.40000 5.00000 7.00000 16.00000    36.00000
       2.50   0.903858490     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.53166    0.40830 0.40000 5.00000 7.00000 16.00000    36.00000
       3.00   0.913991280     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.53000    0.40000 0.40000 5.00000 7.00000 16.00000    36.00000
       4.00   0.913352473     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.53000    0.40000 0.40000 5.00000 7.00000 16.00000    36.00000
       5.00   0.912857283     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09   0.53000    0.40000 0.40000 5.00000 7.00000 16.00000    36.00000
       """)


class ZhaoEtAl2006AscSWISS03(ZhaoEtAl2006Asc):

    """
    --------------------------------------------------------------------
    This class implments an extension of the ZhaoEtAl (2006Asc) model,
    adjusted to be used for the new Swiss Hazard Model [2014].
    1) kappa value
       K-adjustments corresponding to model 03(mid model) - as prepared by Ben Edwards
       K-value for PGA were not provided but infered from SA[]0.01s]
       the model considers a fixed value of vs30=1100m/s
    2) small-magnitude correction
    3) single station sigma - mean inter-event adjustment
    4) single station sigma - inter-event magnitude/distance dependent
    ------------------------------------------------------------------------
    Disclaimer: these equations are modified to be used for the
    new Swiss Seismic Hazard Model [2014].
    The use of these models is the soly responsability of the hazard modeler.
    --------------------------------------------------------------------
    Model implmented by laurentiu.danciu@sed.ethz.ch
    --------------------------------------------------------------------
    """

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):

        C_ADJ = self.COEFFS_FS_ROCK[imt]
        C1_rrup = self._compute_C1_term(C_ADJ, imt, dists)
        phi_ss = self._compute_phi_ss(C_ADJ, rup, C1_rrup, imt)

        mean, stddevs = super(ZhaoEtAl2006AscSWISS03, self).\
            get_mean_and_stddevs(sites, rup, dists, imt, stddev_types)

        #: apply k-correction corresponding to the lower model [01]
        mean_corr = np.exp(
            mean) * C_ADJ['k_adj'] * self._compute_small_mag_correction_term(C_ADJ, rup.mag, dists.rrup)
        mean = np.log(mean_corr)

        std_corr = self._get_corr_stddevs(
            self.COEFFS_ASC[imt], stddev_types, len(sites.vs30), phi_ss)
        stddevs = np.array(std_corr)

        return mean, stddevs

    def _get_corr_stddevs(self, C, stddev_types, num_sites, phi_ss):
        """
        Return standard deviations adjusted for single station sigma
        as proposed to be used in the new Swiss Hazard Model [2014].
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(
                    np.sqrt(
                        C['tauC'] ** 2 +
                        phi_ss ** 2) +
                    np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(C['sigma'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(C['tauC'] + np.zeros(num_sites))
        return stddevs

    def _compute_small_mag_correction_term(self, C, mag, rrup):
        if mag >= 3.00 and mag < 5.5:
            return 1 / np.exp(((5.50 - mag) / C['a1']) ** C['a2'] * (
                C['b1'] + C['b2'] * np.log(np.maximum(np.minimum(rrup, C['Rm']), 10) / 20)))
        elif mag >= 5.50:
            return 1
        else:
            return 1

    def _compute_C1_term(self, C, imt, dists):
        """
        Return C1 coeffs as function of Rrup as proposed by Rodriguez-Marek et al (2013)
        The C1 coeff are used to compute the single station sigma
        """
        C1_rrup = 0.0
        if (dists.rrup < C['Rc11']).any():
            C1_rrup = C['phi_11']
        elif ((dists.rrup >= C['Rc11']).any()
                and (dists.rrup <= C['Rc21']).any()):
            C1_rrup = C['phi_11'] + (C['phi_21'] - C['phi_11']) * \
                ((dists.rrup - C['Rc11']) / (C['Rc21'] - C['Rc11']))
        elif (dists.rrup > C['Rc21']).any():
            C1_rrup = C['phi_21']
        return C1_rrup

    def _compute_phi_ss(self, C, rup, C1_rrup, imt):
        """
        Return C1 coeffs as function of Rrup as proposed by Rodriguez-Marek et al (2013)
        The C1 coeff are used to compute the single station sigma
        """
        phi_ss = 0

        if rup.mag < C['Mc1']:
            phi_ss = C1_rrup
        elif rup.mag >= C['Mc1'] and rup.mag <= C['Mc2']:
            phi_ss = C1_rrup + \
                (C['C2'] - C1_rrup) * \
                ((rup.mag - C['Mc1']) / (C['Mc2'] - C['Mc1']))
        elif rup.mag > C['Mc2']:
            phi_ss = C['C2']
        return phi_ss

    COEFFS_FS_ROCK = CoeffsTable(sa_damping=5, table="""\
       IMT    k_adj           a1              a2              b1              b2              Rm              phi_11      phi_21      C2          Mc1    Mc2    Rc11    Rc21
       pga    1.037040000     1.642085E+00    1.833001E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.58        0.47        0.35        5      7      16      36
       0.05   1.152476093     1.590620E+00    1.617332E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.5520412   0.4490309   0.4059176   5      7      16      36
       0.10   0.995583662     1.478437E+00    1.608714E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.54        0.44        0.43        5      7      16      36
       0.15   0.948713303     1.488927E+00    1.732820E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.580947375 0.47509775  0.400751875 5      7      16      36
       0.20   0.936827687     1.496370E+00    1.826639E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.61        0.5         0.38        5      7      16      36
       0.25   0.941001497     1.473785E+00    1.854356E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.626510191 0.5         0.374496603 5      7      16      36
       0.30   0.951517574     1.455332E+00    1.877313E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.64        0.5         0.37        5      7      16      36
       0.40   0.980951997     1.428722E+00    1.886767E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.617473168 0.488736584 0.37        5      7      16      36
       0.50   0.999448607     1.405794E+00    1.841480E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.6         0.48        0.37        5      7      16      36
       0.60   1.013777169     1.387061E+00    1.805286E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.584217936 0.472108968 0.377891032 5      7      16      36
       0.70   1.022460327     1.371222E+00    1.775240E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.57087439  0.465437195 0.384562805 5      7      16      36
       0.80   1.026784122     1.357502E+00    1.749618E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.559315686 0.459657843 0.390342157 5      7      16      36
       0.90   1.024261640     1.345400E+00    1.727324E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.549120186 0.454560093 0.395439907 5      7      16      36
       1.00   1.025806874     1.334574E+00    1.707623E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.54        0.45        0.4         5      7      16      36
       1.25   1.014356561     1.168218E+00    1.437403E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.53796886  0.439844299 0.4         5      7      16      36
       1.50   1.006424520     1.032296E+00    1.248681E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.536309298 0.431546488 0.4         5      7      16      36
       2.00   0.990611915     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.533690702 0.418453512 0.4         5      7      16      36
       2.50   0.981034792     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.531659562 0.408297812 0.4         5      7      16      36
       3.00   0.978562884     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.53        0.4         0.4         5      7      16      36
       4.00   0.958470267     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.53        0.4         0.4         5      7      16      36
       5.00   0.943169770     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.53        0.4         0.4         5      7      16      36
       """)


class ZhaoEtAl2006AscSWISS08(ZhaoEtAl2006Asc):

    """
    --------------------------------------------------------------------
    This class implments an extension of the ZhaoEtAl (2006Asc) model,
    adjusted to be used for the new Swiss Hazard Model [2014].
    1) kappa value
       K-adjustments corresponding to model 08(upper model) - as prepared by Ben Edwards
       K-value for PGA were not provided but infered from SA[]0.01s]
       the model considers a fixed value of vs30=1100m/s
    2) small-magnitude correction
    3) single station sigma - mean inter-event adjustment
    4) single station sigma - inter-event magnitude/distance dependent
    ------------------------------------------------------------------------
    Disclaimer: these equations are modified to be used for the
    new Swiss Seismic Hazard Model [2014].
    The use of these models is the soly responsability of the hazard modeler.
    --------------------------------------------------------------------
    Model implmented by laurentiu.danciu@sed.ethz.ch
    --------------------------------------------------------------------
    """

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):

        C_ADJ = self.COEFFS_FS_ROCK[imt]
        C1_rrup = self._compute_C1_term(C_ADJ, imt, dists)
        phi_ss = self._compute_phi_ss(C_ADJ, rup, C1_rrup, imt)

        mean, stddevs = super(ZhaoEtAl2006AscSWISS08, self).\
            get_mean_and_stddevs(sites, rup, dists, imt, stddev_types)

        #: apply k-correction corresponding to the lower model [01]
        mean_corr = np.exp(
            mean) * C_ADJ['k_adj'] * self._compute_small_mag_correction_term(C_ADJ, rup.mag, dists.rrup)
        mean = np.log(mean_corr)

        std_corr = self._get_corr_stddevs(
            self.COEFFS_ASC[imt], stddev_types, len(sites.vs30), phi_ss)
        stddevs = np.array(std_corr)

        stddevs = np.array(stddevs)
        return mean, stddevs

    def _get_corr_stddevs(self, C, stddev_types, num_sites, phi_ss):
        """
        Return standard deviations adjusted for single station sigma
        as proposed to be used in the new Swiss Hazard Model [2014].
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(
                    np.sqrt(
                        C['tauC'] ** 2 +
                        phi_ss ** 2) +
                    np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(C['sigma'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(C['tauC'] + np.zeros(num_sites))
        return stddevs

    def _compute_small_mag_correction_term(self, C, mag, rrup):
        if mag >= 3.00 and mag < 5.5:
            return 1 / np.exp(((5.50 - mag) / C['a1']) ** C['a2'] * (
                C['b1'] + C['b2'] * np.log(np.maximum(np.minimum(rrup, C['Rm']), 10) / 20)))
        elif mag >= 5.50:
            return 1
        else:
            return 1

    def _compute_C1_term(self, C, imt, dists):
        """
        Return C1 coeffs as function of Rrup as proposed by Rodriguez-Marek et al (2013)
        The C1 coeff are used to compute the single station sigma
        """
        C1_rrup = 0.0
        if (dists.rrup < C['Rc11']).any():
            C1_rrup = C['phi_11']
        elif ((dists.rrup >= C['Rc11']).any()
              and (dists.rrup <= C['Rc21']).any()):
            C1_rrup = C['phi_11'] + (C['phi_21'] - C['phi_11']) * \
                ((dists.rrup - C['Rc11']) / (C['Rc21'] - C['Rc11']))
        elif (dists.rrup > C['Rc21']).any():
            C1_rrup = C['phi_21']
        return C1_rrup

    def _compute_phi_ss(self, C, rup, C1_rrup, imt):
        """
        Return C1 coeffs as function of Rrup as proposed by Rodriguez-Marek et al (2013)
        The C1 coeff are used to compute the single station sigma
        """
        phi_ss = 0

        if rup.mag < C['Mc1']:
            phi_ss = C1_rrup
        elif rup.mag >= C['Mc1'] and rup.mag <= C['Mc2']:
            phi_ss = C1_rrup + \
                (C['C2'] - C1_rrup) * \
                ((rup.mag - C['Mc1']) / (C['Mc2'] - C['Mc1']))
        elif rup.mag > C['Mc2']:
            phi_ss = C['C2']
        return phi_ss

    COEFFS_FS_ROCK = CoeffsTable(sa_damping=5, table="""\
       IMT    k_adj           a1              a2              b1              b2              Rm              phi_11      phi_21      C2          Mc1    Mc2    Rc11    Rc21
       pga    1.414560000     1.642085E+00    1.833001E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.58        0.47        0.35        5      7      16      36
       0.05   2.012007281     1.590620E+00    1.617332E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.5520412   0.4490309   0.4059176   5      7      16      36
       0.10   1.363140802     1.478437E+00    1.608714E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.54        0.44        0.43        5      7      16      36
       0.15   1.143182969     1.488927E+00    1.732820E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.580947375 0.47509775  0.400751875 5      7      16      36
       0.20   1.039739290     1.496370E+00    1.826639E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.61        0.5         0.38        5      7      16      36
       0.25   0.983550465     1.473785E+00    1.854356E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.626510191 0.5         0.374496603 5      7      16      36
       0.30   0.948764154     1.455332E+00    1.877313E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.64        0.5         0.37        5      7      16      36
       0.40   0.913557081     1.428722E+00    1.886767E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.617473168 0.488736584 0.37        5      7      16      36
       0.50   0.891456331     1.405794E+00    1.841480E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.6         0.48        0.37        5      7      16      36
       0.60   0.881236402     1.387061E+00    1.805286E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.584217936 0.472108968 0.377891032 5      7      16      36
       0.70   0.877081048     1.371222E+00    1.775240E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.57087439  0.465437195 0.384562805 5      7      16      36
       0.80   0.876720492     1.357502E+00    1.749618E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.559315686 0.459657843 0.390342157 5      7      16      36
       0.90   0.875316976     1.345400E+00    1.727324E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.549120186 0.454560093 0.395439907 5      7      16      36
       1.00   0.880663480     1.334574E+00    1.707623E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.54        0.45        0.4         5      7      16      36
       1.25   0.887163571     1.168218E+00    1.437403E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.53796886  0.439844299 0.4         5      7      16      36
       1.50   0.897498847     1.032296E+00    1.248681E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.536309298 0.431546488 0.4         5      7      16      36
       2.00   0.910434478     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.533690702 0.418453512 0.4         5      7      16      36
       2.50   0.918285853     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.531659562 0.408297812 0.4         5      7      16      36
       3.00   0.926223751     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.53        0.4         0.4         5      7      16      36
       4.00   0.922398917     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.53        0.4         0.4         5      7      16      36
       5.00   0.919443026     8.178257E-01    1.000000E+00    1.000000E+00    0.000000E+00    1.000000E+09    0.53        0.4         0.4         5      7      16      36
       """)


class ZhaoEtAl2006AscSWISS05T(ZhaoEtAl2006AscSWISS05):

    """
    --------------------------------------------------------------------
    This class implments an extension of the ZhaoEtAl (2006Asc) model,
    adjusted to be used for the new Swiss Hazard Model [2014].
    1) kappa value
       K-adjustments corresponding to model 05 (lower model) - as prepared by Ben Edwards
       K-value for PGA were not provided but infered from SA[]0.01s]
       the model considers a fixed value of vs30=1100m/s
    2) small-magnitude correction
    3) single station sigma - mean inter-event adjustment
    4) single station sigma - inter-event magnitude/distance dependent
    --------------------------------------------------------------------
    This implmentation of the AkB2010 Model considers the mean inter-event
    adjustement when computing the single station sigma (reported as total
    standard deviation))
    ------------------------------------------------------------------------
    Disclaimer: these equations are modified to be used for the
    new Swiss Seismic Hazard Model [2014].
    The use of these models is the soly responsability of the hazard modeler.
    --------------------------------------------------------------------
    Model implmented by laurentiu.danciu@sed.ethz.ch
    --------------------------------------------------------------------
    """

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        Adjust the meadian value to the soil-type used for
        Swiss hazard Vs30=1100m/s
        """
        C_ADJ = self.COEFFS_PHI_SS[imt]

        mean, stddevs = super(
            ZhaoEtAl2006AscSWISS05T, self).get_mean_and_stddevs(
            sites, rup, dists, imt, stddev_types)

        std_corr = self._get_corr_stddevs(
            self.COEFFS_ASC[imt], stddev_types, len(
                sites.vs30), C_ADJ['phi_ss'])
        stddevs = np.array(std_corr)

        return mean, stddevs

    def _get_corr_stddevs(self, C, stddev_types, num_sites, phi_ss):
        """
        Return standard deviations adjusted for single station sigma
        as proposed to be used in the new Swiss Hazard Model [2014].
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(
                    np.sqrt(
                        C['tauC'] ** 2 +
                        phi_ss ** 2) +
                    np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(C['sigma'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(C['tauC'] + np.zeros(num_sites))
            return stddevs

    COEFFS_PHI_SS = CoeffsTable(sa_damping=5, table="""\
    IMT     phi_ss
    pga     0.460
    0.050    0.453
    0.100    0.450
    0.150    0.468
    0.200    0.480
    0.250    0.480
    0.300    0.480
    0.400    0.469
    0.500    0.460
    0.600    0.457
    0.700    0.455
    0.800    0.453
    0.900    0.452
    1.000    0.450
    1.250    0.442
    1.500    0.435
    2.000    0.425
    2.500    0.417
    3.000    0.410
    4.000    0.410
    5.000    0.410
    """)


class ZhaoEtAl2006AscSWISS03T(ZhaoEtAl2006AscSWISS03):

    """
    --------------------------------------------------------------------
    This class implments an extension of the ZhaoEtAl (2006Asc) model,
    adjusted to be used for the new Swiss Hazard Model [2014].
    1) kappa value
       K-adjustments corresponding to model 03 (mid model) - as prepared by Ben Edwards
       K-value for PGA were not provided but infered from SA[]0.01s]
       the model considers a fixed value of vs30=1100m/s
    2) small-magnitude correction
    3) single station sigma - mean inter-event adjustment
    4) single station sigma - inter-event magnitude/distance dependent
    --------------------------------------------------------------------
    This implmentation of the AkB2010 Model considers the mean inter-event
    adjustement when computing the single station sigma (reported as total
    standard deviation))
    ------------------------------------------------------------------------
    Disclaimer: these equations are modified to be used for the
    new Swiss Seismic Hazard Model [2014].
    The use of these models is the soly responsability of the hazard modeler.
    --------------------------------------------------------------------
    Model implmented by laurentiu.danciu@sed.ethz.ch
    --------------------------------------------------------------------
    """

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        Adjust the meadian value to the soil-type used for
        Swiss hazard Vs30=1100m/s
        """
        C_ADJ = self.COEFFS_PHI_SS[imt]

        mean, stddevs = super(
            ZhaoEtAl2006AscSWISS03T, self).get_mean_and_stddevs(
            sites, rup, dists, imt, stddev_types)

        std_corr = self._get_corr_stddevs(
            self.COEFFS_ASC[imt], stddev_types, len(
                sites.vs30), C_ADJ['phi_ss'])
        stddevs = np.array(std_corr)
        return mean, stddevs

    def _get_corr_stddevs(self, C, stddev_types, num_sites, phi_ss):
        """
        Return standard deviations adjusted for single station sigma
        as proposed to be used in the new Swiss Hazard Model [2014].
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(
                    np.sqrt(
                        C['tauC'] ** 2 +
                        phi_ss ** 2) +
                    np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(C['sigma'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(C['tauC'] + np.zeros(num_sites))
            return stddevs

    COEFFS_PHI_SS = CoeffsTable(sa_damping=5, table="""\
    IMT     phi_ss
    pga     0.460
    0.050    0.453
    0.100    0.450
    0.150    0.468
    0.200    0.480
    0.250    0.480
    0.300    0.480
    0.400    0.469
    0.500    0.460
    0.600    0.457
    0.700    0.455
    0.800    0.453
    0.900    0.452
    1.000    0.450
    1.250    0.442
    1.500    0.435
    2.000    0.425
    2.500    0.417
    3.000    0.410
    4.000    0.410
    5.000    0.410
    """)


class ZhaoEtAl2006AscSWISS08T(ZhaoEtAl2006AscSWISS08):

    """
    --------------------------------------------------------------------
    This class implments an extension of the ZhaoEtAl (2006Asc) model,
    adjusted to be used for the new Swiss Hazard Model [2014].
    1) kappa value
       K-adjustments corresponding to model 08 (upper model) - as prepared by Ben Edwards
       K-value for PGA were not provided but infered from SA[]0.01s]
       the model considers a fixed value of vs30=1100m/s
    2) small-magnitude correction
    3) single station sigma - mean inter-event adjustment
    4) single station sigma - inter-event magnitude/distance dependent
    --------------------------------------------------------------------
    This implmentation of the AkB2010 Model considers the mean inter-event
    adjustement when computing the single station sigma (reported as total
    standard deviation))
    ------------------------------------------------------------------------
    Disclaimer: these equations are modified to be used for the
    new Swiss Seismic Hazard Model [2014].
    The use of these models is the soly responsability of the hazard modeler.
    --------------------------------------------------------------------
    Model implmented by laurentiu.danciu@sed.ethz.ch
    --------------------------------------------------------------------
    """

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        Adjust the meadian value to the soil-type used for
        Swiss hazard Vs30=1100m/s
        """
        C_ADJ = self.COEFFS_PHI_SS[imt]

        mean, stddevs = super(
            ZhaoEtAl2006AscSWISS08T, self).get_mean_and_stddevs(
            sites, rup, dists, imt, stddev_types)

        std_corr = self._get_corr_stddevs(
            self.COEFFS_ASC[imt], stddev_types, len(
                sites.vs30), C_ADJ['phi_ss'])
        stddevs = np.array(std_corr)
        return mean, stddevs

    def _get_corr_stddevs(self, C, stddev_types, num_sites, phi_ss):
        """
        Return standard deviations adjusted for single station sigma
        as proposed to be used in the new Swiss Hazard Model [2014].
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(
                    np.sqrt(
                        C['tauC'] ** 2 +
                        phi_ss ** 2) +
                    np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(C['sigma'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(C['tauC'] + np.zeros(num_sites))
            return stddevs

    COEFFS_PHI_SS = CoeffsTable(sa_damping=5, table="""\
    IMT     phi_ss
    pga     0.460
    0.050    0.453
    0.100    0.450
    0.150    0.468
    0.200    0.480
    0.250    0.480
    0.300    0.480
    0.400    0.469
    0.500    0.460
    0.600    0.457
    0.700    0.455
    0.800    0.453
    0.900    0.452
    1.000    0.450
    1.250    0.442
    1.500    0.435
    2.000    0.425
    2.500    0.417
    3.000    0.410
    4.000    0.410
    5.000    0.410
    """)
