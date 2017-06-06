# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2017 GEM Foundation
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
Module :mod:`openquake.hazardlib.scalerel.stafford2014` implements :class:`Stafford2014`.
"""
from math import log, exp, sin, radians, pi, sqrt
from scipy.stats import norm, chi2
from openquake.hazardlib.scalerel.base import BaseMSRSigma


class Stafford2014(BaseMSRSigma):
    """
    Stafford magnitude -- rupture dimension relationships,
    see 2014, Bull. Seism. Soc. Am., vol. 104, no. 4, pages 1620-1635.

    Implements magnitude-area scaling relationships.
    """
    def get_median_area(self, mag, rake, dip, seismogenic_depth):
        """
        The values are a bilinear function of magnitude, 
        with the break in scaling depending upon the rake, dip and seismogenic depth.

        """
        assert -180 <= rake <= 180
        max_log_width = log(seismogenic_depth / sin(radians(dip)))

        if (-30 <= rake <= 30) or (rake >= 150) or (rake <= -150):
            # strike slip
            # coefficients for the mean rupture width
            beta0 = -2.3000
            beta1 = 0.7167
            # coefficients for the area calculations
            gamma0 = -9.3137
        elif rake > 0:
            # thrust/reverse
            # coefficients for the mean rupture width
            beta0 = -3.8300
            beta1 = 0.9982
            # coefficients for the area calculations
            gamma0 = -9.2749
        else:
            # normal
            # coefficients for the mean rupture width
            beta0 = -4.1055
            beta1 = 1.0370
            # coefficients for the area calculations
            gamma0 = -9.2483

        mag_crit = (max_log_width - beta0) / beta1
        mean_log_area = gamma0 + log(10.0) * mag
        if mag > mag_crit:
            mean_log_area -= log(10.0)*(mag - mag_crit)/4.0

        return exp(mean_log_area)

    def get_std_dev_area(self, mag, rake, dip, seismogenic_depth):
        """
        Standard deviation for Stafford2014. 
        Relatively complicated expression due to effects of censoring upon the rupture width
        """
        assert -180 <= rake <= 180
        max_log_width = log(seismogenic_depth / sin(radians(dip)))

        if (-30 <= rake <= 30) or (rake >= 150) or (rake <= -150):
            # strike slip
            # coefficients for prob of exceeding max width
            alpha0 = -30.8395
            alpha1 = 5.4184
            alpha2 = -0.3044
            # coefficients for the mean rupture width
            beta0 = -2.3000
            beta1 = 0.7167
            # uncensored variability
            sigma_width_ucens = 0.2337
            # coefficients for the area calculations
            sigma_length = 0.3138
            rho_length_width = 0.3104
        elif rake > 0:
            # thrust/reverse
            # coefficients for prob of exceeding max width
            alpha0 = -35.8239
            alpha1 = 5.0680
            alpha2 = -0.0457
            # coefficients for the mean rupture width
            beta0 = -3.8300
            beta1 = 0.9982
            # uncensored variability
            sigma_width_ucens = 0.2285
            # coefficients for the area calculations
            sigma_length = 0.2534
            rho_length_width = 0.1376
        else:
            # normal
            # coefficients for prob of exceeding max width
            alpha0 = -36.9770
            alpha1 = 6.3070
            alpha2 = -0.1696
            # coefficients for the mean rupture width
            beta0 = -4.1055
            beta1 = 1.0370
            # uncensored variability
            sigma_width_ucens = 0.2509
            # coefficients for the area calculations
            sigma_length = 0.3454
            rho_length_width = 0.4336

        # transformed parameters for width prediction
        z = alpha0 + alpha1 * mag + alpha2 * exp(max_log_width)
        p = 1 / (1 + exp(-z))

        # censored and uncensored mean
        mean_log_width_ucens = beta0 + beta1 * mag
        # define a random variate
        psi = (max_log_width - mean_log_width_ucens) / sigma_width_ucens

        # define the chi-squared cdf value for psi**2
        psi2_cdf = chi2.cdf(psi**2, 3)

        # define the uncensored variance of the width
        var_width_ucens = sigma_width_ucens**2
        # get normal cdf for psi
        psi_cdf = norm.cdf(psi)
        # define exponential function of psi
        exp_psi = exp(-psi**2) / psi_cdf

        # compute the censored variance for the rupture width
        if psi < 0:
            var_width_cens = var_width_ucens / (2*pi * psi_cdf) * (pi*(1 - psi2_cdf) - exp_psi)
        else:
            var_width_cens = var_width_ucens / (2*pi * psi_cdf) * (pi*(1 + psi2_cdf) - exp_psi)

        # define the censored standard deviation
        sigma_width_cens = sqrt(var_width_cens)

        # compute the standard deviation of the width
        sigma_width = (1 - p) * sigma_width_cens

        # variance of the area
        var_area = sigma_length**2 + sigma_width**2 + 2 * rho_length_width * sigma_length * sigma_width

        # return the standard deviation
        return sqrt(var_area)
