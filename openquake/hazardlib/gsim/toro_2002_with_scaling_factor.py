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
Module exports
:class:`ToroEtAl2002_SF_0_40`,
:class:`ToroEtAl2002_SF_0_50`,
:class:`ToroEtAl2002_SF_0_75`,
:class:`ToroEtAl2002_SF_1_33`,
:class:`ToroEtAl2002_SF_2_50`,
:class:`ToroEtAl2002ANGRA_SF_0_40`,
:class:`ToroEtAl2002ANGRA_SF_0_50`,
:class:`ToroEtAl2002ANGRA_SF_0_75`,
:class:`ToroEtAl2002ANGRA_SF_1_33`,
:class:`ToroEtAl2002ANGRA_SF_2_50`.
"""
from __future__ import division

import numpy as np
from openquake.hazardlib.gsim.toro_2002 import ToroEtAl2002, ToroEtAl2002ANGRA

class ToroEtAl2002_SF_0_40(ToroEtAl2002):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 0.4
        mean, stddevs = super(ToroEtAl2002_SF_0_40, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs

class ToroEtAl2002_SF_0_50(ToroEtAl2002):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 0.5
        mean, stddevs = super(ToroEtAl2002_SF_0_50, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs
		
class ToroEtAl2002_SF_0_75(ToroEtAl2002):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 0.75
        mean, stddevs = super(ToroEtAl2002_SF_0_75, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs
        
class ToroEtAl2002_SF_1_33(ToroEtAl2002):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 1.33
        mean, stddevs = super(ToroEtAl2002_SF_1_33, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs
        
class ToroEtAl2002_SF_2_50(ToroEtAl2002):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 2.50
        mean, stddevs = super(ToroEtAl2002_SF_2_50, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs

class ToroEtAl2002ANGRA_SF_0_40(ToroEtAl2002ANGRA):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 0.4
        mean, stddevs = super(ToroEtAl2002ANGRA_SF_0_40, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs

class ToroEtAl2002ANGRA_SF_0_50(ToroEtAl2002ANGRA):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 0.5
        mean, stddevs = super(ToroEtAl2002ANGRA_SF_0_50, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs
		
class ToroEtAl2002ANGRA_SF_0_75(ToroEtAl2002ANGRA):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 0.75
        mean, stddevs = super(ToroEtAl2002ANGRA_SF_0_75, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs
        
class ToroEtAl2002ANGRA_SF_1_33(ToroEtAl2002ANGRA):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 1.33
        mean, stddevs = super(ToroEtAl2002ANGRA_SF_1_33, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs
        
class ToroEtAl2002ANGRA_SF_2_50(ToroEtAl2002ANGRA):
    """
    Applies the GMPE with a modified scaling factor that
    is frequency dependent
    """
    def get_mean_and_stddevs(self, sctx, rctx, dctx, imt, full_stddevs):
        """
        """
        scale_factor = 2.50
        mean, stddevs = super(ToroEtAl2002ANGRA_SF_2_50, self).get_mean_and_stddevs(
            sctx, rctx, dctx, imt, full_stddevs)
        return mean + np.log(scale_factor), stddevs
