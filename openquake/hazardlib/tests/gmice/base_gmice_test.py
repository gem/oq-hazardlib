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

import unittest
import numpy as np
from openquake.hazardlib.gsim.base import (SitesContext,
                                           DistancesContext,
                                           RuptureContext)
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, MMI
# Needs a test GMPE from an active shallow region with a magnitude and rrup
# dependence and as little else as possible!
from openquake.hazardlib.gsim.cauzzi_2014 import CauzziEtAl2014NoSOF
from openquake.hazardlib.gsim.allen_2012_ipe import AllenEtAl2012
from openquake.hazardlib.gmice.base import GMICE
from openquake.hazardlib.gmice.worden_2012 import WordenEtAl2012


class FakeGMICE(GMICE):
    DEFINED_FOR_TECTONIC_REGION_TYPE = None
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set()
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = None
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set()
    REQUIRES_SITES_PARAMETERS = set()
    REQUIRES_RUPTURE_PARAMETERS = set()
    REQUIRES_DISTANCES = set()
    INTENSITY_RANGE = (0, 12)
    SOURCE_NAME = "A Fake GMICE"
    SCALE = "Something"
    

    def get_mean_intensity_and_stddevs(self, imls, sites, rup, dists,
                                       imt, stddev_types):
        pass
     
    def get_mean_gm_and_stddevs(self, imls, sites, rup, dists, imt,
                                stddev_types):
        pass

class BaseGMICETestCase(unittest.TestCase):
    """
    General set of tests for the :class: openquake.hazardlib.gmice.base.GMICE
    """
    def setUp(self):

        self.gmice_class = WordenEtAl2012
        self.gmice = self.gmice_class()
        self.rctx = RuptureContext()
        self.rctx.mag = 6.0
        self.dctx = DistancesContext()
        self.dctx.rrup = np.array([5.0, 10.0, 20.0, 50.0, 100.0])
        self.sctx = SitesContext()
        self.sctx.vs30 = 760.0 * np.ones_like(self.dctx.rrup)
        self.imt = PGA()
        self.stddevs = [const.StdDev.TOTAL]

    def _assert_value_error(self, func, error, **kwargs):
        with self.assertRaises(ValueError) as ar:
            func(**kwargs)
        self.assertEqual(str(ar.exception), error)

    def test_correct_property_constructions(self):
        dummy_gmice = FakeGMICE()
        self.assertEqual(dummy_gmice.SCALE, "Something")
        self.assertEqual(dummy_gmice.SOURCE_NAME, "A Fake GMICE")
        self.assertEqual(dummy_gmice.INTENSITY_RANGE, (0, 12))


    def test_compare_gm_to_mmi_gsim(self):
        """
        Compare Ground Motion to Intensity Conversions using the GMPE
        approach and the normal approach
        """
        # Get some ground motion levels from GMPE
        gmpe = CauzziEtAl2014NoSOF()
        imls, _ = gmpe.get_mean_and_stddevs(self.sctx, self.rctx, self.dctx,
                                            self.imt, self.stddevs)
        # Get MMIs from intensity measure levels
        mmis_imls, _ = self.gmice.get_mean_intensity_and_stddevs(np.exp(imls),
                                                                 self.sctx,
                                                                 self.rctx,
                                                                 self.dctx,
                                                                 self.imt,
                                                                 self.stddevs)
        # Get MMIs from the GMPE
        mmis_gmpe, _ = self.gmice.get_mean_intensity_and_stddevs_from_gmpe(
            gmpe,
            self.sctx,
            self.rctx,
            self.dctx,
            self.imt,
            self.stddevs)
        # Check that the results are equal
        np.testing.assert_array_almost_equal(mmis_imls, mmis_gmpe, 7)

    def test_compare_mmi_to_gm_ipe(self):
        """
        Compare MMI to ground motion conversion using an IPE approach
        """
        # Get some intensity measure levels from an IPE
        ipe = AllenEtAl2012()
        mmis, _ = ipe.get_mean_and_stddevs(self.sctx, self.rctx, self.dctx,
                                           MMI(), self.stddevs)
        # Get ground motion levels (PGA) from GMICE
        gmls, _ = self.gmice.get_mean_gm_and_stddevs(mmis, self.sctx,
                                                     self.rctx, self.dctx,
                                                     self.imt, self.stddevs)
        gmls_ipe, _ = self.gmice.get_mean_gm_and_stddevs_from_ipe(ipe,
                                                                  self.sctx,
                                                                  self.rctx,
                                                                  self.dctx,
                                                                  self.imt,
                                                                  self.stddevs)
        np.testing.assert_array_almost_equal(gmls, gmls_ipe, 7)

    def test_call_to_get_mean_and_stddevs(self):
        """
        The method get_mean_and_stddevs is turned off for GMICEs - check
        that the Not Implemented Error is raised
        """
        with self.assertRaises(NotImplementedError) as nie:
            _ = self.gmice.get_mean_and_stddevs(self.sctx,
                                                self.rctx,
                                                self.dctx,
                                                self.imt,
                                                self.stddevs)
        self.assertEqual(
            str(nie.exception),
            "Method 'get_mean_and_stddevs' not available in GMICE object")

