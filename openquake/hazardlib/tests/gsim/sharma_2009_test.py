# The Hazard Library
# Copyright (C) 2012-2015, GEM Foundation
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

from openquake.hazardlib.gsim.sharma_2009 import SharmaEtAl2009
from openquake.hazardlib.tests.gsim.utils import BaseGSIMTestCase


class SharmaEtAl2009TestCase(BaseGSIMTestCase):
    GSIM_CLASS = SharmaEtAl2009

    # Test data were obtained after visual verification of results
    # against versions of Figures 6,7 and 9 digitized using
    # http://arohatgi.info/WebPlotDigitizer/ app/. Note that the
    # correspondance was excellent at 0.04 s for all magitudes,
    # distances, mechanisms, and site conditions, but for other periods
    # the results diverge, up to about 30%. This test data file was
    # generated from the current implementation, not that of the
    # original authors.
    def test_mean(self):
        self.check('SDBK09/SDBK09_MEAN.csv',
                   max_discrep_percentage=0.01)

    # There is no plot of residuals as a function of frequency, so
    # there's absolutely nothing to verify that against.
    def test_std_total(self):
        self.check('SDBK09/SDBK09_STD_TOTAL.csv',
                   max_discrep_percentage=0.01)
