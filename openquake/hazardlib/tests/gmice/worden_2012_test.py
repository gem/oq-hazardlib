# The Hazard Library
# Copyright (C) 2013-2016 GEM Foundation
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
Implements the tests for the Worden et al (2012) GMICE

Test data generated from original Shakemap implementation
"""


from openquake.hazardlib.gmice.worden_2012 import (WordenEtAl2012,
                                                   WordenEtAl2012MagDist)
from openquake.hazardlib.tests.gmice.utils import BaseGMICETestCase

class WordenEtAl2012TestCase(BaseGMICETestCase):
    GMICE_CLASS = WordenEtAl2012

    def test_mean_gmls_to_mmi(self):
        self.check("wgrw12/WORDEN2012_GMLS_MMI_MEAN.csv",
                   max_discrep_percentage=0.1)

    def test_stddev_gmls_to_mmi(self):
        self.check("wgrw12/WORDEN2012_GMLS_MMI_TOTAL_STDDEV.csv",
                   max_discrep_percentage=0.1)

    def test_mean_mmi_to_gmls(self):
        self.check("wgrw12/WORDEN2012_MMI_GMLS_MEAN.csv",
                   max_discrep_percentage=0.1)
    
    def test_stddev_mmi_to_gmls(self):
        self.check("wgrw12/WORDEN2012_MMI_GMLS_TOTAL_STDDEV.csv",
                   max_discrep_percentage=0.1)


class WordenEtAl2012MagDistTestCase(BaseGMICETestCase):
    GMICE_CLASS = WordenEtAl2012MagDist

    def test_mean_gmls_to_mmi(self):
        self.check("wgrw12/WORDEN2012MR_GMLS_MMI_MEAN.csv",
                   max_discrep_percentage=0.1)

    def test_stddev_gmls_to_mmi(self):
        self.check("wgrw12/WORDEN2012MR_GMLS_MMI_TOTAL_STDDEV.csv",
                   max_discrep_percentage=0.1)

    def test_mean_mmi_to_gmls(self):
        self.check("wgrw12/WORDEN2012MR_MMI_GMLS_MEAN.csv",
                   max_discrep_percentage=0.1)
    
    def test_stddev_mmi_to_gmls(self):
        self.check("wgrw12/WORDEN2012MR_MMI_GMLS_TOTAL_STDDEV.csv",
                   max_discrep_percentage=0.1)

