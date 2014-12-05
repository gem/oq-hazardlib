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
import unittest

from openquake.hazardlib.scalerel.bonilla1984 import (
    BonillaEtAl1984Aspect1p0,
    BonillaEtAl1984Aspect0p5,
    BonillaEtAl1984Aspect1p5,
    BonillaEtAl1984Aspect2p0,
    BonillaEtAl1984Aspect2p5,
    BonillaEtAl1984Aspect3p0,
    BonillaEtAl1984Aspect3p5,
    BonillaEtAl1984Aspect4p0,
    BonillaEtAl1984Aspect5p0,
    BonillaEtAl1984Aspect6p0,
    BonillaEtAl1984Aspect7p0,
    BonillaEtAl1984Aspect8p0)

class BonillaEtAl1984Aspect1p0TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 1.0
    """
    MSR = BonillaEtAl1984Aspect1p0
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 101.8591388,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.1020000,
              "Magnitude Sigma": 0.306}


    def setUp(self):
        """
        Instantiates the test function
        """
        self.msr = self.MSR()

    def test_get_median_area(self):
        """
        Tests the median area function
        """
        self.assertAlmostEqual(
            self.msr.get_median_area(self.PARAMS["Magnitude"], None),
            self.PARAMS["Expected Area"],
            places=7)

    def test_get_std_dev_area(self):
        """
        Tests the standard deviation for area
        """
        self.assertAlmostEqual(
            self.msr.get_std_dev_area(None),
            self.PARAMS["Area Sigma"],
            places=7)

    def test_get_median_mag(self):
        """
        Tests the median magnitude from area function
        """
        self.assertAlmostEqual(
            self.msr.get_median_mag(self.PARAMS["Area"], None),
            self.PARAMS["Expected Magnitude"],
            places=7)

    def test_get_std_dev_mag(self):
        """
        Tests the standard deviation for magnitude
        """
        self.assertAlmostEqual(
            self.msr.get_std_dev_mag(None),
            self.PARAMS["Magnitude Sigma"],
            places=7)
    

class BonillaEtAl1984Aspect0p5TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 0p5
    """
    MSR = BonillaEtAl1984Aspect0p5
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 203.7182776,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 6.9954354,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect1p5TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 1.5
    """
    MSR = BonillaEtAl1984Aspect1p5
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 67.9060925,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.1643363,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect2p0TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 2.0
    """
    MSR = BonillaEtAl1984Aspect2p0
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 50.9295694,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.2085646,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect2p5TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 2.5
    """
    MSR = BonillaEtAl1984Aspect2p5
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 40.7436555,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.2428708,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect3p0TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 3.0
    """
    MSR = BonillaEtAl1984Aspect3p0
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 33.9530463,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.2709009,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect3p5TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 3.5
    """
    MSR = BonillaEtAl1984Aspect3p5
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 29.1026111,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.2946001,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect4p0TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 4.0
    """
    MSR = BonillaEtAl1984Aspect4p0
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 25.4647847,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.3151292,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect5p0TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 5.0
    """
    MSR = BonillaEtAl1984Aspect5p0
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 20.3718278,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.3494354,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect6p0TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 6.0
    """
    MSR = BonillaEtAl1984Aspect6p0
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 16.9765231,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.3774655,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect7p0TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 7.0
    """
    MSR = BonillaEtAl1984Aspect7p0
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 14.5513055,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.4011647,
              "Magnitude Sigma": 0.306}


class BonillaEtAl1984Aspect8p0TestCase(unittest.TestCase):
    """
    This tests the implementation of Bonilla et al. (2014) magnitude
    scaling relation for the case when rupture aspect ratio is 8.0
    """
    MSR = BonillaEtAl1984Aspect8p0
    PARAMS = {"Magnitude": 6.0,
              "Expected Area": 12.7323924,
              "Area Sigma": 0.572,
              "Area": 1000.0,
              "Expected Magnitude": 7.4216939,
              "Magnitude Sigma": 0.306}
