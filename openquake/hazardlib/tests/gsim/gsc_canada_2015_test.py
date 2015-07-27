# The Hazard Library
# Copyright (C) 2015, GEM Foundation
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

from openquake.hazardlib.gsim.gsc_canada_2015 import (
    GSCCanada2015WCrustRhypoMed,
    GSCCanada2015WCrustRhypoLow,
    GSCCanada2015WCrustRhypoHigh,
    GSCCanada2015WCrustRjbMed,
    GSCCanada2015WCrustRjbLow,
    GSCCanada2015WCrustRjbHigh,
    GSCCanada2015ENAMed,
    GSCCanada2015ENALow,
    GSCCanada2015ENAHigh,
    GSCCanada2015WOffshoreMed,
    GSCCanada2015WOffshoreLow,
    GSCCanada2015WOffshoreHigh,
    GSCCanada2015SInterMed,
    GSCCanada2015SInterLow,
    GSCCanada2015SInterHigh,
    GSCCanada2015SSlabD30Med,
    GSCCanada2015SSlabD30Low,
    GSCCanada2015SSlabD30High,
    GSCCanada2015SSlabD50Med,
    GSCCanada2015SSlabD50Low,
    GSCCanada2015SSlabD50High
    )

from openquake.hazardlib.tests.gsim.utils import BaseGSIMTestCase


# Western Canada - Distributed Sources
class GSCCanada2015WCrustRhypoMedTestCase(BaseGSIMTestCase):
    """
    Tables generated from the original tables file `Wcrust_med_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WCrustRhypoMed
    MEAN_FILE = "GSC2015/Wcrust_rhypo_med_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Wcrust_rhypo_med_TOTAL.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE,
                   max_discrep_percentage=0.7)

    def test_std_total(self):
        self.check(self.STD_TOTAL_FILE,
                   max_discrep_percentage=0.7)


class GSCCanada2015WCrustRhypoLowTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `Wcrust_low_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WCrustRhypoLow
    MEAN_FILE = "GSC2015/Wcrust_rhypo_low_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Wcrust_rhypo_low_TOTAL.csv"


class GSCCanada2015WCrustRhypoHighTestCase(
        GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `Wcrust_high_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WCrustRhypoHigh
    MEAN_FILE = "GSC2015/Wcrust_rhypo_high_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Wcrust_rhypo_high_TOTAL.csv"


# Western Canada - Fault Sources
class GSCCanada2015WCrustRjbMedTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WcrustFRjb_med_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WCrustRjbMed
    MEAN_FILE = "GSC2015/Wcrust_rjb_med_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Wcrust_rjb_med_TOTAL.csv"


class GSCCanada2015WCrustRjbLowTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WcrustFRjb_low_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WCrustRjbLow
    MEAN_FILE = "GSC2015/Wcrust_rjb_low_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Wcrust_rjb_low_TOTAL.csv"


class GSCCanada2015WCrustRjbHighTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WcrustFRjb_high_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WCrustRjbHigh
    MEAN_FILE = "GSC2015/Wcrust_rjb_high_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Wcrust_rjb_high_TOTAL.csv"


# Eastern North America
class GSCCanada2015ENAMedTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `ENA_med_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015ENAMed
    MEAN_FILE = "GSC2015/ENA_rhypo_med_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/ENA_rhypo_med_TOTAL.csv"


class GSCCanada2015ENALowTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `ENA_low_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015ENALow
    MEAN_FILE = "GSC2015/ENA_rhypo_low_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/ENA_rhypo_low_TOTAL.csv"


class GSCCanada2015ENAHighTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `ENA_high_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015ENAHigh
    MEAN_FILE = "GSC2015/ENA_rhypo_high_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/ENA_rhypo_high_TOTAL.csv"


# Western Canada - Offshore Sources
class GSCCanada2015WOffshoreMedTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `Woffshore_med_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WOffshoreMed
    MEAN_FILE = "GSC2015/Woffshore_rhypo_med_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Woffshore_rhypo_med_TOTAL.csv"


class GSCCanada2015WOffshoreLowTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `Woffshore_low_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WOffshoreLow
    MEAN_FILE = "GSC2015/Woffshore_rhypo_low_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Woffshore_rhypo_low_TOTAL.csv"


class GSCCanada2015WOffshoreHighTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `Woffshore_high_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015WOffshoreHigh
    MEAN_FILE = "GSC2015/Woffshore_rhypo_high_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Woffshore_rhypo_high_TOTAL.csv"


# Interface
class GSCCanada2015SInterMedTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WinterfaceCombo_medclC.txt`
    """
    GSIM_CLASS = GSCCanada2015SInterMed
    MEAN_FILE = "GSC2015/Winterface_rrup_med_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Winterface_rrup_med_TOTAL.csv"


class GSCCanada2015SInterLowTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WinterfaceCombo_lowclC.txt`
    """
    GSIM_CLASS = GSCCanada2015SInterLow
    MEAN_FILE = "GSC2015/Winterface_rrup_low_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Winterface_rrup_low_TOTAL.csv"


class GSCCanada2015SInterHighTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file
    `WinterfaceCombo_highclC.txt`
    """
    GSIM_CLASS = GSCCanada2015SInterHigh
    MEAN_FILE = "GSC2015/Winterface_rrup_high_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/Winterface_rrup_high_TOTAL.csv"

# In-Slab Sources - 30 km Depth
class GSCCanada2015SSlabD30MedTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WinslabD30_med_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015SSlabD30Med
    MEAN_FILE = "GSC2015/WinslabD30_rhypo_med_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/WinslabD30_rhypo_med_TOTAL.csv"


class GSCCanada2015SSlabD30LowTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WinslabD30_low_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015SSlabD30Low
    MEAN_FILE = "GSC2015/WinslabD30_rhypo_low_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/WinslabD30_rhypo_low_TOTAL.csv"


class GSCCanada2015SSlabD30HighTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WinslabD30_high_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015SSlabD30High
    MEAN_FILE = "GSC2015/WinslabD30_rhypo_high_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/WinslabD30_rhypo_high_TOTAL.csv"


# In-Slab Sources - 50 km Depth
class GSCCanada2015SSlabD50MedTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WinslabD50_med_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015SSlabD50Med
    MEAN_FILE = "GSC2015/WinslabD50_rhypo_med_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/WinslabD50_rhypo_med_TOTAL.csv"


class GSCCanada2015SSlabD50LowTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WinslabD50_low_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015SSlabD50Low
    MEAN_FILE = "GSC2015/WinslabD50_rhypo_low_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/WinslabD50_rhypo_low_TOTAL.csv"


class GSCCanada2015SSlabD50HighTestCase(GSCCanada2015WCrustRhypoMedTestCase):
    """
    Tables generated from the original tables file `WinslabD50_high_clC.txt`
    """
    GSIM_CLASS = GSCCanada2015SSlabD50High
    MEAN_FILE = "GSC2015/WinslabD50_rhypo_high_MEAN.csv"
    STD_TOTAL_FILE = "GSC2015/WinslabD50_rhypo_high_TOTAL.csv"
