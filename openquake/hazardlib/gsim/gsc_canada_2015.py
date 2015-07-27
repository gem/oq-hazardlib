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
"""
Full GMPE Set for the 2015 Canadian Seismic Hazard Model
Module exports :class: `GSCCanada2015WCrustRhypoMed`,
    :class:`GSCCanada2015WCrustRhypoLow`,
    :class:`GSCCanada2015WCrustRhypoHigh`,
    :class:`GSCCanada2015WCrustRjbMed`,
    :class:`GSCCanada2015WCrustRjbLow`,
    :class:`GSCCanada2015WCrustRjbHigh`,
    :class:`GSCCanada2015ENAMed`,
    :class:`GSCCanada2015ENALow`,
    :class:`GSCCanada2015ENAHigh`,
    :class:`GSCCanada2015WOffshoreMed`,
    :class:`GSCCanada2015WOffshoreLow`,
    :class:`GSCCanada2015WOffshoreHigh`,
    :class:`GSCCanada2015SInterMed`,
    :class:`GSCCanada2015SInterLow`,
    :class:`GSCCanada2015SInterHigh`,
    :class:`GSCCanada2015SSlabD30Med`,
    :class:`GSCCanada2015SSlabD30Low`,
    :class:`GSCCanada2015SSlabD30High`,
    :class:`GSCCanada2015SSlabD50Med`,
    :class:`GSCCanada2015SSlabD50Low`,
    :class:`GSCCanada2015SSlabD50High`
"""
import os
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGV, PGA, SA
from openquake.hazardlib.gsim.base import GMPETable


# Western Crustal - Rhypo
class GSCCanada2015WCrustRhypoMed(GMPETable):
    """
    Implements the Active Shallow Crustal GMPE for application to distributed
    seismicity sources in Western Canada - Middle Branch

    A complete description of the GMPEs can be found in Atkinson and Adams
    (2012) ...
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Wcrust_rhypo_med.hdf5")


class GSCCanada2015WCrustRhypoLow(GMPETable):
    """
    Implements the Active Shallow Crustal GMPE for application to distributed
    seismicity sources in Western Canada - Low Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Wcrust_rhypo_low.hdf5")


class GSCCanada2015WCrustRhypoHigh(GMPETable):
    """
    Implements the Active Shallow Crustal GMPE for application to distributed
    seismicity sources in Western Canada - High Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Wcrust_rhypo_high.hdf5")


# Western Crustal - Rjb
class GSCCanada2015WCrustRjbMed(GMPETable):
    """
    Implements the Active Shallow Crustal GMPE for application to active
    fault sources in Western Canada - Middle Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Wcrust_rjb_med.hdf5")


class GSCCanada2015WCrustRjbLow(GMPETable):
    """
    Implements the Active Shallow Crustal GMPE for application to active
    fault sources in Western Canada - Low Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Wcrust_rjb_low.hdf5")


class GSCCanada2015WCrustRjbHigh(GMPETable):
    """
    Implements the Active Shallow Crustal GMPE for application to active
    fault sources in Western Canada - High Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Wcrust_rjb_high.hdf5")


# Stable Continental GMPEs - ENA
class GSCCanada2015ENAMed(GMPETable):
    """
    Implements the Stable Shallow Crustal GMPE for application to distributed
    seismicity sources in Central and Eastern Canada - Middle Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "ENA_rhypo_med.hdf5")


class GSCCanada2015ENALow(GMPETable):
    """
    Implements the Stable Shallow Crustal GMPE for application to distributed
    seismicity sources in Central and Eastern Canada - Low Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "ENA_rhypo_low.hdf5")


class GSCCanada2015ENAHigh(GMPETable):
    """
    Implements the Stable Shallow Crustal GMPE for application to distributed
    seismicity sources in Central and Eastern Canada - High Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "ENA_rhypo_high.hdf5")


# Western Offshore GMPEs
class GSCCanada2015WOffshoreMed(GMPETable):
    """
    Implements the Active Shallow Crustal GMPE for application to distributed
    seismicity sources offshore from Western Canada - Middle Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Woffshore_rhypo_med.hdf5")


class GSCCanada2015WOffshoreLow(GMPETable):
    """
    Implements the Active Shallow Crustal GMPE for application to distributed
    seismicity sources offshore from Western Canada - Low Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Woffshore_rhypo_low.hdf5")


class GSCCanada2015WOffshoreHigh(GMPETable):
    """
    Implements the Stable Shallow Crustal GMPE for application to distributed
    seismicity sources in offshore from Western Canada - High Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Woffshore_rhypo_high.hdf5")


# Subduction Interface GMPEs
class GSCCanada2015SInterMed(GMPETable):
    """
    Implements the Subduction Interface GMPE for application to the
    Cascadia Subduction Interface Western Canada - Middle Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Interface_rrup_med.hdf5")


class GSCCanada2015SInterLow(GMPETable):
    """
    Implements the Subduction Interface GMPE for application to the
    Cascadia Subduction Interface Western Canada - Low Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Interface_rrup_low.hdf5")


class GSCCanada2015SInterHigh(GMPETable):
    """
    Implements the Subduction Interface GMPE for application to the
    Cascadia Subduction Interface Western Canada - High Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "Interface_rrup_high.hdf5")


# Subduction In-slab GMPEs (Depth = 30 km)
class GSCCanada2015SSlabD30Med(GMPETable):
    """
    Implements the Subduction Inslab GMPE for application to the Cascadia
    Subduction Slab Western Canada, assuming a hypocentral depth of 30 km
    - Middle Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTRASLAB

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "WinslabD30_rhypo_med.hdf5")


class GSCCanada2015SSlabD30Low(GMPETable):
    """
    Implements the Subduction Inslab GMPE for application to the Cascadia
    Subduction Slab Western Canada, assuming a hypocentral depth of 30 km
    - Low Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "WinslabD30_rhypo_low.hdf5")


class GSCCanada2015SSlabD30High(GMPETable):
    """
    Implements the Subduction Inslab GMPE for application to the Cascadia
    Subduction Slab Western Canada, assuming a hypocentral depth of 30 km
    - High Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "WinslabD30_rhypo_high.hdf5")


# Subduction In-slab GMPEs (Depth = 50 km)
class GSCCanada2015SSlabD50Med(GMPETable):
    """
    Implements the Subduction Inslab GMPE for application to the Cascadia
    Subduction Slab Western Canada, assuming a hypocentral depth of 50 km
    - Middle Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTRASLAB

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "WinslabD50_rhypo_med.hdf5")


class GSCCanada2015SSlabD50Low(GMPETable):
    """
    Implements the Subduction Inslab GMPE for application to the Cascadia
    Subduction Slab Western Canada, assuming a hypocentral depth of 50 km
    - Low Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "WinslabD50_rhypo_low.hdf5")


class GSCCanada2015SSlabD50High(GMPETable):
    """
    Implements the Subduction Inslab GMPE for application to the Cascadia
    Subduction Slab Western Canada, assuming a hypocentral depth of 50 km
    - High Branch
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGV, PGA, SA])

    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    GMPE_TABLE = os.path.join(os.path.dirname(__file__),
                              "gsc_canada_tables",
                              "WinslabD50_rhypo_high.hdf5")
