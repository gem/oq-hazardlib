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

"""
Module :mod:`openquake.hazardlib.gmice.base` defines base classes for
different kinds of :class:`ground motion intensity conversion equations
<GMICE>`.
"""
import abc
import numpy as np
from openquake.hazardlib.imt import MMI
from openquake.hazardlib.gsim.base import GMPE, IPE

class GMICE(GMPE):
    """
    Base class for a type of model known as a Ground Motion Intensity
    Conversion Equation (GMIC). These models represent imt-specific
    equations to convert between macroseismic intensity and known
    ground motion parameters.

    The majority of GMICEs are derived from regression models that are applied
    in both directions:
    GM -> Macroseismic Intensity
    Macroseismic Intensity -> GM

    The base class contains several key methods to apply the conversions.
    """
    # Defines the range of intensity values for which the GMICE applies
    INTENSITY_RANGE = (0.0, 12.0) # Defaulting to full range

    # Shakemap requires a "pretty" source name (e.g. Worden et al., (2012))
    SOURCE_NAME = ""

    # Shakemap requires a scale
    SCALE = ""
    
    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        Doesn't function in this case but is an abstract method of the GMPE
        class - so this is over-ridden to raise a not implemented error
        """
        raise NotImplementedError("Method 'get_mean_and_stddevs' not "
                                  "available in GMICE object")

    @abc.abstractmethod
    def get_mean_intensity_and_stddevs(self, imls, sites, rup, dists, imt,
                                       stddev_types):
        """
        Convert to macroseismic intensity from engineering ground motion
        parameters

        Method must be implemented by subclasses

        :param imls:
            Vector (or array) of ground motion intensity measure levels.
        :param sites:
            Instance of :class:`openquake.hazardlib.gsim.base.SitesContext`
            with parameters of sites collection assigned to respective values
            as numpy arrays. Only those attributes that are listed in class'
            :attr:`openquake.hazardlib.gsim.base.GroundShakingIntensityModel.REQUIRES_SITES_PARAMETERS`
            set are available.
        :param rup:
            Instance of :class:`openquake.hazardlib.gsim.base.RuptureContext`
            with parameters of a rupture assigned to respective values. Only
            those attributes that are listed in class'
            :attr:`openquake.hazardlib.gsim.base.GroundShakingIntensityModel.REQUIRES_RUPTURE_PARAMETERS`
            set are available.
        :param dists:
            Instance of :class:`openquake.hazardlib.gsim.base.DistancesContext`
            with values of distance measures between the rupture and each site
            of the collection assigned to respective values as numpy arrays.
            Only those attributes that are listed in class'
            :attr:`openquake.hazardlib.gsim.base.GroundShakingIntensityModel.REQUIRES_DISTANCES` set are
            available.
        :param imt:
            An instance (not a class) of intensity measure type.
            See :mod:`openquake.hazardlib.imt`.
        :param stddev_types:
            List of standard deviation types, constants from
            :class:`openquake.hazardlib.const.StdDev`.
            Method result value should include
            standard deviation values for each of types in this list.
        :returns:
            Method should return a tuple of two items. The first item is
            the numpy array of floats representing the output (converted)
            macroseismic intensity values. The second item is a list of
            numpy arrays corresponding to each of the standard deviation
            types
        """

    @abc.abstractmethod
    def get_mean_gm_and_stddevs(self, imls, sites, rup, dists, imt,
                                stddev_types):
        """
        Convert from macroseismic intensity to engineering ground motion
        parameters

        :param imls:
            Vector (or array) of macroseismic intensity measure levels.
        :returns:
            Method should return a tuple of two items. The first item is
            the numpy array of floats representing the output (converted)
            ground motion values. The second item is a list of
            numpy arrays corresponding to each of the standard deviation
            types
        """

    def get_mean_intensity_and_stddevs_from_gmpe(self, gmpe, sites, rup, dists,
                                                 imt, stddev_types):
        """
        Same as the ordinary method but fed with the GMPE instead 

        :param gmpe:
            Ground motion prediction equation as subclass of :class:
            `openquake.hazardlib.gsim.base.GMPE`
        """
        imls, _ = gmpe.get_mean_and_stddevs(sites, rup, dists, imt,
                                            stddev_types)
        return self.get_mean_intensity_and_stddevs(np.exp(imls), sites,
                                                   rup, dists,
                                                   imt, stddev_types)

    def get_mean_gm_and_stddevs_from_ipe(self, ipe, sites, rup, dists, imt,
                                         stddev_types, mmi_imt=MMI()):
        """
        Same as the ordinary method but fed with the intensity predicution
        equation instead

        :param ipe:
            Macroseismic intensity prediction equation as subclass of :class:
            `openquake.hazardlib.gsim.base.IPE`
        """
        # Verify that in this case the IMT is macroseismic intensity
        assert isinstance(ipe, IPE)
        imls, _ = ipe.get_mean_and_stddevs(sites, rup, dists, mmi_imt,
                                           stddev_types)
        return self.get_mean_gm_and_stddevs(imls, sites, rup, dists, imt,
                                            stddev_types)
