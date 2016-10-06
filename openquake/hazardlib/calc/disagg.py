import sys
import numpy
import copy

from openquake.baselib.python3compat import raise_

from openquake.hazardlib.calc import filters
from openquake.hazardlib.site import SiteCollection
from openquake.hazardlib.gsim.base import ContextMaker


class DsgMatrix():
    """
    A container for the information needed and the one produced during a
    disaggregation analysis. Some methods for updating information are also
    available.
    """
    def __init__(self, src_ids, trt_ids, mag_bins, dst_bins, eps_bins,
                 lon_bins, lat_bins):
        # Dictionary
        self.pdic = {}
        self.darr = None
        # BINS
        self.mag_bins = numpy.array(mag_bins)
        self.dst_bins = numpy.array(dst_bins)
        self.eps_bins = numpy.array(eps_bins)
        self.lon_bins = numpy.array(lon_bins)
        self.lat_bins = numpy.array(lat_bins)
        # SRC dictionary
        self.src_dict = {}
        for idx, key in enumerate(src_ids):
            self.src_dict[key] = idx
        # TRT dictionary
        self.trt_dict = {}
        for idx, key in enumerate(trt_ids):
            self.trt_dict[key] = idx

    def get_md_mtx(self):
        """
        This method creates the M-R disaggregation matrix
        """
        if self.darr is None:
            self._create_darr()

    def add_contribute(self, mag, dst, eps, trt_id, src_id, lon, lat, pnoe):
        """
        This method adds to the disaggregation dictionary the contribution
        provided by a rupture described in terms of the following parameters:
        :parameter mag:
        :parameter dst:
        :parameter eps:
        :parameter trt_id:
        :parameter src_id:
        :parameter lon:
        :parameter lat:
        """
        # Get indexes
        idx_mag = self._get_idx('mag', mag)
        idx_dst = self._get_idx('dst', dst)
        idx_eps = self._get_idx('eps', eps)
        idx_trt = self._get_trt_key(trt_id)
        idx_src = self._get_src_key(src_id)
        idx_lon = self._get_idx('lon', lon)
        idx_lat = self._get_idx('lat', lat)
        # Create the key for the current set of parameters
        key = ("%04d"*7) % (idx_mag, idx_dst, idx_eps, idx_trt, idx_src,
                            idx_lon, idx_lat)
        # Update the dictionary (note that for simplicity at the moment we
        # consider only the case of independent sources and ruptures).
        # TODO include groups and aggregate probabilities accordingly
        if key in self.pdic:
            self.pdic[key] = 1. - (1. - self.pdic[key]) * pnoe
        else:
            self.pdic[key] = 1. - pnoe

    def get_values_from_key(self, key):
        """
        This method returns the (binned) values of the parameters used to
        populate the disaggregation dictionary.

        :parameter key:
            A key in self.pdic.keys()
        :returns:
            A tuple with the following format:
            (mag, dst, eps, trt, src_id, lon, lat)
        """
        # Check that the key is used in the disaggregation dictionary
        assert key in set(self.pdic.keys())
        # Split the key into various parts. Each one identifies the index of
        # a bin in the histogram describing one property of the rupture
        idx = 0
        mag_str = key[idx:idx+4]
        idx += 4
        dst_str = key[idx:idx+4]
        idx += 4
        eps_str = key[idx:idx+4]
        idx += 4
        trt_str = key[idx:idx+4]
        idx += 4
        src_str = key[idx:idx+4]
        idx += 4
        lon_str = key[idx:idx+4]
        idx += 4
        lat_str = key[idx:idx+4]
        # Get values of the parameters from the various indexes
        mag_val = self._get_val('mag', mag_str)
        print mag_val, mag_str
        dst_val = self._get_val('dst', dst_str)
        print dst_val, dst_str
        eps_val = self._get_val('eps', eps_str)
        print eps_val, eps_str
        lon_val = self._get_val('lon', lon_str)
        print lon_val

        return

    def _get_src_key(self, src):
        if src not in self.src_dict:
            self.src_dict[src] = len(self.src_dict)
        return self.src_dict[src]

    def _get_trt_key(self, trt):
        if trt not in self.trt_dict:
            self.trt_dict[trt] = len(self.trt_dict)
        return self.trt_dict[trt]

    def _get_idx(self, tpe, val):
        """
        This covers dst, eps, mag, lat
        """
        if tpe == 'mag':
            vec = self.mag_bins
        elif tpe == 'dst':
            vec = self.dst_bins
        elif tpe == 'eps':
            vec = self.eps_bins
        elif tpe == 'lon':
            vec = self.lon_bins
        elif tpe == 'lat':
            vec = self.lat_bins
        else:
            raise ValueError('unknown type')

        if val < vec[0]:
            return 0
        elif val > vec[-1]:
            return len(vec) + 1
        else:
            return numpy.min(numpy.nonzero(vec > val)) + 1

    def _get_val(self, tpe, key):
        """
        This supports the calculation of the values of dst, eps, mag, lat given
        a key.
        """
        if tpe == 'mag':
            vec = self.mag_bins
        elif tpe == 'dst':
            vec = self.dst_bins
        elif tpe == 'eps':
            vec = self.eps_bins
        elif tpe == 'lon':
            vec = self.lon_bins
        elif tpe == 'lat':
            vec = self.lat_bins
        else:
            raise ValueError('unknown type')

        idx = int(key)-1
        if idx == len(vec):
            return vec[-1]+(vec[-1]-vec[-2])
        else:
            return vec[idx]


def dsgr(sources, site, imt, iml, gsims, dsgmtx, truncation_level,
         n_epsilons,
         source_site_filter='SourceSitesFilter',
         maximum_distance=None,):
    """
    Disaggregation calculator
    """

    ponet = 1.0
    dsgm = copy.deepcopy(dsgmtx)

    sites = SiteCollection([site])
    sitemesh = sites.mesh

    if source_site_filter == 'SourceSitesFilter':  # default
        source_site_filter = (
            filters.SourceSitesFilter(maximum_distance)
            if maximum_distance else filters.source_site_noop_filter)

    _next_trt_num = 0
    trt_nums = {}
    # here we ignore filtered site collection because either it is the same
    # as the original one (with one site), or the source/rupture is filtered
    # out and doesn't show up in the filter's output
    for src_idx, (source, s_sites) in \
            enumerate(source_site_filter(sources, sites)):
        try:

            tect_reg = source.tectonic_region_type
            gsim = gsims[tect_reg]
            cmaker = ContextMaker([gsim])

            for rupture in source.iter_ruptures():
                # extract rupture parameters of interest
                # mags.append(rupture.mag)
                [jb_dist] = rupture.surface.get_joyner_boore_distance(sitemesh)
                [closest_point] = rupture.surface.get_closest_points(sitemesh)

                src_id = source.source_id
                mag = rupture.mag
                dst = jb_dist
                lon = closest_point.longitude
                lat = closest_point.latitude
                limits = numpy.linspace(-truncation_level,
                                        truncation_level, n_epsilons + 1)
                centers = (limits[1:]+limits[:-1]) /2

                # compute conditional probability of exceeding iml given
                # the current rupture, and different epsilon level, that is
                # ``P(IMT >= iml | rup, epsilon_bin)`` for each of epsilon bins
                sctx, rctx, dctx = cmaker.make_contexts(s_sites, rupture)

                [poes_given_rup_eps] = gsim.disaggregate_poe(
                    sctx, rctx, dctx, imt, iml, truncation_level, n_epsilons
                )

                # collect probability of a rupture causing no exceedances
                pnoes = rupture.get_probability_no_exceedance(poes_given_rup_eps)

                for eps, pnoe in zip(centers, pnoes):
                    if pnoe < (1-1.e-8):
                        dsgm.add_contribute(mag, dst, eps, tect_reg, src_id,
                                            lon, lat, pnoe)
                        # Total probability of non-exceedance
                        ponet *= pnoe

        except Exception as err:
            etype, err, tb = sys.exc_info()
            msg = 'An error occurred with source id=%s. Error: %s'
            msg %= (source.source_id, str(err))
            raise_(etype, msg, tb)

    print 'poet:', 1.-ponet
    return dsgm


"""
def main():
    src_ids = [1]
    trt_ids = [TRT.ACTIVE_SHALLOW_CRUST]
    mag_bins = [4.0, 4.1, 4.2]
    dst_bins = [10., 15., 20., 25., 30.]
    eps_bins = [-3, 0, 3]
    lon_bins = [10.0, 10.1, 10.2]
    lat_bins = [45.1, 45.2, 45.3]

    dm = DsgMatrix(src_ids, trt_ids, mag_bins, dst_bins, eps_bins,
                   lon_bins, lat_bins)

    idx = dm.get_idx('mag', 3.0)
    assert idx == 0
    idx = dm.get_idx('mag', 4.05)
    assert idx == 1
    idx = dm.get_idx('mag', 4.4)
    assert idx == 3

    idx = dm.get_idx('dst', 5.0)
    assert idx == 0
    idx = dm.get_idx('dst', 12.0)
    assert idx == 1
    idx = dm.get_idx('dst', 60.0)
    assert idx == 5

if __name__ == "__main__":
    main()
"""
