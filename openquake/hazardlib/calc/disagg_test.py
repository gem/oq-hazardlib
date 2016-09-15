
import unittest
from disagg import DsgMatrix
from openquake.hazardlib.const import TRT


class TestDsgMtx(unittest.TestCase):

    def setUp(self):
        src_ids = [1]
        trt_ids = [TRT.ACTIVE_SHALLOW_CRUST]
        mag_bins = [4.0, 4.1, 4.2]
        dst_bins = [10., 15., 20., 25., 30.]
        eps_bins = [-3, 0, 3]
        lon_bins = [10.0, 10.1, 10.2]
        lat_bins = [45.1, 45.2, 45.3]
        self.dm = DsgMatrix(src_ids, trt_ids, mag_bins, dst_bins, eps_bins,
                            lon_bins, lat_bins)

    def test_mag_idx_calculation(self):
        idx = self.dm._get_idx('mag', 3.0)
        assert idx == 0
        idx = self.dm._get_idx('mag', 4.05)
        assert idx == 1
        idx = self.dm._get_idx('mag', 4.4)
        assert idx == 3

    def test_dst_idx_calculation(self):
        idx = self.dm._get_idx('dst', 5.0)
        assert idx == 0
        idx = self.dm._get_idx('dst', 12.0)
        assert idx == 1
        idx = self.dm._get_idx('dst', 60.0)
        assert idx == 6

    def test_dst_eps_calculation(self):
        idx = self.dm._get_idx('eps', 2.0)
        assert idx == 1


    def test_add_contribute(self):
        mag = 4.05
        dst = 12.0
        eps = 2.0
        trt_id = TRT.ACTIVE_SHALLOW_CRUST
        src_id = '1'
        pnoe = 0.01

        key = '00010001000100000001'
        self.dm.add_contribute(mag, dst, eps, trt_id, src_id, pnoe)
        assert self.dm.pdic[key] == (1. - pnoe)

        pnoe = 0.03
        self.dm.add_contribute(mag, dst, eps, trt_id, src_id, pnoe)
        prb = 1. - (0.01 * 0.03)
        assert self.dm.pdic[key] == prb

        key = '00000001000100000001'
        mag = 3.00
        self.dm.add_contribute(mag, dst, eps, trt_id, src_id, pnoe)
        assert self.dm.pdic[key] == (1. - pnoe)
