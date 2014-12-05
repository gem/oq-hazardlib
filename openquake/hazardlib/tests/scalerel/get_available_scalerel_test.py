import unittest

from openquake.hazardlib import scalerel


class AvailableScaleRelTestCase(unittest.TestCase):

    def test_get_available_scalerel(self):
        self.assertEqual({'WC1994': scalerel.wc1994.WC1994,
                          'PeerMSR': scalerel.peer.PeerMSR,
                          'PointMSR': scalerel.point.PointMSR,
                          'CEUS2011': scalerel.ceus2011.CEUS2011, 
                          'StrasserInterface': scalerel.strasser2010.StrasserInterface,
                          'StrasserIntraslab': scalerel.strasser2010.StrasserIntraslab,
                          'BonillaEtAl1984Aspect0p5': scalerel.bonilla1984.BonillaEtAl1984Aspect0p5,
                          'BonillaEtAl1984Aspect1p0': scalerel.bonilla1984.BonillaEtAl1984Aspect1p0,
                          'BonillaEtAl1984Aspect1p5': scalerel.bonilla1984.BonillaEtAl1984Aspect1p5,
                          'BonillaEtAl1984Aspect2p0': scalerel.bonilla1984.BonillaEtAl1984Aspect2p0,
                          'BonillaEtAl1984Aspect2p5': scalerel.bonilla1984.BonillaEtAl1984Aspect2p5,
                          'BonillaEtAl1984Aspect3p0': scalerel.bonilla1984.BonillaEtAl1984Aspect3p0,
                          'BonillaEtAl1984Aspect3p5': scalerel.bonilla1984.BonillaEtAl1984Aspect3p5,
                          'BonillaEtAl1984Aspect4p0': scalerel.bonilla1984.BonillaEtAl1984Aspect4p0,
                          'BonillaEtAl1984Aspect5p0': scalerel.bonilla1984.BonillaEtAl1984Aspect5p0,
                          'BonillaEtAl1984Aspect6p0': scalerel.bonilla1984.BonillaEtAl1984Aspect6p0,
                          'BonillaEtAl1984Aspect7p0': scalerel.bonilla1984.BonillaEtAl1984Aspect7p0,
                          'BonillaEtAl1984Aspect8p0': scalerel.bonilla1984.BonillaEtAl1984Aspect8p0
                          },
                         dict(scalerel.get_available_scalerel()))

    def test_get_available_area_scalerel(self):
        self.assertEqual({'WC1994': scalerel.wc1994.WC1994,
                          'StrasserInterface': scalerel.strasser2010.StrasserInterface,
                          'StrasserIntraslab': scalerel.strasser2010.StrasserIntraslab,
                          'BonillaEtAl1984Aspect0p5': scalerel.bonilla1984.BonillaEtAl1984Aspect0p5,
                          'BonillaEtAl1984Aspect1p0': scalerel.bonilla1984.BonillaEtAl1984Aspect1p0,
                          'BonillaEtAl1984Aspect1p5': scalerel.bonilla1984.BonillaEtAl1984Aspect1p5,
                          'BonillaEtAl1984Aspect2p0': scalerel.bonilla1984.BonillaEtAl1984Aspect2p0,
                          'BonillaEtAl1984Aspect2p5': scalerel.bonilla1984.BonillaEtAl1984Aspect2p5,
                          'BonillaEtAl1984Aspect3p0': scalerel.bonilla1984.BonillaEtAl1984Aspect3p0,
                          'BonillaEtAl1984Aspect3p5': scalerel.bonilla1984.BonillaEtAl1984Aspect3p5,
                          'BonillaEtAl1984Aspect4p0': scalerel.bonilla1984.BonillaEtAl1984Aspect4p0,
                          'BonillaEtAl1984Aspect5p0': scalerel.bonilla1984.BonillaEtAl1984Aspect5p0,
                          'BonillaEtAl1984Aspect6p0': scalerel.bonilla1984.BonillaEtAl1984Aspect6p0,
                          'BonillaEtAl1984Aspect7p0': scalerel.bonilla1984.BonillaEtAl1984Aspect7p0,
                          'BonillaEtAl1984Aspect8p0': scalerel.bonilla1984.BonillaEtAl1984Aspect8p0},
                         dict(scalerel.get_available_area_scalerel()))

    def test_get_available_magnitude_scalerel(self):
        self.assertEqual({'PeerMSR': scalerel.peer.PeerMSR,
                          'PointMSR': scalerel.point.PointMSR,
                          'WC1994': scalerel.wc1994.WC1994,
                          'CEUS2011': scalerel.ceus2011.CEUS2011,
                          'StrasserInterface': scalerel.strasser2010.StrasserInterface,
                          'StrasserIntraslab': scalerel.strasser2010.StrasserIntraslab,
                          'BonillaEtAl1984Aspect0p5': scalerel.bonilla1984.BonillaEtAl1984Aspect0p5,
                          'BonillaEtAl1984Aspect1p0': scalerel.bonilla1984.BonillaEtAl1984Aspect1p0,
                          'BonillaEtAl1984Aspect1p5': scalerel.bonilla1984.BonillaEtAl1984Aspect1p5,
                          'BonillaEtAl1984Aspect2p0': scalerel.bonilla1984.BonillaEtAl1984Aspect2p0,
                          'BonillaEtAl1984Aspect2p5': scalerel.bonilla1984.BonillaEtAl1984Aspect2p5,
                          'BonillaEtAl1984Aspect3p0': scalerel.bonilla1984.BonillaEtAl1984Aspect3p0,
                          'BonillaEtAl1984Aspect3p5': scalerel.bonilla1984.BonillaEtAl1984Aspect3p5,
                          'BonillaEtAl1984Aspect4p0': scalerel.bonilla1984.BonillaEtAl1984Aspect4p0,
                          'BonillaEtAl1984Aspect5p0': scalerel.bonilla1984.BonillaEtAl1984Aspect5p0,
                          'BonillaEtAl1984Aspect6p0': scalerel.bonilla1984.BonillaEtAl1984Aspect6p0,
                          'BonillaEtAl1984Aspect7p0': scalerel.bonilla1984.BonillaEtAl1984Aspect7p0,
                          'BonillaEtAl1984Aspect8p0': scalerel.bonilla1984.BonillaEtAl1984Aspect8p0},
                         dict(scalerel.get_available_magnitude_scalerel()))

    def test_get_available_sigma_area_scalerel(self):
        self.assertEqual({'WC1994': scalerel.wc1994.WC1994,
                          'StrasserInterface': scalerel.strasser2010.StrasserInterface,
                          'StrasserIntraslab': scalerel.strasser2010.StrasserIntraslab,
                          'BonillaEtAl1984Aspect0p5': scalerel.bonilla1984.BonillaEtAl1984Aspect0p5,
                          'BonillaEtAl1984Aspect1p0': scalerel.bonilla1984.BonillaEtAl1984Aspect1p0,
                          'BonillaEtAl1984Aspect1p5': scalerel.bonilla1984.BonillaEtAl1984Aspect1p5,
                          'BonillaEtAl1984Aspect2p0': scalerel.bonilla1984.BonillaEtAl1984Aspect2p0,
                          'BonillaEtAl1984Aspect2p5': scalerel.bonilla1984.BonillaEtAl1984Aspect2p5,
                          'BonillaEtAl1984Aspect3p0': scalerel.bonilla1984.BonillaEtAl1984Aspect3p0,
                          'BonillaEtAl1984Aspect3p5': scalerel.bonilla1984.BonillaEtAl1984Aspect3p5,
                          'BonillaEtAl1984Aspect4p0': scalerel.bonilla1984.BonillaEtAl1984Aspect4p0,
                          'BonillaEtAl1984Aspect5p0': scalerel.bonilla1984.BonillaEtAl1984Aspect5p0,
                          'BonillaEtAl1984Aspect6p0': scalerel.bonilla1984.BonillaEtAl1984Aspect6p0,
                          'BonillaEtAl1984Aspect7p0': scalerel.bonilla1984.BonillaEtAl1984Aspect7p0,
                          'BonillaEtAl1984Aspect8p0': scalerel.bonilla1984.BonillaEtAl1984Aspect8p0},
                         dict(scalerel.get_available_sigma_area_scalerel()))

    def test_get_available_sigma_magnitude_scalerel(self):
        self.assertEqual({'PeerMSR': scalerel.peer.PeerMSR,
                          'WC1994': scalerel.wc1994.WC1994,
                          'StrasserInterface': scalerel.strasser2010.StrasserInterface,
                          'StrasserIntraslab': scalerel.strasser2010.StrasserIntraslab,
                          'BonillaEtAl1984Aspect0p5': scalerel.bonilla1984.BonillaEtAl1984Aspect0p5,
                          'BonillaEtAl1984Aspect1p0': scalerel.bonilla1984.BonillaEtAl1984Aspect1p0,
                          'BonillaEtAl1984Aspect1p5': scalerel.bonilla1984.BonillaEtAl1984Aspect1p5,
                          'BonillaEtAl1984Aspect2p0': scalerel.bonilla1984.BonillaEtAl1984Aspect2p0,
                          'BonillaEtAl1984Aspect2p5': scalerel.bonilla1984.BonillaEtAl1984Aspect2p5,
                          'BonillaEtAl1984Aspect3p0': scalerel.bonilla1984.BonillaEtAl1984Aspect3p0,
                          'BonillaEtAl1984Aspect3p5': scalerel.bonilla1984.BonillaEtAl1984Aspect3p5,
                          'BonillaEtAl1984Aspect4p0': scalerel.bonilla1984.BonillaEtAl1984Aspect4p0,
                          'BonillaEtAl1984Aspect5p0': scalerel.bonilla1984.BonillaEtAl1984Aspect5p0,
                          'BonillaEtAl1984Aspect6p0': scalerel.bonilla1984.BonillaEtAl1984Aspect6p0,
                          'BonillaEtAl1984Aspect7p0': scalerel.bonilla1984.BonillaEtAl1984Aspect7p0,
                          'BonillaEtAl1984Aspect8p0': scalerel.bonilla1984.BonillaEtAl1984Aspect8p0},
                         dict(scalerel.get_available_sigma_magnitude_scalerel()))
