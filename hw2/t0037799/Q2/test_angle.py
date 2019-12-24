import unittest
import angle
import numpy as np

class TestDegree(unittest.TestCase):
    def setUp(self):
        self.eps = 1e-6

    def test_zero_length(self):
        self.assertTrue(np.isnan(angle.angle((1,0),(0,0))))
        self.assertTrue(np.isnan(angle.angle((0,0),(1,0))))
        self.assertTrue(np.isnan(angle.angle((0,0),(0,0))))

    def test_0_degree(self):
        self.assertLessEqual(abs(angle.angle((1,0),(1,0)) - 0), self.eps)

    def test_90_degree(self):
        self.assertLessEqual(abs(angle.angle((1,0),(0,1)) - np.pi/2), self.eps)

    def test_45_degree(self):
        self.assertLessEqual(abs(angle.angle((1,0),(1,1)) - np.pi/4), self.eps)

if __name__ == '__main__':
    unittest.main()
