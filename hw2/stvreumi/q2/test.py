from hw2q2 import calculate_angle
import unittest
import numpy as np

class hw2q2Test(unittest.TestCase):
    def test_0length(self):
        self.assertAlmostEqual(calculate_angle((0.0,0.0), (1.0,1.0)), -100.0)

    def test_0angle(self):
        self.assertAlmostEqual(calculate_angle((1.0,1.0), (1.0,1.0)), 0.0)
    
    def test_90angle(self):
        self.assertAlmostEqual(calculate_angle((1.0,0.0), (0.0,1.0)), np.pi/2.0)

if __name__ == '__main__':
    unittest.main()