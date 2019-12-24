import unittest

class Angle2DTestCase(unittest.TestCase):
    '''2D Cartesian coordinate angle test cases'''

    def setUp(self):
        import coord
        import math
        self.angle = coord.angle
        self.pi = math.pi
        self.eps = 1e-7

    def test_angle_zero(self):
        cases = [
            (1, 1),
            (1, -1),
            (-1, 1),
            (-1, -1)
        ]
        for v in cases:
            self.assertTrue(self.angle(v, v) < self.eps)

    def test_angle_pi(self):
        cases = [
            ((1, 1), (2, -2)),
            ((1, -1), (3, 3)),
            ((-1, 1), (1, 1)),
            ((-1, -1), (2, -2))
        ]
        for v1, v2 in cases:
            self.assertTrue(self.angle(v1, v2) - self.pi < self.eps)

    def test_commutivity(self):
        cases = [
            ((1, 2), (3, 4)),
            ((-1, 2), (-3, 4)),
            ((1, -2), (-3, 4)),
            ((-1, 2), (3, -4)),
            ((1, -2), (3, -4)),
            ((-1, -2), (-3, -4)),
        ]
        for v1, v2 in cases:
            self.assertTrue(self.angle(v1, v2) - self.angle(v2, v1) < self.eps)

if __name__ == '__main__':
    unittest.main()
