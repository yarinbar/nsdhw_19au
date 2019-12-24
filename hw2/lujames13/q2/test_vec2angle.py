import unittest
from vec2angle import vec2angle
import math


class TestVec2angle(unittest.TestCase):

	def test_zero_length_2_vector(self):
		self.assertEqual(1, math.isnan(vec2angle((0,0),(0,0))))

	def test_zero_angle(self):
		self.assertAlmostEqual(vec2angle((1,0),(3,0)), 0.0)

	def test_right_angle(self):
		self.assertAlmostEqual(vec2angle((1,0),(0,-1)), math.pi/2)

	def test_one_other_angle(self):
		norm = (3**2 + 4**2)**(0.5) * (7**2 + (-2)**2)**(0.5)
		dot = 3*7 + 4*(-2)
		angle = math.acos(dot / norm)
		self.assertAlmostEqual(vec2angle((3,4),(7,-2)), angle)
	