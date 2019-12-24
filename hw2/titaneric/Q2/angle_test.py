import unittest
import math
import re

from angle import getAngle

epsilon = 1e-6
pi = math.pi
pattern = r"test_(\d+)_deg"


def abs_diff(value, expected):
    return abs(value - expected)


def to_radian(func):
    m = re.match(pattern, func.__name__)
    return int(m.group(1)) / 180 * pi


class TestStringMethods(unittest.TestCase):

    def test_180_deg(self):
        vec1 = (1, 0)
        vec2 = (-1, 0)
        self.assertLessEqual(abs_diff(getAngle(vec1, vec2), to_radian(
            self.test_180_deg)), epsilon)

    def test_150_deg(self):
        vec1 = (1, 0)
        vec2 = (-math.sqrt(3), 1)

        self.assertLessEqual(
            abs_diff(getAngle(vec1, vec2), to_radian(self.test_150_deg)), epsilon)

    def test_135_deg(self):
        vec1 = (1, 0)
        vec2 = (-1, 1)

        self.assertLessEqual(
            abs_diff(getAngle(vec1, vec2), to_radian(self.test_135_deg)), epsilon)

    def test_120_deg(self):
        vec1 = (1, 0)
        vec2 = (-1, math.sqrt(3))
        self.assertLessEqual(
            abs_diff(getAngle(vec1, vec2), to_radian(self.test_120_deg)), epsilon)

    def test_90_deg(self):
        vec1 = (1, 0)
        vec2 = (0, 1)
        self.assertLessEqual(
            abs_diff(getAngle(vec1, vec2), to_radian(self.test_90_deg)), epsilon)

    def test_60_deg(self):
        vec1 = (1, 0)
        vec2 = (1, math.sqrt(3))
        self.assertLessEqual(
            abs_diff(getAngle(vec1, vec2), to_radian(self.test_60_deg)), epsilon)

    def test_45_deg(self):
        vec1 = (1, 0)
        vec2 = (1, 1)
        self.assertLessEqual(
            abs_diff(getAngle(vec1, vec2), to_radian(self.test_45_deg)), epsilon)

    def test_30_deg(self):
        vec1 = (1, 0)
        vec2 = (math.sqrt(3), 1)
        self.assertLessEqual(
            abs_diff(getAngle(vec1, vec2), to_radian(self.test_30_deg)), epsilon)

    def test_0_deg(self):
        vec1 = (1, 0)
        vec2 = (1, 0)
        self.assertLessEqual(
            abs_diff(getAngle(vec1, vec2), to_radian(self.test_0_deg)), epsilon)


if __name__ == '__main__':
    unittest.main()
