#!/bin/env python3.6
import unittest
import random, sys, math
from vec2 import vec2

"""
python implementation of vector computation
"""

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def equal_float(v1, v2):
    return abs(v1 - v2) < sys.float_info.epsilon

class VecTC(unittest.TestCase):

    def test_construction(self):

        with self.assertRaisesRegex(
            RuntimeError,
            "cannot not calculate angle of zero length vec2",
        ):
            vec2(0, 0).angle(vec2(12, 3))
            vec2(12, 3).angle(vec2(0, 0))
            vec2(0, 0).angle(vec2(0, 0))

    def test_constructor(self):
        """test constructor without any paramemter"""
        v = vec2()
        self.assertEqual(v.x(), 0)
        self.assertEqual(v.y(), 0)

        """test constructor with initial x and y value"""

        x, y = (random.random(), random.random())
        w = vec2(x, y)
        self.assertEqual(w.x(), x)
        self.assertEqual(w.y(), y)

    def test_setter_getter(self):
        """test setter"""
        v = vec2()
        x, y = (random.random(), random.random())
        v.x(x), v.y(y)
        self.assertEqual(v.x(), x)
        self.assertEqual(v.y(), y)

    def test_copy(self):
        """test copy"""
        v = vec2(random.random(), random.random())
        m = v.copy()
        self.assertEqual(m, v)
        self.assertNotEqual(repr(m), repr(v))
        m.x(m.x() + 1), m.y(m.y() + 1)
        self.assertNotEqual(m, v)

    def test_dot(self):
        """test dot"""
        v = vec2(random.random(), random.random())
        w = vec2(random.random(), random.random())
        self.assertTrue(
                equal_float(w.dot(v),
                            dotproduct([w.x(), w.y()], [v.x(), v.y()])))

    def test_length(self):
        """test length"""
        w = vec2(random.random(), random.random())
        self.assertTrue(
                equal_float(w.length(), length([w.x(), w.y()])))

    def test_angle(self):
        """test angle"""
        v = vec2(random.random(), random.random())
        w = vec2(random.random(), random.random())
        self.assertTrue(
                equal_float(
                    angle([v.x(), v.y()], [w.x(), w.y()]),
                    v.angle(w)))

    def test_operators(self):
        """test operators"""
        v = vec2(random.random(), random.random())
        w = vec2(random.random(), random.random())
        x, y = random.random(), random.random()
        self.assertTrue(vec2(x, y) == vec2(x, y))
        self.assertTrue(vec2(x, y) != vec2(x + 1, y + 1))
        self.assertNotEqual(vec2(x, y), vec2(x + 1, y + 1))

        self.assertEqual(v + w, vec2(v.x() + w.x(), v.y() + w.y()))
        self.assertEqual(v - w, vec2(v.x() - w.x(), v.y() - w.y()))
        self.assertEqual(+v, vec2(v.x(), v.y()))
        self.assertEqual(-v, vec2(-v.x(), -v.y()))

        m = random.random()
        self.assertEqual(v * m, vec2(v.x() * m, v.y() * m))
        self.assertEqual(m * v, vec2(m * v.x(), m * v.y()))

        m = v.copy()
        m += w
        self.assertEqual(m, v + w)

        m = v.copy()
        r = random.random()
        m *= r
        self.assertEqual(m, v * r)

if __name__ == '__main__':
    unittest.main()
