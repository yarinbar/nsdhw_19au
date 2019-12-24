from angle import *
import pytest
import math

PI = 3.14159265358979323846

def test_same_point():
    p1 = Point(0, 0)
    p2 = Point(0, 0)
    r_angle = calc_angle(p1, p2)
    assert math.isclose(r_angle, -1, abs_tol=1e-06)

def test_0_degree():
    p1 = Point(0, 1)
    p2 = Point(0, 2)
    r_angle = calc_angle(p1, p2)
    assert math.isclose(r_angle, 0, abs_tol=1e-06)

def test_45_degree():
    p1 = Point(1, 1)
    p2 = Point(0, 1)
    r_angle = calc_angle(p1, p2)
    assert math.isclose(r_angle, PI/4, abs_tol=1e-06)

def test_90_degree():
    p1 = Point(1, 0)
    p2 = Point(0, 1)
    r_angle = calc_angle(p1, p2)
    assert math.isclose(r_angle, PI/2, abs_tol=1e-06)

def test_135_degree():
    p1 = Point(1, 0)
    p2 = Point(-1, 1)
    r_angle = calc_angle(p1, p2)
    assert math.isclose(r_angle, 3*PI/4, abs_tol=1e-06)

def test_180_degree():
    p1 = Point(1, 0)
    p2 = Point(-1, 0)
    r_angle = calc_angle(p1, p2)
    assert math.isclose(r_angle, PI, abs_tol=1e-06)