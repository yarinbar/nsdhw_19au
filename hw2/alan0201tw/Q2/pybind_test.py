import numpy as np
from vector import compute_angle

# import nose for unit-test

c_PI = np.pi
epsilon = np.finfo(np.float64).eps

################################################
##
##  Setting up
##
################################################

def setup_module(module):
    print ("") # this is to get a newline after the dots
    print ("========================================================================")
    print ("NOSE_SETUP : Setting up nosetest environment...")
    print ("NOSE_SETUP : using constant PI from numpy, PI = ", c_PI)
    print ("NOSE_SETUP : using constant epsilon for type 'double' from numpy, epsilon = ", epsilon)
    print ("========================================================================")

################################################
##
##  Helper functions
##
################################################

def expect_difference_within_epsolon(a, b):
    return ( np.abs(a - b) ) <= epsilon

def expect_difference_within_error(a, b, error):
    return ( np.abs(a - b) ) <= error

################################################
##
##  Unit tests
##
################################################

def test_zero_to_zero():
    assert np.isnan(compute_angle((0, 0), (0, 0)))

def test_unitX_to_unitNegX():
    assert expect_difference_within_epsolon(
        compute_angle((1, 0), (-1, 0)),
        c_PI)

def test_unitX_to_unitX():
    assert expect_difference_within_epsolon(
        compute_angle((1, 0), (1, 0)),
        0)

def test_unitX_to_unitY():
    assert expect_difference_within_epsolon(
        compute_angle((1, 0), (0, 1)),
        c_PI / 2)

def test_unitX_to_unitXY():
    assert expect_difference_within_epsolon(
        compute_angle((1, 0), (1, 1)),
        c_PI / 4)

# this helps find the bug when returning angles in radian might be negative
# when returning the value from atan2, it should have an additional abs()
def test_unitY_to_unitXY():
    assert expect_difference_within_epsolon(
        compute_angle((0, 1), (1, 1)),
        c_PI / 4)

def test_unitXY_to_unitNegXY():
    assert expect_difference_within_epsolon(
        compute_angle((1, 1), (-1, -1)),
        c_PI)

def test_identical_vector():
    assert expect_difference_within_epsolon(
        compute_angle((0.78221, 1.12313231), (0.78221, 1.12313231)),
        0)

def test_identical_vector_to_neg():
    assert expect_difference_within_epsolon(
        compute_angle((0.78221, 1.12313231), (-0.78221, -1.12313231)),
        c_PI)

################################################
##
##  Unit tests - following cases are from online 
##  tools for computing vector angles.
##
################################################

def test_random_test_case0():
    assert expect_difference_within_error(
        compute_angle((21, 1), (0, 13)),
        1.5232, 0.0001)

def test_random_test_case1():
    assert expect_difference_within_error(
        compute_angle((2.312, 2.1312), (-32.1232, 0.3123123)),
        2.3871, 0.0001)

def test_random_test_case2():
    assert expect_difference_within_error(
        compute_angle((1.23123, 7.1312321), (123.123, -4.1231321)),
        1.4333, 0.0001)

def test_random_test_case3():
    assert expect_difference_within_error(
        compute_angle((-9.17734, 7.5648962), (5.134146, -3.652341)),
        3.0705, 0.0001)