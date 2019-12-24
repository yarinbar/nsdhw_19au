import mylib
import math

EPS = 10**-7

def testParallel():
    assert mylib.angle([1,0], [1,0]) == 0.0

def testRightAngle():
    assert abs(mylib.angle([1, 0], [0, 1]) - math.pi / 2) < EPS

def testOtherAngle():
    assert abs(mylib.angle([3, 4], [4, 3]) - 0.28379419445991516) < EPS

