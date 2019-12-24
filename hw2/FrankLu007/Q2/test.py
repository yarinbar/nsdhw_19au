import math
import libcpp

eps = 10**-6

def isParallel():
    assert libcpp.angle([1,0], [1,0]) == 0.0

def testAngle():
    assert abs(libcpp.angle([1, 0], [0, 1]) - 3.1415926 / 2) < eps

def testAngle2():
    assert abs(libcpp.angle([3, 4], [4, 3]) - 0.2837941) < eps