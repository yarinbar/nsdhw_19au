import pybind_vector2d as vd

import math
import pytest

# helper function
def dot(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def norm(v):
    return math.sqrt(dot(v, v))

def angle(v1, v2):
    return math.acos(dot(v1, v2) / (norm(v1) * norm(v2)))

# test angle from CPP
def compute_angle(v1: vd.Vector2D, v2: vd.Vector2D) -> float:
    return vd.angle(v1, v2)

def close_to(a: float, b: float) -> bool:
        return math.isclose(a, b, rel_tol=1e-5)

class TestClass:
    def test_dot(self):
        [v1, v2] = [vd.Vector2D(1, 2), vd.Vector2D(3, 4)]
        assert(close_to(v1 * v2,
                        dot([v1.getX(), v1.getY()], [v2.getX(), v2.getY()])))
        
    def test_norm(self):
        v1 = vd.Vector2D(1, 2)
        assert(close_to(v1.norm(), norm([v1.getX(), v1.getY()])))
        
    def test_angle(self):
        [v1, v2] = [vd.Vector2D(1, 2), vd.Vector2D(3, 4)]
        assert(close_to(compute_angle(v1, v2),
                        angle([v1.getX(), v1.getY()], [v2.getX(), v2.getY()])))
