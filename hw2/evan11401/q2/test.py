import unittest
import math
from angle import calculateAngle

class AngleUnitTest(unittest.TestCase):
    epsilon = 1e-8
    def test_degree0(self):
        vec1 = (0.,1.)
        vec2 = (0.,1.)
        test_angle = calculateAngle(vec1,vec2)
        true_angle = math.radians(0)
        print(true_angle,test_angle)
        self.assertLess(abs(test_angle-true_angle),self.epsilon)
    def test_degree30(self):
        vec1 = (1.,0.)
        vec2 = (math.sqrt(3),1.)
        test_angle = calculateAngle(vec1,vec2)
        true_angle = math.radians(30)
        print(true_angle,test_angle)
        self.assertLess(abs(test_angle-true_angle),self.epsilon)
    def test_degree45(self):
        vec1 = (1.,0.)
        vec2 = (math.sqrt(2),math.sqrt(2))
        test_angle = calculateAngle(vec1,vec2)
        true_angle = math.radians(45)
        print(true_angle,test_angle)
        self.assertLess(abs(test_angle-true_angle),self.epsilon)
    def test_degree60(self):
        vec1 = (1.,0.)
        vec2 = (1.,math.sqrt(3))
        test_angle = calculateAngle(vec1,vec2)
        true_angle = math.radians(60)
        print(true_angle,test_angle)
        self.assertLess(abs(test_angle-true_angle),self.epsilon)
    def test_degree90(self):
        vec1 = (1.,0.)
        vec2 = (0.,1.)
        test_angle = calculateAngle(vec1,vec2)
        true_angle = math.radians(90)
        print(true_angle,test_angle)
        self.assertLess(abs(test_angle-true_angle),self.epsilon)
    def test_degree180(self):
        vec1 = (0.,1.)
        vec2 = (0.,-1.)
        test_angle = calculateAngle(vec1,vec2)
        true_angle = math.radians(180)
        print(true_angle,test_angle)
        self.assertLess(abs(test_angle-true_angle),self.epsilon)
    def test_degree270(self):
        vec1 = (0.,1.)
        vec2 = (1.,0.)
        test_angle = calculateAngle(vec1,vec2)
        true_angle = math.radians(270)
        print(true_angle,test_angle)
        self.assertLess(abs(test_angle-true_angle),self.epsilon)

if __name__ == '__main__':
    unittest.main()
