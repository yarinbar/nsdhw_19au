from PointVec import PointVec, angle_btw
import pytest
import math

MAGIC_X, MAGIC_Y = 27734.23, 91823.99
MAGIC_PT = (MAGIC_X, MAGIC_Y)
MAGIC_SCALE = 98.038977663

def float_equal( v1:float, v2:float, abs_tol=1e-05):
  return math.isclose(v1,v2, abs_tol=abs_tol)

def scale_and_rotate( xy, times, rad):
  x, y = xy
  x, y = x*times, y*times
  nx =  x * math.cos(rad) + y * math.sin(rad)
  ny = -x * math.sin(rad) + y * math.cos(rad)
  return (nx,ny)

def test_ptvec_init():
  x,y = MAGIC_PT
  p1 = PointVec(x,y)
  assert float_equal(x,p1.m_x)
  assert float_equal(y,p1.m_y)

  p2 = PointVec()
  assert float_equal(0,p2.m_x)
  assert float_equal(0,p2.m_y)

def test_ptvec_angle_btw_180():
  rad = math.pi # 180 deg
  assert float_equal(rad, angle_btw(PointVec(1,0), PointVec(-1,0)))

  x,y = MAGIC_PT
  u,v = scale_and_rotate( (x,y), MAGIC_SCALE, rad)
  p1 = PointVec(x,y)
  p2 = PointVec(u,v)
  assert float_equal(rad, angle_btw(p1,p2))

def test_ptvec_angle_btw_0():
  rad = 0
  assert float_equal(rad, angle_btw(PointVec(1,0), PointVec(1,0)))

  x,y = MAGIC_PT
  u,v = scale_and_rotate( (x,y), MAGIC_SCALE, rad)
  p1 = PointVec(x, y)
  p2 = PointVec(u,v)
  assert float_equal(rad, angle_btw(p1,p2))

def test_ptvec_angle_btw_90():
  rad = math.pi/2 # 90 deg
  assert float_equal(rad, angle_btw(PointVec(-1,0), PointVec(0,-1)))

  x,y = MAGIC_PT
  u,v = scale_and_rotate( (x,y), MAGIC_SCALE, rad)
  p1 = PointVec(x,y)
  p2 = PointVec(u,v)
  assert float_equal( rad, angle_btw(p1,p2) )

def test_ptvec_angle_btw_magic():
  rad = math.pi * (3.78**9)
  rad = rad % (2*math.pi)
  x,y = MAGIC_PT
  u,v = scale_and_rotate( (x,y), MAGIC_SCALE, rad)
  p1 = PointVec(x,y)
  p2 = PointVec(u,v)
  assert float_equal( rad, angle_btw(p1,p2) )