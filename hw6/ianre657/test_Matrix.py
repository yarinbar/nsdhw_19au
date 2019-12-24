from _matrix import Matrix
from _matrix import multiply_naive, multiply_mkl
import pytest

def compareMat(mat1, mat2, m, n):
    for i in range(m):
        for j in range(n):
            if mat1[i,j]!= mat2[i,j]:
                print(f'm1[i,j]:{mat1[i,j]}')
                return False
    return True

def test_constructor():
    m,n =20,10
    m1 = Matrix(m,n)
    m2 = Matrix(m,n)
    for i in range(m):
        for j in range(n):
            m1[i,j] = 8.0
            m2[i,j] = 8.0
    assert(compareMat(m1,m2,m,n))
    assert(m1==m2)

def test_mat_info():
    r,c = 13,19
    m1 = Matrix(r,c)
    assert(r==m1.nrow)
    assert(c==m1.ncol)

def test_multiply():
    m1 = Matrix(10,20)
    m2 = Matrix(20,10)

    for i in range(10):
        for j in range(20):
            m1[i,j] = 8.0

    for i in range(20):
        for j in range(10):
            m2[i,j] = 6.0

    m3 = multiply_naive(m1,m2)
    m4 = multiply_mkl(m1,m2)
    assert(compareMat(m3,m4,10,10))
