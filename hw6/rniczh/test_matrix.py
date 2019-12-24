import unittest
import random
import _matrix

class Test(unittest.TestCase):
    def create_matrix(self, m, n):
        # setup matrix
        A = _matrix.Matrix(m, n)
        for i in range(m):
            for j in range(n):
                A[i, j] = random.randint(0, 1000)
        return A

    def test_size(self):
        for i in range(10, 100, 10):
            m = self.create_matrix(i, i)
            self.assertEqual(m.nrow, i)
            self.assertEqual(m.ncol, i)

    def test_naive_mul(self):
        m1 = self.create_matrix(2, 3)
        m2 = self.create_matrix(3, 4)
        c1 = _matrix.multiply_naive(m1, m2)
        self.assertEqual(c1.nrow, 2)
        self.assertEqual(c1.ncol, 4)
        for i in range(2):
            for j in range(4):
                v = 0.0
                for k in range(3):
                    v += m1[i,k] * m2[k,j]
                self.assertEqual(v, c1[i,j])
        
    def test_mkl_mul(self):
        m1 = self.create_matrix(2, 3)
        m2 = self.create_matrix(3, 4)
        c2 = _matrix.multiply_mkl(m1, m2)
        self.assertEqual(c2.nrow, 2)
        self.assertEqual(c2.ncol, 4)
        for i in range(2):
            for j in range(4):
                v = 0.0
                for k in range(3):
                    v += m1[i,k] * m2[k,j]
                self.assertEqual(v, c2[i,j])

    def test_tile_mul(self):
        m1 = self.create_matrix(2, 3)
        m2 = self.create_matrix(3, 4)
        c3 = _matrix.multiply_tile(m1, m2, 64)
        self.assertEqual(c3.nrow, 2)
        self.assertEqual(c3.ncol, 4)
        for i in range(2):
            for j in range(4):
                v = 0.0
                for k in range(3):
                    v += m1[i,k] * m2[k,j]
                self.assertEqual(v, c3[i,j])                
        
    def test_equal(self):
        m1 = self.create_matrix(3, 4)
        m2 = self.create_matrix(4, 5)
        c1 = _matrix.multiply_naive(m1, m2)
        c2 = _matrix.multiply_mkl(m1, m2)
        c3 = _matrix.multiply_tile(m1, m2, 64)
        self.assertEqual(c1, c2)
        self.assertEqual(c2, c3)
