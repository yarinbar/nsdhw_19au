import random
import _matrix
import unittest
import numpy as np
import sys

def equal_float(v1, v2):
    if not abs(v1 - v2) < sys.float_info.epsilon:
        print(abs(v1 - v2), sys.float_info.epsilon)
    return abs(v1 - v2) < sys.float_info.epsilon

class MatrixTest(unittest.TestCase):

    def make_matrices(self, row, col):

        matrix = _matrix.Matrix(row, col)
        darray = np.zeros((row, col))
        for i in range(row):
            for j in range(col):
                v = random.randint(0, 10000)
                matrix[i, j] = v
                darray[i, j] = v
        return matrix, darray

    def test_basic(self):

        row, col = random.randrange(50, 200), random.randrange(50, 200)
        m, a = self.make_matrices(row, col)
        self.assertEqual((m.nrow, m.ncol), a.shape)
        self.assertEqual((row, col), (m.nrow, m.ncol))
        for i in range(row):
            for j in range(col):
                self.assertEqual(m[i, j], a[i, j])

    def test_multiply(self):

        d1, d2, d3 = random.randrange(50, 200), \
                     random.randrange(50, 200), \
                     random.randrange(50, 200)

        d1, d2, d3 = random.randrange(3, 4), \
                     random.randrange(3, 4), \
                     random.randrange(3, 4)

        m, a = self.make_matrices(d1, d2)
        n, b = self.make_matrices(d2, d3)

        k = _matrix.multiply_tile(m, n, 64)
        l = _matrix.multiply_mkl(m, n)
        c = np.dot(a, b)
        for i in range(d1):
            for j in range(d3):
                self.assertEqual(k[i, j], l[i, j])
                self.assertEqual(l[i, j], c[i, j])

# vim: set fenc=utf8 ff=unix et sw=4 ts=4 sts=4:
