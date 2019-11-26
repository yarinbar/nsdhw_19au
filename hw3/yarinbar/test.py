from unittest import TestCase
import _matrix
import random

class TestMat(TestCase):

    def create_random(self, rows, cols):
        A = _matrix.Matrix(rows, cols)
        for i in range(rows):
            for j in range(cols):
                A[i, j] = random.randint(0, 100)
        return A

    def create_zero(self, rows, cols):
        A = _matrix.Matrix(rows, cols)
        for i in range(rows):
            for j in range(cols):
                A[i, j] = 0
        return A

    def create_id(self, rows, cols):
        A = _matrix.Matrix(rows, cols)
        for i in range(rows):
            for j in range(cols):
                if i == j:
                    A[i, j] = 1
                else:
                    A[i, j] = 0
        return A

    def test_matrix_size(self):
        for i in range(10, 200, 10):
            m = self.create_zero(i, i)
            self.assertEqual(m.nrow, i)
            self.assertEqual(m.ncol, i)

    def test_naive_mul(self):
        zero = self.create_zero(10, 10)
        id   = self.create_id(10, 10)
        rand = self.create_random(10, 10)
        self.assertEqual(_matrix.multiply_naive(zero, id), zero)
        self.assertEqual(_matrix.multiply_naive(id, zero), zero)
        self.assertEqual(_matrix.multiply_naive(id, id), id)
        self.assertEqual(_matrix.multiply_naive(rand, id), rand)
        self.assertEqual(_matrix.multiply_naive(id, rand), rand)

    def test_mkl_mul(self):
        zero = self.create_zero(10, 10)
        id   = self.create_id(10, 10)
        rand = self.create_random(10, 10)

        self.assertEqual(_matrix.multiply_mkl(zero, id), zero)
        self.assertEqual(_matrix.multiply_mkl(id, zero), zero)
        self.assertEqual(_matrix.multiply_mkl(id, id), id)
        self.assertEqual(_matrix.multiply_mkl(rand, id), rand)
        self.assertEqual(_matrix.multiply_mkl(id, rand), rand)


    def test_equal(self):
        A = self.create_random(5, 10)
        B = self.create_random(10, 20)

        self.assertEqual(_matrix.multiply_mkl(A, B), _matrix.multiply_naive(A, B))
