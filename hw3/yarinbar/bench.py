import _matrix
from random import random
import time

if __name__ == '__main__':

    A = _matrix.Matrix(1024, 1024)
    B = _matrix.Matrix(1024, 1024)

    for i in range(1024):
        for j in range(1024):
            num = random()
            A[i, j] = num
            B[i, j] = num

    start = time.time()
    _matrix.multiply_naive(A, B)
    end = time.time()
    naive_time = end - start

    start = time.time()
    _matrix.multiply_mkl(A, B)
    end = time.time()
    mkl_time = end-start

    speed_up = naive_time / mkl_time

    print('Start multiply_naive take min = {} seconds'.format(naive_time))
    print('Start multiply_mkl, take min = {} seconds'.format(mkl_time))
    print('MKL speed-up over naive: {} x'.format(speed_up))
