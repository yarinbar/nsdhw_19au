import timeit

import numpy as np
import _matrix

def different_tile_size():
    matrix_size = 512
    A_content = np.eye(matrix_size)
    A = _matrix.Matrix(A_content)
    # print(A)

    B_content = np.arange(1, matrix_size ** 2 + 1).reshape(matrix_size, matrix_size)
    B = _matrix.Matrix(B_content)

    # print(B)

    MKL_MM = _matrix.multiply_mkl(A, B)

    for p in range(1, 9):
        tiled_MM = _matrix.multiply_tile(A, B, 2 ** p)
        assert tiled_MM == MKL_MM
    
    for p in [14, 22, 56]:
        tiled_MM = _matrix.multiply_tile(A, B, p)
        assert tiled_MM == MKL_MM

def performance():
    matrix_size = 512
    A_content = np.eye(matrix_size)
    A = _matrix.Matrix(A_content)
    # print(A)

    B_content = np.arange(1, matrix_size ** 2 + 1).reshape(matrix_size, matrix_size)
    B = _matrix.Matrix(B_content)


    for p in range(4, 8):
        tsize = 2 ** p
        ns = dict(_matrix=_matrix, mat1=A, mat2=B, tsize=tsize)
        t_tile = timeit.Timer('_matrix.multiply_tile(mat1, mat2, tsize)', globals=ns)
        t_naive = timeit.Timer('_matrix.multiply_naive(mat1, mat2)', globals=ns)

        time_tile = min(t_tile.repeat(10, 1))
        time_naive = min(t_naive.repeat(10, 1))
        ratio = time_tile/time_naive
        print(tsize, ratio)

def odd_tile_size():
    matrix_size = 3
    A_content = np.eye(matrix_size)
    A = _matrix.Matrix(A_content)
    # print(A)

    B_content = np.arange(1, matrix_size ** 2 + 1).reshape(matrix_size, matrix_size)
    B = _matrix.Matrix(B_content)

    print(_matrix.multiply_tile(A, B, 2))


def nd_array_test():
    size = 100
    mat = _matrix.Matrix(size,size)
    for it in range(size):
        for jt in range(size):
            mat[it, jt] = it * size + jt + 1
    np_array = np.array(mat, copy=False)
    np_array[0, 0] = 10
    print(np_array[0, 0], mat[0, 0])
    np_array[1, 0] = 200

    print(np_array[1, 0], mat[1, 0])
    np_array.fill(0)
    print(np_array[0, 0], mat[0, 0])

    mat.array.fill(-1)
    print(mat.array[0,0], mat[0, 0])


if __name__ == "__main__":
    nd_array_test()
