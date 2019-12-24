import time
import _matrix
import random

SIZE = 1000
REPEAT = 5
FILE_NAME = "performance.txt"

def benchmark(func, A, B, *args):
    start = time.time()
    C = func(A, B, *args)
    end = time.time()
    return end - start

if __name__ == "__main__":
    # set up matrix
    mat1 = _matrix.Matrix(SIZE, SIZE)
    mat2 = _matrix.Matrix(SIZE, SIZE)
    for i in range(SIZE):
        for j in range(SIZE):
            mat1[i,j] = random.uniform(0, 1000)
            mat2[i,j] = random.uniform(0, 1000)
    # start bench
    min_mkl_time = float('inf')
    min_tile_time = float('inf')
    for i in range(REPEAT):
        temp = benchmark(_matrix.multiply_mkl, mat1, mat2)
        if min_mkl_time > temp: min_mkl_time = temp
        temp = benchmark(_matrix.multiply_tile, mat1, mat2, 17)
        if min_tile_time > temp: min_tile_time = temp
    # output to file
    with open(FILE_NAME, "w") as f:
        f.write(
"""Start multiply_mkl (repeat={}), take min = {} seconds
Start multiply_tile (repeat={}), take min = {} seconds
tile speed-up over mkl: {} x
""".format(REPEAT, min_mkl_time, REPEAT, min_tile_time, min_mkl_time / min_tile_time))
