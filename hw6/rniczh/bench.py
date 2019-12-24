import time
import random
import _matrix

def benchmark(mul, A, B):
    start = time.time()
    C = mul(A, B)
    end = time.time()
    return (C, end - start)

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

if __name__ == "__main__":
    # benchmark configuration
    # n:      the matrix with n * n
    # tsize:  the tiling size of tile version 
    # repeat: number of testing, take the `min` executation time
    n = 1000
    tsize = 64
    repeat = 5

    # setup matrix
    print("setup random matrix with {} * {}...".format(n, n))
    A = _matrix.Matrix(n, n)
    B = _matrix.Matrix(n, n)
    for i in range(n):
        for j in range(n):
            A[i, j] = random.randint(0, 1000)
            B[i, j] = random.randint(0, 1000)

    # doing benchmark
    print("doing benchmark ... repeat {} times".format(repeat))

    min_naive_time = float("inf")
    min_dgemm_time = float("inf")
    min_tile_time  = float("inf")
    for i in range(repeat):
        # naive
        (C1, naive_time) = benchmark(_matrix.multiply_naive, A, B)
        if naive_time < min_naive_time:
            min_naive_time = naive_time
        # mkl
        (C2, dgemm_time) = benchmark(_matrix.multiply_mkl, A, B)
        if dgemm_time < min_dgemm_time:
            min_dgemm_time = dgemm_time
        # tile
        start_tile = time.time()
        C3 = _matrix.multiply_tile(A, B, tsize)
        end_tile = time.time()
        tile_time = end_tile - start_tile
        if tile_time < min_tile_time:
            min_tile_time = tile_time
        
    filename = "performance.txt"
    f = open(filename,"w+")
    print("writing to file: {}".format(bcolors.UNDERLINE + filename + bcolors.ENDC))
    f.write(
"""Start multiply_naive (repeat={}), take min = {} seconds
Start multiply_mkl (repeat={}), take min = {} seconds
Start multiply_tile (repeat={}), take min = {} seconds
MKL speed-up over naive: {} x
MKL speed-up over tile: {} x
tile speed-up over naive: {} x
""".format(repeat, min_naive_time,
           repeat, min_dgemm_time,
           repeat, min_tile_time,
           min_naive_time/min_dgemm_time,
           min_tile_time/min_dgemm_time,
           min_naive_time/min_tile_time))
    f.close()

    ok_str = bcolors.OKGREEN + "Successful" + bcolors.ENDC
    fail_str = bcolors.FAIL + "Failed" + bcolors.ENDC
    equality = ok_str if C1 == C2 and C2 == C3 else fail_str

    print(
"""
================== benchmark =================
- Multiply {} * {} matrix: 
  [Naive] {:.4f} seconds
  [Dgemm] {:.4f} seconds
  [Tile]  {:.4f} seconds

- Equality test (Naive == Dgemm == Tile): {}
"=============================================
""".format(n, n, min_naive_time, min_dgemm_time, min_tile_time, equality))
        
