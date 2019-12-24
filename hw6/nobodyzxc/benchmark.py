import time
import random
import _matrix

times = 5
naive_base = _matrix.Matrix(1024, 1024)
tile_base = _matrix.Matrix(1024, 1024)
cblas_base = _matrix.Matrix(1024, 1024)
multiplier = _matrix.Matrix(1024, 1024)

for i in range(1024):
    for j in range(1024):
        base, muler = random.random(), random.random()
        naive_base[i, j] = base
        tile_base[i, j] = base
        cblas_base[i, j] = base
        multiplier[i, j] = muler

print("run naive multiplication...")
naive_times = []
for i in range(times):
    beg = time.time()
    naive_base = _matrix.multiply_naive(naive_base, multiplier)
    end = time.time()
    naive_times.append(end - beg)

min_naive = min(naive_times)

print("run tiled multiplication...")
tile_times = []
for i in range(times):
    beg = time.time()
    tile_base = _matrix.multiply_tile(tile_base, multiplier, 64)
    end = time.time()
    tile_times.append(end - beg)

min_tile = min(tile_times)


print("run cblas multiplication...")
cblas_times = []
for i in range(times):
    beg = time.time()
    cblas_base = _matrix.multiply_mkl(cblas_base, multiplier)
    end = time.time()
    cblas_times.append(end - beg)

min_cblas = min(cblas_times)

record = open("performance.txt", 'w')
record.write(
"""Start multiply_naive (repeat={}), take min = {} seconds
Start multiply_tile (repeat={}), take min = {} seconds
Start multiply_mkl (repeat={}), take min = {} seconds
MKL speed-up over naive: {:.4f} x
tile speed-up over naive: {:.4f} x""".format(
    times, min_naive,
    times, min_tile,
    times, min_cblas,
    min_naive / min_cblas,
    min_naive / min_tile))
