#! /usr/bin/env python3
import argparse
import time
import random
import itertools
from sys import stdout


import _matrix

parser = argparse.ArgumentParser(
    description='benchmark for matrix multiplcation')
parser.add_argument('--times',
    nargs = 1,
    default=[3],
    help="Number of repetition for a test case")

args = parser.parse_args()
times = int(args.times[0])

def time_delta(func):
    def time_wrap(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        delta = time.time() - start
        return result, delta
    return time_wrap


@time_delta
def MKL_MM(A, B):
    return _matrix.multiply_mkl(A, B)


@time_delta
def naive_MM(A, B):
    return _matrix.multiply_naive(A, B)

def copy_to_mat(mat, row, col, li):
  for i in range(row):
    for j in range(col):
      mat[i,j] = li[ i*mat.ncol+ j]

def log(s, file ,screen=stdout):
    if file != None:
        print( s, file=file)
    if screen != None:
        print(s, file=screen)

if __name__ == "__main__":
    f = open("performance.txt", 'w')

    test_ranges = [1000,1250]
    print("Begin benchmark...")
    test_cases = list(itertools.product(test_ranges,test_ranges,test_ranges))

    naive_runtime = [0]*len(test_cases)
    mkl_runtime = [0]*len(test_cases)
    mkl_speedup = [0]*len(test_cases)

    for idx,(m,n,k) in enumerate(test_cases):
        f.flush()
        log(f'test case {idx+1:2}/{len(test_cases):2}: , dim= {m:4} x {n:4} x {k:4} (repeat={times})',f)
        total_size = m * k

        population_size = 100
        random_population = [int(random.randrange(-500, 1000))
                            for _ in range(int(population_size))]

        naive_list = [0 for _ in range(times)]
        MKL_list = [0 for _ in range(times)]

        for i in range(times):
            A = _matrix.Matrix(m, n)
            A_content = random.choices(random_population, k=m * n)
            copy_to_mat(A,m, n, A_content )

            B = _matrix.Matrix(n, k)
            B_content = random.choices(random_population, k=n * k)
            copy_to_mat(B,n, k, B_content )

            naive_result, naive_time = naive_MM(A, B)
            MKL_result, MKL_time = MKL_MM(A, B)
            assert naive_result == MKL_result

            naive_list[i] = naive_time
            MKL_list[i] = MKL_time

        naive_runtime[idx] = sum(naive_list) / len(naive_list)
        mkl_runtime[idx] = sum(MKL_list) / len(MKL_list)
        mkl_speedup[idx] = naive_runtime[idx]/mkl_runtime[idx]

        log(f'    naive:{naive_runtime[idx]:10f}sec,  mkl:{mkl_runtime[idx]:10f}sec, mkl-speedup:{mkl_speedup[idx]:5.2f} x'
            , f)

    naive_time_total = sum(naive_runtime)
    mkl_time_total = sum(mkl_runtime)

    raw_speedup = 1
    for i in mkl_speedup:
        raw_speedup = raw_speedup * i
    avg_speedup = pow(raw_speedup, 1.0/len(mkl_speedup))

    log(f'MKL speed-up over naive (normalized): {avg_speedup} x', f)
    f.close()
    print("Finish benchmark...")
