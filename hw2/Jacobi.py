import numpy as np
import time
import math

n = int(input())
def jacobi (A , f , x ):
    xnew = n * [0]
    for i in range(n):
        s = 0
        for j in range(0, i):
            s = s + A [ i ][ j ] * x [ j ]
        for j in range(i + 1, n):
            s = s + A [ i ][ j ] * x [ j ]
        xnew[ i ] = ( f [ i ] - s ) / A [ i ][ i ]
    return xnew

def diff(x, xnew):
    s = float(0)
    for i in range(n):
        s += float((x[i] - xnew[i]) ** 2) 
    return float(math.sqrt(s))

def solve (A , f):
    eps = 0.01
    xnew = np.zeros(n)
    while True:
        x = xnew
        xnew = jacobi (A , f , x )
        if ( diff(x, xnew ) < eps ):
            return x


A = np.random.rand(n, n)
f = np.random.rand(n)
for i in range(n):
    for j in range(n):
        if i != j:
            A[i][i] += A[i][j]
print("\nMy algorithm:\n")
my_time = time.time()
solve(A, f)
my_time = time.time() - my_time
for el in solve(A, f):
    print(el)
linalg_time = time.time()
print("\nLinalg:\n")
for el in np.linalg.solve(A, f):
    print(el)
linalg_time = time.time() - linalg_time

print("\nMy time:", my_time)
print("Linalg time:", linalg_time)
