import time
import numpy as np
import math

n = int(input())
def seidel (A , f , x ):
    xnew = np.zeros(n)
    for i in range(n):
        s = float(0)
        for j in range(i):
            s = float(s + A [ i ][ j ] * xnew [ j ])
        for j in range(i + 1, n):
            s = float(s + A[ i ][ j ] * x[ j ])
        xnew[ i ] = float(( f[ i ] - s) / A [ i ][ i ])
    return xnew

def diff(x, xnew):
    s = float(0)
    for i in range(0, n):
        s += float((x[i] - xnew[i]) ** 2)
    return float(math.sqrt(s))

def solve (A , f ):
    xnew = np.zeros(n)
    eps = 0.01
    while True:
        x = np.array(xnew)
        xnew = seidel (A , f , x )
        if(diff(x, xnew) < eps):
            break
    return x

A = np.random.rand(n, n)
f = np.random.rand(n)
for i in range(n):
    for j in range(n):
        if i != j:
            A[i][i] += A[i][j]

print("\nMy algorithm:\n")
my_time = time.time()
for i in solve(A, f):
    print(i)
my_time = time.time() - my_time
print("\nLinalg:\n")
linalg_time = time.time()
for i in np.linalg.solve(A, f):
    print(i)
linalg_time = time.time() - linalg_time

print("\nMy time:", my_time)
print("Linalg time:", linalg_time)
