#Gauss method
import numpy as np
import time

n = int(input())
A = np.random.rand(n, n)
f = np.random.rand(n)
A1 = A.copy()
f1 = f.copy()
for i in range(n):
    for j in range(n):
        if i != j:
            A[i][i] += A[i][j] 

linalg_time = time.time()
print("LINALG\n")
x1 = np.linalg.solve(A,f)
linalg_time = time.time() - linalg_time
for i in range(n):
    print(x1[i])

my_time = time.time()
for k in range(n):
        for j in range(k + 1, n):
            A[k, j] = A[k, j] / A[k, k]
        f[k] /= A[k][k]    
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = A[i][j] - A[i][k] * A[k][j]
            f[i] = f[i] - A[i][k] * f[k]
            A[i][k] = 0
x = np.zeros(n)
for i in range(n - 1, -1, -1):
    x[i] = f[i]
    for j in range(i + 1, n):
        x[i] -= A[i][j] * x[j]
my_time = time.time() - my_time
print("\nMy solution\n")
for i in range(n):
    print(x[i])
print("my time:", my_time)
print("linalg time", linalg_time)
