#sweep method
import numpy as np
from scipy.linalg import solve_banded
import time

def sweep ( a, b, c, f, n):
    alpha = (n + 1) * [0]
    beta = (n + 1) * [0]
    x = np.random.rand(n)
    a[0] = 0
    c[n -  1] = 0
    alpha[0] = 0
    beta[0] = 0
    for i in range(0, n):  
        d = float(a[i] * alpha[i] + b[i])
        alpha [i + 1] = float(-c[i] / d)
        
        beta [i + 1] = float((f[i] - (a[i] * beta[i])) / (d))
    x[n - 1] = float(beta[n])
    for i in range(n - 2, -1, -1):
        x[i] = float(alpha[i + 1] * x[i + 1] + beta[i + 1])
    return x

n = int(input()) 

a = np.random.rand(n)
a[0] = 0
b = np.random.rand(n)
c = np.random.rand(n)
c[n - 1] = 0
for i in range(n):
    b[i] += a[i] + b[i]

f = np.random.rand(n)
a_linalg = np.array([0] + c[:-1].tolist())
b_linalg = b
c_linalg = np.array(a[1:].tolist() + [0])
A = np.array([a_linalg, b_linalg, c_linalg])
my_time = time.time()
print("\nMy answer:", sweep(a, b, c, f, n))
my_time = time.time() - my_time
linalg_time = time.time()
print("\nLinalg answer:", solve_banded((1, 1), A, f))
linalg_time = time.time() - linalg_time
print("\nMy time:", my_time)
print("Linalg time:", linalg_time)
