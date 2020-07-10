import string
import numpy as np
from scipy import interpolate
import scipy
import matplotlib.pyplot as plt

def interpolate(t, m, n, x, A, B, C, D, handle):
	f = np.zeros(m)
	for j in range(m):
		for i in range(0, n - 1):
			if (x[i] <= t[j] and t[j] <= x[i+1]):
				f[j] = A[i]*(t[j]-x[i])**3 + B[i]*(t[j]-x[i])**2 + C[i]*(t[j]-x[i]) + D[i]
	return f


def sweep ( a, b, c, f, n):
    alpha = (n + 1) * [0]
    beta = (n + 1) * [0]
    x = n * [0]
    a[0] = 0
    c[n -  1] = 0
    alpha[0] = 0
    beta[0] = 0
    for i in range(0, n):  
        d = a[i] * alpha[i] + b[i]
        alpha [i + 1] = -c[i] / d
        
        beta [i + 1] = (f[i] - (a[i] * beta[i])) / (d)
    x[n - 1] = beta[n]
    for i in range(n - 2, -1, -1):
        x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1]
    return x

def generateSpline (x , y):
    n = x.shape[0] - 1
    h = ( x[ n ] - x[0]) / n
    a = np.array ([0] + [1] * ( n - 1) + [0])
    b = np.array ([1] + [4] * ( n - 1) + [1])
    c = np.array ([0] + [1] * ( n - 1) + [0])

    f = np.zeros (n + 1)
    for i in range (1 , n):
        f[ i ] = 3 * ( y [i -1] - 2 * y[i ] + y[i + 1]) / h**2
    s = sweep( a , b , c , f, n + 1)
    B = [0] * (n + 1)   
    A = [0] * (n + 1)
    C = [0] * (n + 1)
    D = [0] * (n + 1)
    for i in range (n):
        B [ i ] = s [ i ]
        D [ i ] = y [ i ]
    for i in range(n):
        A[i] = (B[i + 1] - B[i]) / (3 * h)
        C[i] = ((y[i + 1] - y[i]) / h) - ((B[i + 1] + 2 * B[i]) * h) / 3
    return A , B , C , D

handle = open("train.dat")

x = []
x = [float(i) for i in handle.readline().split()]
handle.close()
handle = open("train.ans")
y = []
y = [float(i) for i in handle.readline().split()]
handle.close()
handle = open('test.dat')
z = []
z = [float(i) for i in handle.readline().split()]
handle.close()
x = np.array(x)
print('X:', x)
y = np.array(y)
print('Y:', y)
z = np.array(z)
print('Z:', z)
n = int(len(x))
m = int(len(z))
A, B ,C , D = generateSpline(x, y)
handle = open('test.ans', 'w')
for el in A:
    handle.write(str(el) + " ")
handle.write("\n")
for el in B:
    handle.write(str(el) + " ")
handle.write("\n")
for el in C:
    handle.write(str(el) + " ")
handle.write("\n")
for el in D:
    handle.write(str(el) + " ")
handle.write("\n")

print("A:", A)
print("B:", B)
print("C:", C)
print("D:", D)  

function = interpolate(z, m,n,  x, A, B, C, D, handle)
print("My function:", function)
print("Linalg function:", scipy.interpolate.InterpolatedUnivariateSpline(x, y)(z))

min_xz = min( np.min(x), np.min(z) )
max_xz =  max( np.max(x), np.max(z) )

xnew = np.linspace(min_xz , max_xz, 50 )
ynew = interpolate(xnew, len(xnew), n, x, A, B, C, D, handle)


plt.plot(x, y, 'o', xnew, ynew)
plt.plot(z, function, 'o')
plt.grid(True)
plt.show()

handle.close()
