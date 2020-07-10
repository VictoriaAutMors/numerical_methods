# numerical_methods

## Content
* [Definitions and Basics](#basics)
* [Homework 1 - Exact Solution of Linear Systems](#hw1)
  * [Gaussian elimination](#gauss)
  * [Tridiagonal matrix or Thomas algorithm](#sweep)
  * [Cholesky decomposition](#cholesky)
* [Homework 2 - Iterative methods](#hw2)
  * [Seidel method](#seidel)
  * [Jacobi method](#jacobi)
* [Homework 3 - Interpolations](#hw3)
  * [Linear interpolation](#linear)
  * [Polynomial or Lagrange interpolation](#lagrange)
  * [Spline interpolation](#spline)
* [Homework 4 - Problems of mathematical physics](#hw4)
  * [Numerical methods for diffusion equations](#heat)
  * [Numerical methods for transfer equations](#transfer)
* [Dependencies](#dependencies)
  * [numpy](#numpy)
  * [scipy](#scipy)
  * [matplotlib](#plt)
  * [pygame](#pygame)
  * [ffmpeg](#ffmpeg)
* [How to run programs](#run)

# <a name="hw1"></a> Exact Solution of Linear Systems

## <a name="gauss"></a> Gaussian elimination
### Asymptotics: ![image1](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/Gauss1.svg)

Gaussian elimination, also known as row reduction, is an algorithm in linear algebra for solving a system of linear equations. It is usually understood as a sequence of operations performed on the corresponding matrix of coefficients. This method can also be used to find the rank of a matrix, to calculate the determinant of a matrix, and to calculate the inverse of an invertible square matrix. 

### code realisation in Python3: 
```
def forward(A, f, n):
    for k in range(n):
        A[k] = A[k] / A[k][k]
        f[k] = f[k] / A[k][k] 
        
        for i in range(k + 1, n):
            A[i] = A[i] - A[k] * A[i][k]
            f[i] = f[i] - f[k] * A[i][k]
            A[i][k] = 0
    return A, f
```
```
def backward(A, f, n):
    myAnswer = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = f[i]
        for j in range(i + 1, n):
            x[i] = x[i] - A[i][j] * x[j]
    return np.array(x)
```
### Plot of the data comparing linalg and my programm resolved time:

![image10](https://github.com/VictoriaAutMors/numerical_methods/blob/master/hw1/plots/Gauss.png)

vertical - time in seconds; horizontal - matrix size

## <a name="sweep"></a> Tridiagonal matrix or Thomas algorithm
### Asymptotics: ![image11](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/Tma3.svg)

In numerical linear algebra, the tridiagonal matrix algorithm, also known as the Thomas algorithm (named after Llewellyn Thomas), is a simplified form of Gaussian elimination that can be used to solve tridiagonal systems of equations. A tridiagonal system for n unknowns may be written as 

![image12](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/Tma1.svg)

![image13](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/Tma2.svg)

Thomas' algorithm is not stable in general, but is so in several special cases, such as when the matrix is diagonally dominant.

### code realisation in Python3: 
```
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
```
### Plot of the data comparing linalg and my programm resolved time:

![image14](https://github.com/VictoriaAutMors/numerical_methods/blob/master/hw1/plots/sweep.png)

## <a name="cholesky"></a> Cholesky decomposition

### Asymptotics: ![image15](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/Gauss1.svg)

The Cholesky decomposition of a Hermitian positive-definite matrix A is a decomposition of the form 
A = LL where L is a lower triangular matrix with real and positive diagonal entires, and L denotes the conjugate transpose of L. Every Hermitian positive-definite matrix (and thus also every real-valued symmetric positive-definite matrix) has a unique Cholesky decomposition.

### code realisation in Python3:
```
def Cholesky(A):
    n = len(A)
    L = [[0.0] * n for i in range(n)]
    for i in range(n):
        for k in range(i + 1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in range(k))
            if (i == k): 
                L[i][k] = math.sqrt(abs(A[i][i] - tmp_sum))
            else:
                L[i][k] = (1.0 / L[k][k] * (A[i][k] - tmp_sum))
    return L
```

### Plot of the data comparing linalg and my programm resolved time:

![image16](https://github.com/VictoriaAutMors/numerical_methods/blob/master/hw1/plots/Cholesky.png)

## <a name="hw3"></a> Homework 3 - Interpolations
In the mathematical field of numerical analysis, interpolation is a type of estimation, a method of constructing new data points within the range of a discrete set of known data points.

## <a name="linear"></a> Linear interpolation
One of the simplest methods is linear interpolation (sometimes known as lerp). Consider the above example of estimating f(2.5). Since 2.5 is midway between 2 and 3, it is reasonable to take f(2.5) midway between f(2) = 0.9093 and f(3) = 0.1411, which yields 0.5252. Linear interpolation is quick and easy, but it is not very precise. Generally, linear interpolation takes two data points, say (xa,ya) and (xb,yb), and the interpolant is given by: 

![image2](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/1LI.svg)

![image3](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/2LI.svg)

![image4](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/3LI.svg)

The error is proportional to the square of the distance between the data points. The error in some other methods, including polynomial interpolation and spline interpolation, is proportional to higher powers of the distance between the data points.The following error estimate shows that linear interpolation is not very precise.

### code realisation in Python3: 
```
def cofficients(x ,y):
    a = np.zeros(n)
    b = np.zeros(n)
    for i in range(n-1):
        a[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
        b[i] = y[i]
    return a, b


def get_answer(test_ans, x, y, z, m):
    answer = np.zeros(m)
    for j in range(0, m):
        for i in range(0 ,n - 1):
            if (z[j] < x[i + 1] and x[i] <= z[j]):
                answer[j] = a[i] * (z[j] - x[i]) + b[i]

    test_ans.write(str(answer[j]) + ' ')
    return answer
```
### Plot of the data with linear interpolation superimposed: 

![image5](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/5Li.png)

## <a name="lagrange"></a> Polynomial or Lagrange interpolation
Polynomial interpolation is a generalization of linear interpolation. Note that the linear interpolant is a linear function. We now replace this interpolant with a polynomial of higher degree. Consider again the problem given in linear interpolation. The following sixth degree polynomial goes through all the seven points: 
![image6](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/1La.svg)

Generally, if we have n data points, there is exactly one polynomial of degree at most n−1 going through all the data points. The interpolation error is proportional to the distance between the data points to the power n. Furthermore, the interpolant is a polynomial and thus infinitely differentiable. So, we see that polynomial interpolation overcomes most of the problems of linear interpolation. However, polynomial interpolation also has some disadvantages. Calculating the interpolating polynomial is computationally expensive compared to linear interpolation. Furthermore, polynomial interpolation may exhibit oscillatory artifacts, especially at the end points

### code realisation in Python3: 
```
def phi (i, x, y, t):
    val = 1
    for j in range(n):
        if j != i:
            val = val * ((t - x[j]) / ( x [ i ] - x [ j ]))
        if j == i:
            val *= 1
    return val

def Lagrange (x, y, t):
    ans = 0
    for i in range(n):
        ans = ans + y[i] * phi(i, x, y, t)
    return ans
```
### Plot of the data with polynomial interpolation applied:

![image7](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/2La.png)

## <a name="spline"></a> Spline interpolation

Remember that linear interpolation uses a linear function for each of intervals [xk,xk+1]. Spline interpolation uses low-degree polynomials in each of the intervals, and chooses the polynomial pieces such that they fit smoothly together. The resulting function is called a spline.

For instance, the natural cubic spline is piecewise cubic and twice continuously differentiable. Furthermore, its second derivative is zero at the end points. The natural cubic spline interpolating the points in the table above is given by

![image8](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/1Si.svg)

In this case we get f(2.5) = 0.5972.
Like polynomial interpolation, spline interpolation incurs a smaller error than linear interpolation, while the interpolant is smoother and easier to evaluate than the high-degree polynomials used in polynomial interpolation. 

### code realisation in Python3: 
```
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
```
### Plot of the data with spline interpolation applied:

![image9](https://github.com/VictoriaAutMors/numerical_methods/blob/master/Images/2Si.png)

# <a name="dependencies"></a> Dependencies
### <a name="numpy"></a> numpy

``` sudo apt-get install python3-numpy ```

``` pip3 install numpy ```

### <a name="scipy"></a> scipy

``` sudo apt-get install python3-scipy ```
``` pip3 install scipy ```

### <a name="plt"></a> matplotlib

``` sudo apt-get install python3-matplotlib ```
``` pip3 install matplotlib ```

### <a name="pygame"></a> pygame

``` sudo apt-get install python3-pygame ```

``` pip3 install pygame ```

### <a name="ffmpeg"></a> ffmpeg

``` sudo apt-get install python3-ffmpeg ```
``` pip3 install ffmpeg ```

# <a name="run"></a> How to run programs

``` python3 programName.py```

**Example:**   ``` python3 Lagrange_interpolate.py ```
