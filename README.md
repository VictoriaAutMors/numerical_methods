# numerical_methods

## Content
* [Definitions and Basics](#basics)
* [Exact Solution of Linear Systems](#hw1)
  * [Gaussian elimination](#gauss)
  * [Tridiagonal matrix algorithm](#sweep)
  * [Cholesky decomposition](#cholesky)
* [Iterative methods](#hw2)
  * [Seidel method](#seidel)
  * [Jacobi method](#jacobi)
* [Homework 3 - Interpolation](#hw3)
  * [Linear interpolation](#linear)
  * [Polynomial or Lagrange interpolation](#lagrange)
  * [Spline interpolation](#spline)
* [Problems of mathematical physics](#hw4)
  * [Numerical methods for diffusion equations](#heat)
  * [Numerical methods for transfer equations](#transfer)
* [Dependencies](#dependencies)
  * [numpy](#numpy)
  * [scipy](#scipy)
  * [matplotlib](#plt)
  * [pygame](#pygame)
  * [ffmpeg](#ffmpeg)
* [How to run programs](#run)
* [Questions ans suggestions](#questions)

## <a name="hw3"></a> Homework 3 - Interpolation
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
