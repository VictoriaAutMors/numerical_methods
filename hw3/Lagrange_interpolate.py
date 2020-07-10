import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange


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


def reading(x, y, z):
	train_dat = open('train2.dat')
	for i in train_dat.readline().split():
		x.append(float(i))
	train_dat.close()

	train_ans = open('train2.ans')
	for i in train_ans.readline().split():	
		y.append(float(i))
	train_ans.close()

	test_dat = open('test2.dat')
	for i in test_dat.readline().split():
		z.append(float(i))
	test_dat.close()

	return x, y, z

x, y, z = [], [], []
x, y, z = reading(x, y, z)

n = len(x)
m = len(y)

print('x:', x)
print('y:', y)

answer = []
test_ans = open('test2.ans', 'w')
for i in z:
	answer.append(Lagrange(x, y, i))
for i in answer:
	test_ans.write(str(i) + " ")
test_ans.write("\n")

test_ans.close()
print('z:', z)
print("My answer:")
for i in answer:
	print(i)


print("Linalg answer:")
for i in lagrange(x, y)(z):
	print(i)

#Build My interpolate
min_xz = min( np.min(x), np.min(z) )
max_xz =  max( np.max(x), np.max(z) )

xnew = np.linspace(min_xz , max_xz, 50 )
ynew = []
for el in xnew:
	ynew.append(Lagrange(x,y,el))

plt.plot(x, y, 'o', xnew, ynew)
plt.plot(z, answer, 'o')
plt.grid(True)
plt.show()


#Build Linalg interpolate 
min_xz = min( np.min(x), np.min(z) )
max_xz =  max( np.max(x), np.max(z) )

xnew = np.linspace(min_xz , max_xz, 50 )
ynew = []
for el in xnew:
	ynew.append(lagrange(x,y)(el))

plt.plot(x, y, 'o', xnew, ynew)
plt.plot(z, answer, 'o')
plt.grid(True)
plt.show()
