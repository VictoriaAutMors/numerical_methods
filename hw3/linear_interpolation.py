import numpy as np
import matplotlib.pyplot as plt

def reading(x, y, z):
    train_dat = open('train1.dat')
    for i in train_dat.readline().split():
        x.append(float(i))
    train_dat.close()

    train_ans = open('train1.ans')
    for i in train_ans.readline().split():	
        y.append(float(i))
    train_ans.close()

    test_dat = open('test1.dat')
    for i in test_dat.readline().split():
        z.append(float(i))
    test_dat.close()

    return x, y, z

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

x, y, z = [], [], []

x, y, z = reading(x, y, z)

n = int(len(x))
m = int(len(z))
print('x:', x)
print('y:', y)
print('z:', z)

a, b = cofficients(x, y)

test_ans = open('test1.ans', 'w')
answer = get_answer(test_ans, x, y, z, m)
for i in answer:
    test_ans.write(str(i) + " ")
test_ans.write("\n")

test_ans.close()

#Build My interpolate
print('Mine answer:')
for i in answer:
    print(i)

plt.plot(x, y) 
plt.plot(z, answer) 
plt.show()

#Build Linalg interpolate 
print('Linalg answer:')
for i in np.interp(z, x, y):
    print(i)

plt.plot(x, y)
plt.plot(z, np.interp(z, x, y))
plt.show()
