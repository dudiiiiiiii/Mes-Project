import numpy as np
from numpy.linalg import solve
from matplotlib.pyplot import plot, show
from math import pi
from scipy.integrate import quad
from scipy.special.orthogonal import p_roots

n = int(input("n : "))
#n = 7
h = 3 / (n-1)
a = 0
b = 3

G = 6.67 * pow(10, -11)

def en(n):
    def ex(x):
        if h*n >= x >= h*(n-1):
            return (x/h) -n +1
        elif x > h*n and x <= h*(n+1):
            return n - (x/h) + 1
        return 0
    return ex


B = [[0 for _ in range(n - 2)] for _ in range(n - 2)]

#https://austingwalters.com/gaussian-quadrature/
def gauss(f, a, b, n=10):
    [x, w] = p_roots(n + 1)
    return 0.5 * (b - a) * sum(w * f(0.5 * (b - a) * x + 0.5 * (b + a)))

def fillB(w, v):

    # pochodna fi razy pochodna v zawsze będzie równa +/- 1/h^2
    # (całkuję w granicach danej funkcji testowej e)

    if w == v:
        f = lambda x: -1/pow(h,2)
        res = gauss(f, h*(w-1), h*(w+1))
    else:
        f = lambda x: 1/pow(h, 2)
        p = min(v,w)
        res = gauss(f, h*p, h*(p+1))

    return res

for w in range(0, n-2):
    for v in range(max(w-1, 0), min(n-2, w+2)):
        B[w][v] = fillB(w, v)

L = [0 for i in range(n-2)]

def fillL(v):
    #return 4 * pi * G * gauss(en(v), 1, 2) error
    return 4 * pi * G * quad(en(v), 1, 2)[0]

for v in range(n-2):
    L[v] = fillL(v)

X = solve(B, L)

u = lambda x: - x/ 3 + 5

# def w(x):
#     res = 0
#     for q in range(n - 2):
#         res += (X[q] * en(q + 1)(x))
#         return  res

def w(x):
    return sum(X[q] * en(q + 1)(x) for q in range(n - 2))

# res function
fi = lambda x: u(x) + w(x)

p = 150
xlist = [3 * i / (p - 1) for i in range(p)]
ylist = [fi(i) for i in xlist]
# print(xlist)
# print(ylist)
plot(xlist, ylist)
show()


