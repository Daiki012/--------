import numpy as np
from numpy import pi
import scipy as sp
import sympy as sy
from sympy import *
I = sy.I
from sympy.plotting import plot
from sympy.plotting import plot3d
from scipy.integrate import dblquad
from scipy import integrate
import matplotlib.pyplot as plt
import time
from pylab import *

start=time.time()

a = 1
def q1(x, y):
    return a*(x+sqrt(3)*y)
def q2(x, y):
    return a*(x-sqrt(3)*y)
t = 1
def H(q1, q2):
    return -t*np.array([    [0, 1+exp(1j*q1) ,1+exp(-1j*q2)],
                            [1+exp(-1j*q1), 0, 1+exp(-1j*(q1+q2))],
                            [1+exp(1j*q2), 1+exp(1j*(q1+q2)), 0]])
def eigvals(x, y):
    return sp.linalg.eigh(H(q1(x, y), q2(x, y)))[0]
N = 2
x_ary = y_ary = linspace(0, 0.5, N)

g1 = 2*pi*np.array([sqrt(3)/2, 1/2])/(a*sqrt(3))
g2 = 2*pi*np.array([sqrt(3)/2, -1/2])/(a*sqrt(3))
g3 = g1+g2

# vec = np.zeros((3, N**2), complex)
# weight = zeros(N**2)
# for i in (0, 1, 2):
#     j = 0
#     for x in x_ary:
#         for y in y_ary:
#             k = g1*x+g2*y
#             # vals, vecs = sp.linalg.eigh(H(q1(x, y), q2(x, y)))
#             # vec = sp.linalg.eigh(H(q1(x, y), q2(x, y)))[2][All,3]
#             vec[i, j] = sp.linalg.eigh(H(q1(x, y), q2(x, y)))[1][2][i]
#             j = j+1
# def weight(i, j):
#     return np.abs(vec[i, j])**2
#print(weight(0,0))

#pDOS

eps = 0.01
def d(x):
    return eps/(pi*(x**2+eps**2))
def diffdos(x, y, E):
    return d(E-eigvals(x, y)[0])+d(E-eigvals(x, y)[1])+d(E-eigvals(x, y)[2])
def diffdosvec(x, y, E):
    return array([d(E-eigvals(x, y)[0]), d(E-eigvals(x, y)[1]), d(E-eigvals(x, y)[2])])
N = 256
vec = np.zeros((3, N**2), complex)
def midpoint(E, k):
    g1h, g3h, ans = g1/N, g3/N, 0
    for i in range(N):
        gi = (2*i+1)/2*g1h
        for j in range(N):
            gj = (2*j+1)/2*g3h
            x = gi[0]+gj[0]
            y = gi[1]+gj[1]
            vec = sp.linalg.eigh(H(q1(x, y), q2(x, y)))[1][k]
            weight = np.abs(vec)**2
            #print(weight)
            ans = ans+(1/2)*(2*pi/(3*N))*(2*pi/(3*N))*(weight.dot(diffdosvec(x, y, E)))
    return ans
E = linspace(-4.5, 2.5, N)
fig, ax = plt.subplots()
ax.set_xlim([-4.5, 2.5])
ax.set_ylim([0, 1.5])
ax.set_title("PDOS of Kagomé Lattice(Numerical)")

plot(E, 6*midpoint(E, 0)/(2*pi**2), '#377eb8', label = 'E1')
plot(E, 6*midpoint(E, 1)/(2*pi**2), '#ff7f00', label = 'E2')
plot(E, 6*midpoint(E, 2)/(2*pi**2), '#4daf4a', label = 'E3')

plt.xlabel('E')
plt.ylabel('D(E)')
plt.legend()
print(float(time.time())-start)
plt.show()


# eigvals, eigvecs = np.array([[np.sort(sp.linalg.eigh(H(q1(x, y), q2(x, y)))) for x in x_ary] for y in y_ary])
# sample = [[sp.linalg.eigh(H(q1(x, y), q2(x, y))) for x in x_ary] for y in y_ary]

# print(sample)
# print(len(sample))

# def eigvals(x, y):
#     return np.sort(sp.linalg.eigh(H(q1(x, y), q2(x, y))))



# Berry phase
# q1 = sy.Symbol('q1')
# q2 = sy.Symbol('q2')
# H2 = sy.Matrix([        [0, 1+sy.exp(I*q1) ,1+sy.exp(-I*q2)],
#                         [1+sy.exp(-I*q1), 0, 1+sy.exp(-I*(q1+q2))],
#                         [1+sy.exp(I*q2), 1+sy.exp(I*(q1+q2)), 0]])
# vec = H2.eigenvects()[0][2][0]
# print(vec[0])
# diffq1 = sy.Matrix([[sy.diff(vec[0], q1)], [sy.diff(vec[1], q1)], [sy.diff(vec[2], q1)]])
# diffq2 = sy.Matrix([[sy.diff(vec[0], q2)], [sy.diff(vec[1], q2)], [sy.diff(vec[2], q2)]])
# print(diffq1.dot(vec)-diffq2.dot(vec))
# plot3d(diffq1.dot(vec)-diffq2.dot(vec), (q1, -1, 1), (q2, -2*pi, 2*pi))