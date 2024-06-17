import numpy as np
from numpy import pi
import scipy as sp
import sympy as sy
from sympy.plotting import plot
from sympy.plotting import plot3d
from scipy.integrate import dblquad
from scipy import integrate
import matplotlib.pyplot as plt
import time
from pylab import *

start=time.time()

#Analytical
'''
N = 100
x1 = np.linspace(-np.pi, np.pi, N)
x2 = np.linspace(-np.pi, np.pi, N)

X1, X2 = np.meshgrid(x1, x2)
X = np.c_[np.ravel(X1), np.ravel(X2)]

def Y1(X1, X2):
    return -1+np.sqrt(2*(np.cos(X1)+np.cos(X2)+np.cos(X1+X2))+3)
def Y2(X1, X2):
    return -1-np.sqrt(2*(np.cos(X1)+np.cos(X2)+np.cos(X1+X2))+3)
def Y3(X1, X2):
    return 2+X1*0

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X1, X2, Y1(X1, X2), cmap='bwr', linewidth=0)
surf = ax.plot_surface(X1, X2, Y2(X1, X2), cmap='bwr', linewidth=0)
surf = ax.plot_surface(X1, X2, Y3(X1, X2), cmap='bwr', linewidth=0)
ax.set_title("Band Structure of Kagomé Lattice")

plt.legend()
ax.set_xlabel('θ1')
ax.set_ylabel('θ2')
ax.set_zlabel('E')
plt.show()

'''

eps = 0.01
t = 1
def d(x):
    return eps/(pi*(x**2+eps**2))

#def E4(x, y):
    #return -2*t*(np.cos(x)+np.cos(y))
#plot3d(E1(k1, k2), (k1, -3, 3), (k2, -3, 3))
#enumerate



fig, ax = plt.subplots()
for a in (1,):
    def q1(x, y):
        return a*(x+sqrt(3)*y)
    def q2(x, y):
        return a*(x-sqrt(3)*y)
    def E1(x, y):
        return -t+t*sqrt(2*(cos(q1(x, y))+cos(q2(x, y))+cos(q1(x, y)+q2(x, y)))+3)
    def E2(x, y):
        return -t-t*sqrt(2*(cos(q1(x, y))+cos(q2(x, y))+cos(q1(x, y)+q2(x, y)))+3)
    def E3(x, y):
        return 2*t
    def diffdos(x, y, E):
        return (d(E-E1(q1(x, y), q2(x, y)))+d(E-E2(q1(x, y), q2(x, y))))/(2*pi**2)
    
    n = 128
    E_ary = linspace(-4.5, 2.5, n)
    val1 = val2 = val3 = zeros(n)   
    idx = 0
    diffdos2 = lambda y, x: diffdos(x, y, E)
    for E in E_ary:
        val1[idx], err = integrate.dblquad(diffdos2, -2*pi/(3*a), -pi/(3*a), lambda x: -sqrt(3)*(x+(2*pi/(3*a))), lambda x: sqrt(3)*(x+(2*pi/(3*a))), epsabs = 0.05)
        val2[idx], err = integrate.dblquad(diffdos2, -pi/(3*a), pi/(3*a), -pi/(sqrt(3)*a), pi/(sqrt(3)*a), epsabs = 0.05)
        val3[idx], err = integrate.dblquad(diffdos2, pi/(3*a), 2*pi/(3*a), lambda x: sqrt(3)*(x-(2*pi/(3*a))), lambda x: -sqrt(3)*(x-(2*pi/(3*a))), epsabs = 0.05)
        idx = idx + 1
    ax.plot(E_ary, val1+val2+val3, label='a=%a' %a)
ax.set_xlim([-4.5, 2.5])
ax.set_ylim([0, 0.5])
plt.vlines(2, 0, 1.5, '#377eb8')
plt.xlabel('E')
plt.ylabel('D(E)')
#plt.legend()
ax.set_title("DoS of Kagomé Lattice")
print(float(time.time())-start)
plt.show()

'''

t = 1
a = 1
def q1(x, y):
    return a*(x+sqrt(3)*y)
def q2(x, y):
    return a*(x-sqrt(3)*y)
def E1(x, y, a):
    return -t+t*sqrt(2*(cos(q1(x, y))+cos(q2(x, y))+cos(q1(x, y)+q2(x, y)))+3)
def E2(x, y, a):
    return -t-t*sqrt(2*(cos(q1(x, y))+cos(q2(x, y))+cos(q1(x, y)+q2(x, y)))+3)
def E3(x, y, a):
    return 2*t

fig, ax = plt.subplots()

ax.set_xlim([-0.1, 4*pi/3+0.1])
ax.set_ylim([-4.1, 2.1])
plt.hlines(2, 0, 4*pi/3, "#4daf4a")
x = linspace(0, 2*pi/(3*a), 100)
plot(x, E1(x, 0, a), "#377eb8")
plot(x, E2(x, 0, a), "#ff7f00")
x = linspace(2*pi/(3*a), pi/(2*a), 100)
plot(4*pi/(3*a)-x, E1(x, -sqrt(3)*(x-(2*pi/(3*a))), a), "#377eb8")
plot(4*pi/(3*a)-x, E2(x, -sqrt(3)*(x-(2*pi/(3*a))), a), "#ff7f00")
x = linspace(pi/(2*a), 0, 100)
plot(4*pi/(3*a)-x, E1(x, (1/sqrt(3))*x, a), "#377eb8")
plot(4*pi/(3*a)-x, E2(x, (1/sqrt(3))*x, a), "#ff7f00")

plt.vlines(0, -4.1, 4, "k", linewidth = 0.25)
plt.vlines(2*pi/3, -4.1, 4, "k", linewidth = 0.25)
plt.vlines(5*pi/6, -4.1, 4, "k", linewidth = 0.25)
plt.vlines(4*pi/3, -4.1, 4, "k", linewidth = 0.25)
plt.hlines(-1, -0.1, 4*pi/3+0.1, "m", linestyle = ":")
plt.text(0, -4.1, "Γ", fontsize=14)
plt.text(2*pi/3, -4.1, "K", fontsize=14)
plt.text(5*pi/6, -4.1, "M", fontsize=14)
plt.text(4*pi/3, -4.1, "Γ", fontsize=14)
ax.tick_params(labelbottom = False)
ax.set_title("Band Structure of Kagomé Lattice")
plt.xlabel('k')
plt.ylabel('E')
plt.show()


#Numerical
a = 1
def q1(x, y):
    return a*(x+sqrt(3)*y)
def q2(x, y):
    return a*(x-sqrt(3)*y)
def H(q1, q2):
    return -np.array([  [0, 1+exp(1j*q1) ,1+exp(-1j*q2)],
                        [1+exp(-1j*q1), 0, 1+exp(-1j*(q1+q2))],
                        [1+exp(1j*q2), 1+exp(1j*(q1+q2)), 0]])

N = 64
x_ary = y_ary = linspace(-4, 4, N)
eigvals = np.array([[np.sort(sp.linalg.eigh(H(q1(x, y), q2(x, y)), eigvals_only = True)) for x in x_ary] for y in y_ary] )

def eigvals(x, y):
    return np.sort(sp.linalg.eigh(H(q1(x, y), q2(x, y)), eigvals_only = True))



n = 128
fig, ax = plt.subplots()
ax.set_xlim([-0.1, 4*pi/3+0.1])
ax.set_ylim([-4.1, 2.1])

val0 = val1 = val2 = zeros(n)

xx1 = linspace (0, 2*pi/(3*a), n)
for i in range(n):
    x = i*2*pi/(3*a*n)
    val0[i] = eigvals(x, 0)[0]
plot(xx1, val0, "#ff7f00")
for i in range(n):
    x = i*2*pi/(3*a*n)
    val1[i] = eigvals(x, 0)[1]
plot(xx1, val1, "#377eb8")
for i in range(n):
    x = i*2*pi/(3*a*n)
    val2[i] = eigvals(x, 0)[2]
plot(xx1, val2, "#4daf4a", linewidth = 3)
print(val0, val1 ,val2)


xx2 = linspace(2*pi/(3*a), pi/(2*a), n)
def eigvals2(x):
    return eigvals(x, -sqrt(3)*(x-(2*pi/(3*a))))
i = 0
for x in xx2:
    val0[i] = eigvals2(x)[0]
    i = i+1
plot(4*pi/(3*a)-xx2, val0, "#ff7f00")
i = 0
for x in xx2:
    val1[i] = eigvals2(x)[1]
    i = i+1
plot(4*pi/(3*a)-xx2, val1, "#377eb8")
i = 0
for x in xx2:
    val2[i] = eigvals2(x)[2]
    i = i+1
plot(4*pi/(3*a)-xx2, val2, "#4daf4a")

xx3 = linspace(pi/(2*a), 0, n)
def eigvals3(x):
    return eigvals(x, (1/sqrt(3))*x)
i = 0
for x in xx3:
    val0[i] = eigvals3(x)[0]
    i = i+1
plot(4*pi/(3*a)-xx3, val0, "#ff7f00")
i = 0
for x in xx3:
    val1[i] = eigvals3(x)[1]
    i = i+1
plot(4*pi/(3*a)-xx3, val1, "#377eb8")
i = 0
for x in xx3:
    val2[i] = eigvals3(x)[2]
    i = i+1
plot(4*pi/(3*a)-xx3, val2, "#4daf4a")

plt.vlines(0, -4.1, 4, "k", linewidth = 0.25)
plt.vlines(2*pi/3, -4.1, 4, "k", linewidth = 0.25)
plt.vlines(5*pi/6, -4.1, 4, "k", linewidth = 0.25)
plt.vlines(4*pi/3, -4.1, 4, "k", linewidth = 0.25)
plt.hlines(-1, -0.1, 4*pi/3+0.1, "m", linestyle = ":")
plt.text(0, -4.1, "Γ", fontsize=14)
plt.text(2*pi/3, -4.1, "K", fontsize=14)
plt.text(5*pi/6, -4.1, "M", fontsize=14)
plt.text(4*pi/3, -4.1, "Γ", fontsize=14)
ax.tick_params(labelbottom=False)
ax.set_title("Band Structure of Kagomé Lattice(Numerical)")
plt.xlabel('k')
plt.ylabel('E')

plt.show()


eps = 0.01
def d(x):
    return eps/(pi*(x**2+eps**2))
n = 256


def diffdos(x, y, E):
    return d(E-eigvals(x, y)[0])+d(E-eigvals(x, y)[1])+d(E-eigvals(x, y)[2])
def diffdos1(x, y, E):
    return d(E-eigvals(x, y)[0])
def diffdos2(x, y, E):
    return d(E-eigvals(x, y)[1])
def diffdos3(x, y, E):
    return d(E-eigvals(x, y)[2])
diffdoss = (diffdos, diffdos1, diffdos2, diffdos3)
def midpoint(E, j):
    xh, ans = 4*pi/(3*a*n), 0
    for i in range(n):
        xi = -2*pi/(3*a)+xh*(2*i+1)/2
        if -2*pi/(3*a) <= xi <= -pi/(3*a):
            for i in range(n):
                yh = 2*sqrt(3)*(xi+(2*pi/(3*a)))/n
                yi = -sqrt(3)*(xi+(2*pi/(3*a)))+i*yh
                ans = ans+xh*yh*diffdoss[j](xi, yi, E)
        if -pi/(3*a) <= xi <= pi/(3*a):
            for i in range(n):
                yh = 2*pi/(sqrt(3)*a*n)
                yi = -pi/(sqrt(3)*a)+i*yh
                ans = ans+xh*yh*diffdoss[j](xi, yi, E)
        if pi/(3*a) <= xi <= 2*pi/(3*a):
            for i in range(n):
                yh = -2*sqrt(3)*(xi-(2*pi/(3*a)))/n
                yi = sqrt(3)*(xi-(2*pi/(3*a)))+i*yh
                ans = ans+xh*yh*diffdoss[j](xi, yi, E)
    return ans

E = linspace(-4.5, 2.5, 128)
fig, ax = plt.subplots()
ax.set_xlim([-4.5, 2.5])
ax.set_ylim([0, 1.5])
plot(E, midpoint(E, 0)/(2*pi**2))

#plot(E, midpoint(E, 1)/(2*pi**2), label = 'E1')
#plot(E, midpoint(E, 2)/(2*pi**2), label = 'E2')
#plot(E, midpoint(E, 3)/(2*pi**2), label = 'E3')

plt.xlabel('E')
plt.ylabel('D(E)')
plt.legend()
ax.set_title("DoS of Kagomé Lattice(Numerical)")

plt.show()
'''
