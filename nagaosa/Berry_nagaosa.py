import numpy as np
from numpy import pi
from numpy import *
import scipy as sp
import matplotlib.pyplot as plt
import sympy as sy
from sympy.plotting import plot
from sympy.plotting import plot3d
from mpl_toolkits.mplot3d import Axes3D
import time
from pylab import *

start=time.time()

kb = 8.6*10**(-5)
e = 1.6*10**(-19)

es = ep = 1
ts = tp = tsp = 1/np.sqrt(2)

a = 1
g1 = 2*pi*np.array([1, 0])/a
g2 = 2*pi*np.array([0, 1])/a

s1 = np.array([ [0, -1j],
                [1j, 0]])
s2 = np.array([ [0, 1],
                [1, 0]])
def h(kx, ky):
    return np.array([   [es-2*ts*(np.cos(kx)+np.cos(ky)), np.sqrt(2)*tsp*(1j*np.sin(kx)+np.sin(ky))],
                        [np.sqrt(2)*tsp*(-1j*np.sin(kx)+np.sin(ky)), ep+tp*(np.cos(kx)+np.cos(ky))]], dtype = complex)
def eigvals(kx, ky):
    return sp.linalg.eigh(h(kx, ky))[0]
def eigvecs(kx, ky):
    return sp.linalg.eigh(h(kx, ky))[1]

# Band Structure
# N = 64
# fig, ax = plt.subplots()
# q_ary = linspace(0, 1, N)
# for l in range(2):
#     val = zeros(N)
#     offset = zeros(4)
#     # Γ-X
#     i = 0
#     for q in q_ary:
#         K = 0*g1+q*g2/2
#         val[i] = eigvals(K[0], K[1])[l]
#         i = i+1
#     offset[0] = 0
#     plot(sum(offset)+q_ary*norm(g2/2), val, '#377eb8')
#     # X-M
#     i = 0
#     for q in q_ary:
#         K = q*g1/2+g2/2
#         val[i] = eigvals(K[0], K[1])[l]
#         i = i+1
#     offset[1] = norm(g2/2)
#     plot(sum(offset)+q_ary*(norm(g1/2)), val, '#377eb8')
#     # M-Γ
#     i = 0
#     for q in q_ary:
#         K = (1/2-q/2)*g1+(1/2-q/2)*g2
#         val[i] = eigvals(K[0], K[1])[l]
#         i = i+1
#     offset[2] = norm(g1/2)
#     plot(sum(offset)+q_ary*norm(g1/2+g2/2), val, '#377eb8')
#     offset[3] = norm(g1/2+g2/2)
# ax.set_xlim([0, sum(offset)])
# ax.set_ylim([-2.5, 4])
# plt.xticks(
#     [sum(offset[0:1]), sum(offset[0:2]), sum(offset[0:3]), sum(offset)],
#     ["Γ", "X", "M", "Γ"])
# for i in range(4):
#     plt.vlines(sum(offset[0:i+1]), -2.5, 4, "k", linewidth = 0.25)
# plt.hlines(0, 0, sum(offset), "m", linestyle = ":")
# ax.set_title("Band Structure of Nagaosa Model")
# print(float(time.time())-start)
# plt.show()

# DOS
# eps = 0.01
# def d(x):
#     return eps/(pi*(x**2+eps**2))
# def diffdos(kx, ky, E):
#     return d(E-eigvals(kx, ky)[0])+d(E-eigvals(kx, ky)[1])
# N = 256
# vec = np.zeros((3, N**2), complex)
# def midpoint(E):
#     g1h, g2h, ans = g1/N, g2/N, 0
#     for i in range(N):
#         gi = (2*i+1)/2*g1h
#         for j in range(N):
#             gj = (2*j+1)/2*g2h
#             kx = gi[0]+gj[0]
#             ky = gi[1]+gj[1]
#             ans = ans+(1/2)*(2*pi/(3*N))*(2*pi/(3*N))*diffdos(kx, ky, E)
#     return ans
# E = linspace(-2, 4, N)
# fig, ax = plt.subplots()
# ax.set_xlim([-2, 4])
# ax.set_ylim([0, 1.5])
# ax.set_title("DOS of Nagaosa Model")

# plot(E, 6*midpoint(E)/(2*pi**2), '#377eb8')

# plt.xlabel('E')
# plt.ylabel('D(E)')
# plt.legend()
# print(float(time.time())-start)
# plt.show()


def hx(kx, ky):
    return np.array([   [2*ts*np.sin(kx), np.sqrt(2)*1j*tsp*np.cos(kx)],
                        [-np.sqrt(2)*1j*tsp*np.cos(kx), -tp*np.sin(kx)] ])
def hy(kx, ky):
    return np.array([   [2*ts*np.sin(ky), np.sqrt(2)*tsp*np.cos(ky)],
                        [np.sqrt(2)*tsp*np.cos(ky), -tp*np.sin(ky)] ])
def berry(n, kx, ky):
    delta = 10**(-12)
    B = 0
    for m in (0, 1):
        B = B + np.conjugate(np.transpose(eigvecs(kx, ky)[n]))@hx(kx, ky)@eigvecs(kx, ky)[m]*np.conjugate(np.transpose(eigvecs(kx, ky)[m]))@hy(kx, ky)@eigvecs(kx, ky)[n]/((eigvals(kx, ky)[n]-eigvals(kx, ky)[m])**2+delta)
    return -2*np.imag(B)

#Berry Curvature
N = 128
fig, ax = plt.subplots()
q_ary = linspace(0, 1, N)
for n in range(2):
    val = zeros(N)
    offset = zeros(4)
    # Γ-X
    i = 0
    for q in q_ary:
        K = 0*g1+q*g2/2
        val[i] = berry(n, K[0], K[1])
        i = i+1
    offset[0] = 0
    plot(sum(offset)+q_ary*norm(g2/2), val, '#377eb8')
    # X-M
    i = 0
    for q in q_ary:
        K = q*g1/2+g2/2
        val[i] = berry(n, K[0], K[1])
        i = i+1
    offset[1] = norm(g2/2)
    plot(sum(offset)+q_ary*(norm(g1/2)), val, '#377eb8')
    # M-Γ
    i = 0
    for q in q_ary:
        K = (1/2-q/2)*g1+(1/2-q/2)*g2
        val[i] = berry(n, K[0], K[1])
        i = i+1
    offset[2] = norm(g1/2)
    plot(sum(offset)+q_ary*norm(g1/2+g2/2), val, '#377eb8')
    offset[3] = norm(g1/2+g2/2)
ax.set_xlim([0, sum(offset)])
ax.set_ylim([-6, 6])
plt.xticks(
    [sum(offset[0:1]), sum(offset[0:2]), sum(offset[0:3]), sum(offset)],
    ["Γ", "X", "M", "Γ"])
for i in range(4):
    plt.vlines(sum(offset[0:i+1]), -6, 6, "k", linewidth = 0.25)
plt.hlines(0, 0, sum(offset), "m", linestyle = ":")
ax.set_ylabel('Ω_z')
ax.set_title("z component of Berry Curvature of Nagaosa Model")


N = 64
# def density(n, mu, T):
#     beta = 1/(kb*T)
#     g1h, g2h, ans = g1/N, g2/N, 0
#     for i in range(N):
#         gi = (2*i+1)/2*g1h
#         for j in range(N):
#             gj = (2*j+1)/2*g2h
#             kx = gi[0]+gj[0]
#             ky = gi[1]+gj[1]
#             ans = ans+((2*pi/(a*N))**2)*(1/(exp(beta*(eigvals(kx, ky)[n]-mu))+1))
#     return ans/(2*pi**2)
def ahe(n, mu, T):
    beta = 1/(kb*T)
    g1h, g2h, ans = g1/N, g2/N, 0
    for i in range(N):
        gi = (2*i+1)/2*g1h
        for j in range(N):
            gj = (2*j+1)/2*g2h
            kx = gi[0]+gj[0]
            ky = gi[1]+gj[1]
            ans = ans+((2*pi/(a*N))**2)*berry(n, kx, ky)*(1/(exp(beta*(eigvals(kx, ky)[n]-mu))+1))
    return e**2*ans/(2*pi**2)

#fig, ax = plt.subplots()
mu_array = linspace(-2.5, 4.5, N)
# dens = zeros((2, N))
# for n in range(2):
#     i = 0
#     for mu in mu_array:
#         dens[n][i] = density(n, mu, 300)
#         i = i+1
#     plot(mu_array, dens[n], label='n_%i'%n)


# ahes = zeros((2, N))
# for n in range(2):
#     i = 0
#     for mu in mu_array:
#         ahes[n][i] = ahe(n, mu, 300)
#         i = i+1
#     plot(mu_array, ahes[n], label='σ_%i'%n)

# ax.set_xlabel('μ')
# #ax.set_ylabel('n(μ)')
# ax.set_ylabel('σ_{xy}')
# plot(mu_array, ahes[0]+ahes[1], label='σ_0+σ_1')
# plt.legend()
print(float(time.time())-start)
plt.show()