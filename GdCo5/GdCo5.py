import numpy as np
from numpy import pi
import scipy as sp
from sympy.plotting import plot
from scipy.integrate import dblquad
from scipy import integrate
import matplotlib.pyplot as plt
import pandas as pd
import time
from pylab import *

start=time.time()

a = b = 4.91543215799018
c = 3.94232418211356
r1 = a*np.array([1, 0, 0])
r2 = b*np.array([-1, sqrt(3), 0])/2
r3 = c*np.array([0, 0, 1])
g1 = 2*pi*np.array([1, 1/sqrt(3), 0])/a
g2 = 2*pi*np.array([0, 2/sqrt(3), 0])/b
g3 = 2*pi*np.array([0, 0, 1])/c

df = pd.read_table('/Users/ikedadaiki/Python3/物理工学特別輪講/GdCo5/tb_fe_10_nozeros.dat', header = None)
# print(df[:1])
# print(df.iloc[0, 5])

def Hij(x, kx, ky, kz):
    R = df.iloc[x, 0]*r1+df.iloc[x, 1]*r2+df.iloc[x, 2]*r3
    K = np.array([kx, ky, kz])
    return df.iloc[x, 5]*np.exp(-1j*np.dot(R, K))


def H(kx, ky, kz):
    A = zeros((25, 25), complex)
    for x in range (len(df)):
        i = df.iloc[x, 4]-1
        j = df.iloc[x, 3]-1
        A[i][j] += Hij(x, kx, ky, kz)
    return A
def eigvals(kx, ky, kz):
    return sp.linalg.eigh(H(kx, ky, kz))[0]


N = 2
q_ary = linspace(0, 1, N)
fig, ax = plt.subplots()

for l in range(25):
    val = zeros(N)
    offset = zeros(8)
    # Γ-M
    i = 0
    for q in q_ary:
        K = q*g1/2+0*g2+0*g3
        val[i] = eigvals(K[0], K[1], K[2])[l]
        i = i+1
    offset[0] = 0
    plot(sum(offset)+q_ary*norm(g1/2), val, '#377eb8')
    # M-K
    i = 0
    for q in q_ary:
        K = (1/2-q/6)*g1+q*g2/3+0*g3
        val[i] = eigvals(K[0], K[1], K[2])[l]
        i = i+1
    offset[1] = norm(g1/2)
    plot(sum(offset)+q_ary*(norm(g1/6-g2/3)), val, '#377eb8')
    # K-Γ
    i = 0
    for q in q_ary:
        K = (1/3-q/3)*g1+(1/3-q/3)*g2+0*g3
        val[i] = eigvals(K[0], K[1], K[2])[l]
        i = i+1
    offset[2]= norm(g1/6-g2/3)
    plot(sum(offset)+q_ary*norm(g1/3+g2/3), val, '#377eb8')
    # Γ-A
    i = 0
    for q in q_ary:
        K = 0*g1+0*g2+q*g3/2
        val[i] = eigvals(K[0], K[1], K[2])[l]
        i = i+1
    offset[3] = norm(g1/3+g2/3)
    plot(sum(offset)+q_ary*norm(g3/2), val, '#377eb8')
    # A-L
    i = 0
    for q in q_ary:
        K = q*g1/2+0*g2+g3/2
        val[i] = eigvals(K[0], K[1], K[2])[l]
        i = i+1
    offset[4] = norm(g3/2)
    plot(sum(offset)+q_ary*norm(g1/2), val, '#377eb8')
    # L-H
    i = 0
    for q in q_ary:
        K = (1/2-q/6)*g1+q*g2/3+g3/2
        val[i] = eigvals(K[0], K[1], K[2])[l]
        i = i+1
    offset[5] = norm(g1/2)
    plot(sum(offset)+q_ary*norm(g1/6-g2/3), val, '#377eb8')
    # H-A
    i = 0
    for q in q_ary:
        K = (1/3-q/3)*g1+(1/3-q/3)*g2+g3/2
        val[i] = eigvals(K[0], K[1], K[2])[l]
        i = i+1
    offset[6] = norm(g1/6-g2/3)
    plot(sum(offset)+q_ary*norm(g1/3+g2/3), val, '#377eb8')
    offset[7] = norm(g1/3+g2/3)

ax.set_xlim([0, sum(offset)])
ax.set_ylim([-1.5, 1])
plt.xticks(
    [sum(offset[0:1]), sum(offset[0:2]), sum(offset[0:3]), sum(offset[0:4]), sum(offset[0:5]), sum(offset[0:6]), sum(offset[0:7]), sum(offset)],
    ["Γ", "M", "K", "Γ", "A", "L", "H", "A"])
for i in range(8):
    plt.vlines(sum(offset[0:i+1]), -1.5, 1, "k", linewidth = 0.25)
plt.hlines(0, 0, sum(offset), "m", linestyle = ":")
ax.set_title("Band Structure of GdCo5")
print(float(time.time())-start)
plt.show()

