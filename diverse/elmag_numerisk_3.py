import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange

def create_wire(L, N):
    return np.linspace(-L, L, N)[:-1]

def calc_magneticfield(x, y, L):
    dL = L[1] - L[0]
    B = np.zeros_like(x)
    n = x.shape[0]

    for i in range(len(wire)):
        r = np.sqrt((L[i]-x)**2+y**2)
        B += dL*y/r**3
    B[B>4] = 4
    B[B<-4] = -4
    return B

L = 1
n = 100
x = np.linspace(-2, 2, n)
xi, yi = np.meshgrid(x, x)
wire = create_wire(L, n)
B = calc_magneticfield(xi, yi, wire)
plt.imshow(B, 'hot')
plt.show()
