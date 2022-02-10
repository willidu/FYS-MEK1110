import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(-8, 8, 41)
x, y = np.meshgrid(t, t)

def h(x, y, h0=1200, a=250, R=2):
    return h0 + a/(R**2) * x * y

h_ = h(x, y)

plt.contour(x, y, h_)
plt.show()
