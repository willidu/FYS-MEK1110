import numpy as np
import matplotlib.pyplot as plt

m = 5       # kg
k = 500     # N/m
l0 = 0.5    # m
h = 0.3     # m

x = np.linspace(-0.75, 0.75, 100)
Fx = -k*x*(1-l0/np.sqrt(x**2 + h**2))

plt.plot(x, Fx, label='Horizontal spring force')
plt.legend()
plt.show()
