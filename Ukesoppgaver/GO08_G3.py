import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

x0 = 1.
F0 = 1.

x = np.linspace(-100*x0, 100*x0, 1001)
x = x[x!=0]
F = F0 * np.sin(x/x0)/(x/x0)

plt.plot(x, F, label='F(x)')
plt.xlabel('Posisjon x [m]')
plt.ylabel('Kraft F [N]')
plt.legend()
plt.show()

W = np.trapz(F, x)
print(f'Arbeid: {W:.2f} J')

U = cumtrapz(-F, x)
plt.plot(x[:-1], U, label='U(x)')
plt.show()
