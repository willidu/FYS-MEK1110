from time import time
t1 = time()
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0,10,100_001) # Jeg valgte å ha 100 001 tidsintervaller ettersom utregningen ikke tok mer enn noen få sekunder og dette sikrer nøyaktige resultater   
a = np.zeros_like(t)
v = np.zeros_like(t)
y = np.zeros_like(t)
A0 = .45
fc = 488
tc = .67
p = 1.293
Cd = 1.2
F = 400
m = 80
fv = 25.8
w = 0
dt = t[1] - t[0]
a[0] = 400/80
v[0] = 0
y[0] = 0

# def A(i):
#     return A0*(1-.25*np.exp(-(t[i]/tc)**2))

# def D(i):
#     return 1/2 * p * Cd * A0*(1-.25*np.exp(-(t[i]/tc)**2)) * (v[i] - w)**2  

# def FV(i):
#     return -fv*v[i]

# def Fc(i):
#     return fc*np.exp(-(t[i]/tc)**2)


for i in range(len(t) - 1):
    a[i] =  (F + fc*np.exp(-(t[i]/tc)**2) -fv*v[i] - 1/2 * p * Cd * A0*(1-.25*np.exp(-(t[i]/tc)**2)) * (v[i] - w)**2  ) / m
    v[i+1] = v[i] + a[i]*dt
    y[i+1] = y[i] + v[i+1]*dt 
    tol = 10E-4 / 2
    if abs(y[i] - 100) < tol:
        print(f'The runner used {t[i]} seconds to run 100m')

plt.plot(t,a,t,v,t,y)
plt.grid()
plt.xlabel('i (s)')
plt.legend(['Acceleration', 'Velocity', 'Position'])
print(time()-t1)
plt.show()
