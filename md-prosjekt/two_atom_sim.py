import numpy as np
import matplotlib.pyplot as plt

def a(r, t):
    ri = r[0]
    rj = r[1]
    r_vec = ri-rj
    r_norm = np.linalg.norm(r_vec)

    term = 24*(2*r_norm**-12 - r_norm**-6)
    vec = term / (r_norm**2) * r_vec
    return np.asarray([vec, -vec])

def Euler(T, dt, r0, v0):
    n = int(np.ceil(T/dt))
    t = np.linspace(0, T, n)
    x = np.zeros((n, 2))
    v = np.zeros_like(x)
    x[0] = r0; v[0] = v0

    for i in range(n-1):
        v[i+1] = v[i] + a(x[i], t[i]) * dt
        x[i+1] = x[i] + v[i] * dt

    return t, x, v

def EulerCromer(T, dt, r0, v0):
    n = int(np.ceil(T/dt))
    t = np.linspace(0, T, n)
    x = np.zeros((n, 2))
    v = np.zeros_like(x)
    x[0] = r0; v[0] = v0

    for i in range(n-1):
        v[i+1] = v[i] + a(x[i], t[i]) * dt
        x[i+1] = x[i] + v[i+1] * dt

    return t, x, v

def VelocityVerlet(T, dt, r0, v0):
    n = int(np.ceil(T/dt))
    t = np.linspace(0, T, n)
    x = np.zeros((n, 2))
    v = np.zeros_like(x)
    x[0] = r0; v[0] = v0

    a_ = a(x[0], t[0])
    for i in range(n-1):
        x[i+1] = x[i] + v[i]*dt + 0.5*a_*dt**2
        a_2 = a(x[i+1], t[i+1])
        v[i+1] = v[i] + 0.5*(a_ + a_2)*dt
        a_ = a_2

    return t, x, v

def U_marked(r):
    ri = r[0]
    rj = r[1]
    r_norm = abs(np.linalg.norm(rj-ri))
    return 4*(r_norm**(-12) - r_norm**(-6))
