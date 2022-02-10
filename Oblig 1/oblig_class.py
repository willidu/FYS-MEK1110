import numpy as np
import matplotlib.pyplot as plt

class Sprint:
    def __init__(self, rho, cd, a0, w, m, F, dt, x0, v0, acc0):
        self.rho = rho
        self.cd = cd
        self.a0 = a0
        self.w = w
        self.m = m
        self.F = F
        self.dt = dt
        self.n = int(np.ceil(10/dt))
        self.x0 = x0
        self.v0 = v0
        self.acc0 = acc0

    def calculate_100m(self):
        t = np.linspace(0, 10, self.n)
        x = np.zeros(self.n); x[0] = self.x0
        v = np.zeros(self.n); v[0] = self.v0
        a = np.zeros(self.n); a[0] = self.acc0
 
        for i in range(self.n-1):
            if x[i-1] > 100:
                self.i_max = i
                t = t[:i]; x = x[:i]; v = v[:i]; a = a[:i]
                break

            a[i] = (2*self.F - self.rho*self.cd*self.a0*(v[i]**2))/(2*self.m)
            v[i+1] = v[i] + a[i] * self.dt
            x[i+1] = x[i] + v[i+1] * self.dt
            t = np.append(t, (i+1)*self.dt)
        
        return t, x, v, a
    
    def plot(self, t, x, v, a):
        plt.clf()

        plt.subplot(311)
        plt.plot(t, x, 'r', label='x(t)')
        plt.ylabel('x [m]')
        plt.ylim(0)
        plt.xlim(0)

        plt.subplot(312)
        plt.plot(t, v, 'g', label='v(t)')
        plt.ylabel('v [ms^-1]')
        plt.ylim(0)
        plt.xlim(0)

        plt.subplot(313)
        plt.plot(t, a, 'b', label='a(t)')
        plt.ylabel('a [ms^-2]')
        plt.xlabel('t [s]')
        plt.xlim(0)

        plt.tight_layout()

    def terminal_velocity(self):
        t, x, v, a = self.calculate_100m()
        self.dt = 1  # Scaling dt to reduce computing time
        i = len(t)-1
        while a[i] > 0.01:
            a = np.append(a, (2*self.F - self.rho*self.cd*self.a0*(v[i]**2))/(2*self.m))
            v = np.append(v, v[i] + a[i] * self.dt)
            x = np.append(x, x[i] + v[i+1] * self.dt)
            t = np.append(t, (i+1)*self.dt)
            i += 1

        return v[-1]

class Sprint2(Sprint):
    def __init__(self, rho, cd, a0, w, m, F, dt, x0, v0, acc0, fc, tc, fv):
        super().__init__(rho, cd, a0, w, m, F, dt, x0, v0, acc0)
        self.fc = fc
        self.tc = tc
        self.fv = fv

    def calculate_100m(self):
        t = np.array([0])
        x = np.zeros(self.n); x[0] = self.x0
        v = np.zeros(self.n); v[0] = self.v0
        a = np.zeros(self.n); a[0] = self.acc0

        for i in range(self.n-1):
            if x[i-1] > 100:
                self.i_max = i
                t = t[:i]; x = x[:i]; v = v[:i]; a = a[:i]
                break

            D = -0.5*self.rho*self.cd*self.a0*(1-0.25*np.exp(-(t[i]/self.tc)**2))*(v[i]-self.w)**2
            F_net = self.F + self.fc*np.exp(-(t[i]/self.tc)**2) - self.fv*v[i] + D
            
            a[i] = (F_net/self.m)
            v[i+1] = v[i] + a[i] * self.dt
            x[i+1] = x[i] + v[i+1] * self.dt
            t = np.append(t, (i+1)*self.dt)
        
        return t, x, v, a

class Forces(Sprint2):
    def forces(self):
        t = np.array([0])
        x = np.zeros(self.n); x[0] = self.x0
        v = np.zeros(self.n); v[0] = self.v0
        a = np.zeros(self.n); a[0] = self.acc0

        F = np.full(self.n, self.F)
        Fc = np.zeros(self.n)
        Fv = np.zeros(self.n)
        D = np.zeros(self.n)

        for i in range(self.n-1):
            if x[i-1] > 100:
                self.i_max = i
                t = t[:i]; F = F[:i]; Fc = Fc[:i]; Fv = Fv[:i]; D = D[:i]
                break
            
            Fc[i] = self.fc*np.exp(-(t[i]/self.tc)**2)
            Fv[i] = -self.fv*v[i]
            D[i] = -0.5*self.rho*self.cd*self.a0*(1-0.25*np.exp(-(t[i]/self.tc)**2))*(v[i]-self.w)**2

            F_net = self.F + Fc[i] + Fv[i] + D[i]
            
            a[i] = (F_net/self.m)
            v[i+1] = v[i] + a[i] * self.dt
            x[i+1] = x[i] + v[i+1] * self.dt
            t = np.append(t, (i+1)*self.dt)
        
        return t, F, Fc, Fv, D

rho = 1.293     # kg/m^3
cd = 1.2
a0 = 0.45       # m^2
w = 0           # m/s
m = 80          # kg
F = 400         # N
x0 = 0          # m
v0 = 0          # m/s
acc0 = 0        # m/s^2


bolt = Sprint(rho, cd, a0, w, m , F, 0.01, x0, v0, acc0)
t, x, v, a = bolt.calculate_100m()
bolt.plot(t, x, v, a)
plt.show()
print(f'Time to run {x[-1]:.2f}m is {t[-1]:.2f}s')

vt = bolt.terminal_velocity()
print(f'Terminal velocity: {vt:.2f}m/s')

fc = 488       # N
tc = 0.67      # s
fv = 25.8      # Ns/m

bolt2 = Sprint2(rho, cd, a0, w, m , F, 0.01, x0, v0, acc0, fc, tc, fv)
t, x, v, a = bolt2.calculate_100m()
bolt2.plot(t, x, v, a)
plt.show()
print(f'Time to run {x[-1]:.2f}m is {t[-1]:.2f}s')

bolt3 = Forces(rho, cd, a0, w, m , F, 0.01, x0, v0, acc0, fc, tc, fv)
t, F_, Fc, Fv, D = bolt3.forces()
plt.plot(t, F_, t, Fc, t, Fv, t, D)
plt.legend(['F - konstant kraft', 'Fc - sammenkrypning', 'Fv - fysisk begrensning', 'D - luftmotstand'], loc='best', bbox_to_anchor=(0.5, 0.3, 0.5, 0.5))
plt.xlim(0)
plt.ylabel('Force [N]')
plt.xlabel('Time [s]')
plt.show()

for wind in (-1.0, 1.0):
    runner = Sprint2(rho, cd, a0, wind, m , F, 0.01, x0, v0, acc0, fc, tc, fv)
    t, x, v, a = runner.calculate_100m()
    print(f'Time to run {x[-1]:.2f}m with {wind}m/s wind is {t[-1]:.2f}s')
