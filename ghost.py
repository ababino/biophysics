from scipy.integrate import odeint, ode
import numpy as np
import pylab as pl
from models import dendrite, soma, nmda, HH, LR




def nmda_impulse(y, t, T=10, t_in=20, t_end=23):
    rb = 0.005
    dy = nmda(y)
    if t_in <= t < t_end:
        dy[0] += -rb * T * y[0]
        dy[1] += rb * T * y[0] - rb * y[1] * T
        dy[2] += rb * y[1]
    return dy

LRt = lambda t, y: LR(y)
y0 = np.array([-50.0, 0.7, 0.0530])
#time = np.linspace(0, 500, 1000000)
solver = 'dop853'#'vode'#'dop853'
r = ode(LRt).set_integrator(solver)
r.set_initial_value(y0, 0)
t1 = 20
dt = 0.01
v = []
time = []
while r.t < t1:#r.successful() and
    r.integrate(r.t+dt)
    v.append(r.y[0])
    time.append(r.t)
#out = odeint(HH_input, y0, time)
pl.figure(1)
pl.plot(time, v)
pl.grid(True)
pl.show()
"""

def HH_input(t, y):
    #y[0] -= y[7]*2*10**4
    Mg = 1
    gNMDA = 10**4
    dy = np.zeros_like(y)
    dy[:4] = HH(y[:4])
    dy[4:] = nmda_impulse(y[4:], t)
    B = 1.0 / (1 + np.exp(-0.062*y[0])*Mg/3.57)
    ICa = gNMDA * B * (y[7])
    dy[0] += ICa
    return dy


y0 = np.array([-65.000,  0.3177, 0.0530, 0.5959, 1, 0, 0, 0])
#time = np.linspace(0, 500, 1000000)
solver = 'vode'#'dop853'
r = ode(HH_input).set_integrator(solver)
r.set_initial_value(y0, 0)
t1 = 500
dt = 0.01134579
I = []
v = []
time = []
while r.t < t1:#r.successful() and
    r.integrate(r.t+dt)
    v.append(r.y[0])
    I.append(r.y[7])
    time.append(r.t)
#out = odeint(HH_input, y0, time)
pl.figure(1)
pl.plot(time, v)
pl.grid(True)
pl.figure(2)
pl.plot(time, I)
pl.grid(True)
pl.show()



y0 = np.array([1, 0, 0, 0])
time = np.linspace(0, 1000, 10000)
out = odeint(nmda_impulse, y0, time)
pl.figure(1)
pl.plot(time, out[:, 3])
pl.grid(True)
pl.show()

def ghost(y, t, I0=9, t_in=0, t_end=2000):
    gc = 1
    kappa = 0.4
    dy = np.zeros_like(y)
    dy[:2] = soma(y[:2])
    dy[2:] = dendrite(y[2:])
    if t_in <= t < t_end:
        dy[0] += I0
    dy[0] += gc * (y[2] - y[0]) / kappa
    dy[2] += gc * (y[0] - y[2]) / (1 - kappa)
    return dy

y0 = np.array([ -70.0, 4.5*10**(-5), -70.0, 0.97, 0.0024754, 0.69686])
time = np.linspace(0, 650, 650000)
out = odeint(ghost, y0, time)
pl.figure(1)
pl.plot(time, out[:, 0])
pl.xlim([400, 650])
pl.grid(True)

pl.figure(2)
pl.plot(time, out[:, 5])
pl.xlim([400, 650])
pl.ylim([0, 0.3])
pl.grid(True)

pl.show()
"""
