from scipy.integrate import odeint
import numpy as np
import pylab as pl

def HH(y, t, I0, t_in, t_end):
    gna = 120.0
    gk = 36.0
    gl = 0.3
    Vna = 115.0
    Vk = -12.0
    Vl = 10.6
    Cm = 1.0
    dydt = np.zeros_like(y)
    V = y[0]
    n = y[1]
    m = y[2]
    h = y[3]
    Iap = 0.0
    if t_in < t < t_end:
        Iap = I0
    alpha_m = 0.1 * (25 - V) / (np.exp((25 - V) / 10) - 1)
    beta_m = 4 * np.exp(-V / 18)
    alpha_h = 0.07 * np.exp(-V / 20)
    beta_h = 1 / (1 + np.exp((30 -V) / 10))
    alpha_n = 0.01 * (10 - V) / (np.exp((10 - V) / 10) - 1)
    beta_n = 0.125 * np.exp(-V / 80)
    term1 = gna * m**3 * h * (V - Vna)
    term2 = gk * n**4 * (V - Vk)
    term3 = gl * (V - Vl)
    dydt[0] = -(term1 + term2 + term3 + Iap)/ Cm
    dydt[1] = alpha_n * (1 - n) - beta_n * n
    dydt[2] = alpha_m * (1 - m) - beta_m * m
    dydt[3] = alpha_h * (1 - h) - beta_h * h
    return dydt

def llinas(y, t, V):
    Vm = V - 70
    k10 = 2.0
    k20 = 1.0
    Z1 = 1.0
    Z2 = 0.0
    R = 0.025
    k1 = k10 * np.exp(R * Z1 * Vm)
    k2 = k20 * np.exp(R * Z2 * Vm)
    dydt = k1 * (1 - y) - k2 * y
    return dydt

def current(O, V):
    Vm = V - 70
    PCa = 10**(-10)
    Ci = 10**(-4)
    Ce = 40
    R = 0.025
    j = PCa * Vm * (Ci -Ce * np.exp(-2 * R * Vm)) / (1 - np.exp(-2 * R * Vm))
    ICa = (j * O**5) / 5
    return ICa



def  LR(x, t):
    v =x[1]  # potencial de membrana
    n =x[2]  # Subunidad n
    Ca = x[3] # calcio

    # Parametros
    gKCa=0.02
    gK=3
    gCa=3.2
    Kd=1
    VK=-75
    VCa=100
    VL=-40
    Cm=.5
    gL=0.012
    va=30
    vp=50
    f=0.007
    kc=0.02
    k1=0.0275

    An=0.01*((10-(v+va))/(np.exp((10-(v+va))*0.1)-1))
    Bn=0.125*np.exp(-(v+va)/80)
    Am=0.1*((25-(v+vp))/(np.exp((25-(v+vp))*0.1)-1))
    Bm=4*np.exp(-(v+vp)/18)
    Ah=0.07*np.exp(-(v+vp)/20)
    Bh=1/(1+np.exp((30-(v+vp))*0.1))

    m_inf=Am/(Am+Bm)
    h_inf=Ah/(Ah+Bh)

    #ecuaiones
    y[1] = -( gCa*m_inf^3*h_inf*(v-VCa) + (gK*n^4+gKCa*Ca/(Ca+Kd))*(v-VK) + gL*(v-VL)) / Cm
    y[2] = An*(1-n)-Bn*n
    y[3] = f*(-k1*gCa*m_inf^3*h_inf*(v-VCa)-kc*Ca)

    #y = y'*1000
    return y

func = lambda y, t: HH(y, t, -5, 2, 4)
y0 = np.array([0.000, 0.3177, 0.0530, 0.5959])
time = np.linspace(0, 20, 1000)
y = odeint(func, y0, time)
pl.figure(1)
pl.plot(time, y[:, 0])

func = lambda y, t: HH(y + np.array([20, 0, 0, 0]), t, -5, 2, 4)
y0 = np.array([0.000, 0.3177, 0.0530, 0.5959])
y = odeint(func, y0, time)
pl.plot(time, y[:, 0])


pl.xlabel('t')
pl.ylabel('V')
pl.grid(True)

pl.show()
