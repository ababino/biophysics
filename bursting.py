from scipy.integrate import ode, odeint
import numpy as np
import pylab as pl


def  LR(x, t):
    v = x[0]  # potencial de membrana
    n = x[1]  # Subunidad n
    Ca = 0*x[2] # calcio

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
    y = np.zeros((3,))
    y[0] = -( gCa*m_inf**3*h_inf*(v-VCa) + (gK*n**4+gKCa*Ca/(Ca+Kd))*(v-VK) + gL*(v-VL)) / Cm
    y[1] = An*(1-n)-Bn*n
    y[2] = f*(-k1*gCa*m_inf**3*h_inf*(v-VCa)-kc*Ca)

    y = y*1000
    return y


time = np.linspace(0, 20, 2000)
ci = [-50, 0.7, 0.0530]   #estos son los valores estacionarios
#ci = [vini, 0.7, 0.0530]   #estos son los valores estacionarios
x = odeint(LR, ci, time)

v = x[:, 0] # variable v
n = x[:, 1] # variable w
Ca = x[:, 2] # variable w

pl.figure(1)
pl.plot(time, v)
pl.xlabel('tiempo')
pl.ylabel('V')
"""
pl.figure(2)
pl.plot(time, Ca)
pl.xlabel('tiempo')
pl.ylabel('Ca')
"""

pl.show()
