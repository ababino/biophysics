import numpy as np


def cable(V, I=0):
    dV = np.zeros_like(V)
    N = dV.shape[0]
    dV[1:N-1] = V[:N-2] - 3 * V[1:N-1] + V[2:]
    dV[0] = V[1] + dV[1] + I - V[0]
    dV[N-1] = V[N-1] - V[N-2] - dV[N-2]
    return dV


def na_current(V, m, h):
    gna = 120.0
    Vna = 115.0
    INa = gna * m**3 * h * (V - Vna)
    return INa


def HH(y):
    gk = 36.0
    gl = 0.3
    Vk = -12.0
    Vl = 10.6
    Cm = 1.0
    dy = np.zeros_like(y)
    V = y[0]
    n = y[1]
    m = y[2]
    h = y[3]
    alpha_m = 0.1 * (25.0 - V) / (np.exp((25.0 - V) / 10.0) - 1.0)
    beta_m = 4.0 * np.exp(-V / 18.0)
    alpha_h = 0.07 * np.exp(-V / 20.0)
    beta_h = 1.0 / (1.0 + np.exp((30.0 - V) / 10.0))
    alpha_n = 0.01 * (10.0 - V) / (np.exp((10.0 - V) / 10.0) - 1.0)
    beta_n = 0.125 * np.exp(-V / 80.0)
    INa = na_current(V, m, h)
    IK = gk * n**4 * (V - Vk)
    IL = gl * (V - Vl)
    dy[0] = -(INa + IK + IL) / Cm
    dy[1] = alpha_n * (1 - n) - beta_n * n
    dy[2] = alpha_m * (1 - m) - beta_m * m
    dy[3] = alpha_h * (1 - h) - beta_h * h
    return dy


def soma(y):
    gna = 55.0
    gk = 20.0
    gl = 0.18
    Vna = 40.0
    Vk = -88.5
    Vl = -70.0
    Cm = 1.0
    V12 = -40
    K = 3
    tau = 0.39
    dy = np.zeros_like(y)
    V = y[0]
    n = y[1]
    m_inf = 1 / (1 + np.exp(-(V - V12) / K))
    n_inf = 1 / (1 + np.exp(-(V - V12) / K))
    term1 = gna * m_inf**2 * (1 - n) * (Vna - V)
    term2 = gk * n**2 * (Vk - V)
    term3 = gl * (Vl - V)
    dy[0] = (term1 + term2 + term3)/ Cm
    dy[1] = (n_inf - n) / tau
    return dy


def dendrite(y):
    gna = 5.0
    gk = 15.0
    gl = 0.18
    Vna = 40.0
    Vk = -88.5
    Vl = -70.0
    Cm = 1.0
    V12 = -40
    V12h = -52
    V12p = -65
    K = 5
    Kh = -5
    Kp = -6
    tau_h = 1
    tau_n = 0.9
    tau_p = 5
    dy = np.zeros_like(y)
    V = y[0]
    h = y[1]
    n = y[2]
    p = y[3]
    m_inf = 1 / (1 + np.exp(-(V - V12) / K))
    n_inf = 1 / (1 + np.exp(-(V - V12) / K))
    h_inf = 1 / (1 + np.exp(-(V - V12h) / Kh))
    p_inf = 1 / (1 + np.exp(-(V - V12p) / Kp))
    INa = gna * m_inf**2 * h * (Vna - V)
    IK = gk * n**2 * p * (Vk - V)
    Il = gl * (Vl - V)
    dy[0] = (INa + IK + Il)/ Cm
    dy[1] = (h_inf - h) / tau_h
    dy[2] = (n_inf - n) / tau_n
    dy[3] = (p_inf - p) / tau_p
    return dy


def nmda(y):
    ru = .0129
    rd = 0.0084
    rr = 0.0068
    r0 = 0.0465
    rc = 0.0738
    C1 = y[1]
    C2 = y[2]
    O = y[3]
    D = 1 - sum(y)
    dy = np.zeros_like(y)
    dy[0] = ru * C1
    dy[1] = -ru * C1 + ru * C2
    dy[2] = -ru * C2 -rd * C2 +rr * D -r0 * C2 + rc * O
    dy[3] = r0 * C2 - rc * O
    return dy

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



def  LR(y, ICa=0):
    v = y[0]  # potencial de membrana
    n = y[1]  # Subunidad n
    Ca = y[2] # calcio

    # Parametros
    gKCa = 0.02
    gK = 3.0
    gCa = 3.2
    Kd = 1.0
    VK = -75.0
    VCa = 100.0
    VL = -40.0
    Cm = 0.5
    gL = 0.012
    va = 30.0
    vp = 50.0
    f = 0.007
    kc = 0.02
    k1 = 0.0275

    An = 0.01 * ((10.0 - (v + va)) / (np.exp((10.0 - (v + va)) * 0.1) - 1.0))
    Bn = 0.125 * np.exp(-(v + va) / 80.0)
    Am = 0.1 * ((25.0 - (v + vp)) / (np.exp((25.0 - (v + vp)) * 0.1) - 1.0))
    Bm = 4.0 * np.exp(-(v + vp) / 18.0)
    Ah = 0.07 * np.exp(-(v + vp) / 20.0)
    Bh = 1.0 / (1.0 + np.exp((30.0 - (v + vp)) * 0.1))

    m_inf = Am / (Am + Bm)
    h_inf = Ah / (Ah + Bh)

    ICa += gCa * m_inf**3 * h_inf * (v - VCa)
    IK = (gK * n**4 + gKCa * Ca / (Ca + Kd)) * (v - VK)
    IL = gL * (v - VL)

    #ecuaiones
    dy = np.zeros_like(y)
    dy[0] = -( ICa + IK + IL) / Cm
    dy[1] = An * (1 - n) - Bn * n
    dy[2] = f * (-k1 * ICa - kc * Ca)

    #dy = dy * 1000
    return dy


def  LR_ICa(y, ICa):
    v = y[0]  # potencial de membrana
    n = y[1]  # Subunidad n
    Ca = y[2] # calcio

    # Parametros
    gKCa = 0.02
    gK = 3.0
    Kd = 1.0
    VK = -75.0
    VL = -40.0
    Cm = 0.5
    gL = 0.012
    va = 30.0
    f = 0.007
    kc = 0.02
    k1 = 0.0275

    An = 0.01 * ((10.0 - (v + va)) / (np.exp((10.0 - (v + va)) * 0.1) - 1.0))
    Bn = 0.125 * np.exp(-(v + va) / 80.0)

    IK = (gK * n**4 + gKCa * Ca / (Ca + Kd)) * (v - VK)
    IL = gL * (v - VL)

    #ecuaiones
    dy = np.zeros_like(y)
    dy[0] = -( ICa + IK + IL) / Cm
    dy[1] = An * (1 - n) - Bn * n
    dy[2] = f * (-k1 * ICa - kc * Ca)

    #paso el tiempo a segundos
    #dy = dy * 1000
    return dy
