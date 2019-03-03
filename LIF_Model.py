import numpy as np
import random as rnd
import matplotlib.pyplot as plt

# Parameters
G = 0 # Membrane conductance
C = 3 # Capacitance
R = 1 # Membrane Resistance
E_L = 1 # Resting Potential

Theta_Inf = 0 # Instantaneous threshold
deltaTheta_s =0 #  update to the spike-dependent component of the threshold
deltaV = 0 # Voltage Addition following spike
b_s = 0 # Spike induced threshold time constant
deltaI_j = 0 # Afterspike current amplitudes A_j
k_j = 0 # Afterspike current time constants
f_v = 1 # voltage fraction following rest
a = 0 # Adaptation index of thereshold
b_v = 0 # Voltage induced threshold time constant
f_j = 0 #

# Time
T = 0.100 # Simulation length
dt = 0.00002 # Step size
time = np.arange(0, T+dt, dt) # Step values array

# Voltage
V = np.empty(len(time)) # array for saving Voltage history
V[0] = E_L # set initial to resting potential

# External current
Ie = np.zeros(len(time))

# Afterspike current
Ij = np.empty(len(time))
Ij[0] = 0 

# Spike-dependent threshold component
Theta_s = np.empty(len(time))
Theta_s[0] = 0

# Voltage-dependent threshold component
Theta_v = np.empty(len(time))
Theta_v[0] = 0

# Exicitatory conductance
g_e = np.empty(len(time))
g_e[0] = 0

# Inhibitatory conductance
g_i = np.empty(len(time))
g_i[0] = 0

def spike_trains(rate, big_t, tau_ref):
    if 1 <= rate * tau_ref:
        print("firing rate not possible given refractory period f/p")
        return []

    exp_rate = rate / (1 - tau_ref * rate)

    spike_train = []

    t = rnd.expovariate(exp_rate)

    while t < big_t:
        spike_train.append(t)
        t += tau_ref + rnd.expovariate(exp_rate)

    return spike_train



def LIFRASCAT():

    for i in range(1, len(time)):
        dTheta_s = -b_s * Theta_s[i - 1]
        Theta_s[i] = Theta_s[i - 1] + dTheta_s * dt

        dTheta_v = a * (V[i - 1] - E_L) - b_v * (Theta_v[i - 1] - Theta_Inf)
        Theta_v[i] = Theta_v[i - 1] + dTheta_v * dt

        dIj = -k_j * Ij[i -1]
        Ij[i] = Ij[i - 1] + dIj * dt

        dg_e = = -g_e[i -1]/tau_e
        Ij[i] = Ij[i - 1] + dIj * dt


        
        dV = 1/C * (Ij[i - 1] + G*(V[i - 1] - E_L) + g_e[i - 1] * (V[i - 1]) + g_i[i - 1] * (V[i - 1] - E_L))
        V[i] = V[i - 1] + dV * dt

        # Check for a spike and update voltage etc.
        if V[i] > Theta_v[i] + Theta_s[i]:
            Ij[i] = f_j * Ij[i] + deltaI_j
            V[i] = E_L + f_v * (V[i] - E_L) - deltaV
            Theta_s[i] = Theta_s[i] + deltaTheta_s
        
        # check if time has passed an event time.
        if 

    return

