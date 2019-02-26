import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# PARAMETERS
G = 0 # Membrane conductance
C = 3 # Capacitance
R = 1 # Membrane Resistance
E_L = 1 # Resting Potential

Theta_Inf = 0 # Instantaneous threshold

deltaV = 0 # Voltage Addition following spike
b_s = 0 # Spike induced threshold time constant
deltaI_j = 0 # Afterspike current amplitudes
k_j = 0 # Afterspike current time constants
f_v = 1 # Current fraction following spike
a = 0 # Adaptation index of thereshold
b_v = 0 # Voltage induced threshold time constant


# Time
T = 0.100 # Simulation length
dt = 0.00002 # Step size
time = np.arange(0, T+dt, dt) # Step values array

# Voltage
V = np.empty(len(time)) # array for saving Voltage history
V[0] = E_L # set initial to resting potential

# External current
Ie = np.zeros(len(time))
# add some external current stuff

# Afterspike current
Ij = np.empty(len(time))
Ij[0] = 0 

# Spike-dependent threshold component
Theta_s = np.empty(len(time))
Theta_s[0] = 0

# Voltage-dependent threshold component
Theta_v = np.empty(len(time))
Theta_v[0] = 0


def LIFRASCAT():

    for i in range(1, len(time)):
        dTheta_s = -b_s * Theta_s[i - 1]
        Theta_s[i] = Theta_s[i - 1] + dTheta_s * dt

        dTheta_v = a * (V[i - 1] - E_L) - b_v * (Theta_v[i - 1] - Theta_Inf)
        Theta_v[i] = Theta_v[i - 1] + dTheta_v * dt

        dIj = -k_j * Ij[i -1]
        Ij[i] = I[j - 1] + dIj * dt
        
        dV = 1/C * (I[i - 1] + G*(V[i - 1] - E_L))
        V[i] = V[i - 1] + dV * dt

        # Check for a spike and update voltage etc.
        if V[i] > Theta_v[i] + Theta_s[i]:
            Ij[i] = 0
            V[i] = E_L + f_v * (V[i] - E_L) - deltaV

