import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from quantities import Hz, s, ms
# from elephant.spike_train_generation import homogeneous_poisson_process
# from elephant.statistics import isi, cv

# Parameters
C = 5.05826489016877e-10 # Capacitance
R = 64762107.67144212 # Membrane Resistance
G = 1/R # Membrane conductance
e_L = -0.0659410438537598 # Resting Potential

Theta_Inf = 0.01724190991525249 # Instantaneous threshold
deltaTheta_s = 0.0024954940387932347 #  update to the spike-dependent component of the threshold - a_spike
deltaV = 0.024641792569422692 # Voltage Addition following spike
b_s = 77.13772508509386 # Spike induced threshold time constant
deltaI = np.array([-2.2756493847005494e-11, -1.718122169268734e-10]) # Afterspike current amplitudes A_j
k = np.array([[1/0.3333333333333333, 1/0.03333333333333334]]) # Afterspike current time constants
f_v = -0.34237330414075656 # voltage fraction following rest
a = 0.0024954940387932347 # Adaptation index of thereshold
b_v = 80.89140483348395 # Voltage induced threshold time constant
f = np.array([1., 1.]) #

# Time
T = 2.00 # Simulation length
dt = 0.00002 # Step size
time = np.arange(0, T+dt, dt) # Step values array

# Voltage
V = np.empty(len(time)) # array for saving Voltage history
V[0] = e_L + 0.001 # set initial to resting potential

# External current
Ie = np.zeros(len(time))

# Afterspike current
I = np.ndarray(shape=(len(time), 2))
I[0] = [0., 0.] 

# Spike-dependent threshold component
Theta_s = np.empty(len(time))
Theta_s[0] = 0.

# Voltage-dependent threshold component
Theta_v = np.empty(len(time))
Theta_v[0] = 0.

# Exicitatory conductance
g_e = np.empty(len(time))
g_e[0] = 0.
tau_e = R * C
rate_e = 35.

# Inhibitatory conductance
g_i = np.empty(len(time))
g_i[0] = 0.
tau_i = R * C
rate_i = 35.

def poisson_spike_trains(rate, big_t, tau_ref):
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

e_spike_train = poisson_spike_trains(rate_e, T, 0)
i_spike_train = poisson_spike_trains(rate_i, T, 0)

def LIFRASCAT():

    for i in range(1, len(time)):
        dTheta_s = -b_s * Theta_s[i - 1]
        Theta_s[i] = Theta_s[i - 1] + dTheta_s * dt

        dTheta_v = a * (V[i - 1] - e_L) - b_v * (Theta_v[i - 1] - Theta_Inf)
        Theta_v[i] = Theta_v[i - 1] + dTheta_v * dt

        dI = np.multiply(-1 * k, I[i -1])
        I[i,:] = I[i - 1,:] + dI * dt

        dg_e = -g_e[i -1]/tau_e
        g_e[i] = g_e[i - 1] + dg_e * dt

        dg_i = -g_i[i -1]/tau_i
        g_i[i] = g_i[i - 1] + dg_i * dt
        
        dV = 1/C * (np.sum(I[i - 1,:]) + G*(V[i - 1] - e_L) + g_e[i - 1] * (V[i - 1]) + g_i[i - 1] * (V[i - 1] - e_L))
        V[i] = V[i - 1] + dV * dt

        # Check for a spike and update voltage etc.
        if V[i] > Theta_v[i] + Theta_s[i]:
            print("spike")
            I[i] = f * I[i] + deltaI
            V[i] = e_L + f_v * (V[i] - e_L) - deltaV
            Theta_s[i] = Theta_s[i] + deltaTheta_s
        
        # check if time has passed an event time.
        if np.isin(i * dt, e_spike_train):
            
            return
        
    return

LIFRASCAT()
print(V)
plt.plot(time, V)
plt.show()

# spiketrain_list = [homogeneous_poisson_process(rate=rate_e*Hz, t_start=0.0*s, t_stop=100.0*s, as_array=False) for i in range(1)]

# for i, spiketrain in enumerate(spiketrain_list):
#         t = spiketrain.rescale(ms)
#         plt.plot(t, i * np.ones_like(t), 'k.', markersize=2)
# plt.axis('tight')
# plt.xlim(0, 1000)
# plt.xlabel('Time (ms)', fontsize=16)
# plt.ylabel('Spike Train Index', fontsize=16)
# plt.gca().tick_params(axis='both', which='major', labelsize=14)
# plt.show()