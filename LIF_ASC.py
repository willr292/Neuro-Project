import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from quantities import Hz, s, ms
from elephant.spike_train_generation import homogeneous_poisson_process
from elephant.statistics import mean_firing_rate
import seaborn as sns
import pandas as pd
import json


with open('neuronal_model_529894099/neuron_config.json') as json_file:  
    data = json.load(json_file)

# Parameters
C = data['C'] #8.587009771685807e-11 # Capacitance
R = data['R_input'] #340045253.17350113 # Membrane Resistance
G = 1/R # Membrane conductance
e_L = data['El_reference'] #-0.07139990234375 # Resting Potential
Theta_Inf = data['th_inf'] #0.01724190991525249 # Instantaneous threshold
deltaI = np.array(data['asc_amp_array']) # np.array([-2.2756493847005494e-11, -1.718122169268734e-10]) # Afterspike current amplitudes A_j
k = np.array(data['asc_tau_array']) # np.array([[1/0.3333333333333333, 1/0.03333333333333334]]) # Afterspike current time constants
R = np.array(data['AScurrent_reset_method']['params']['r']) # np.array([1., 1.]) #
V_r = e_L

# Time
T = 1.00 # Simulation length
dt = 1e-4 # Step size
time = np.arange(0, T+dt, dt) # Step values array

# Voltage
V = np.empty(len(time)) # array for saving Voltage history
V[0] = e_L # set initial to resting potential

# Afterspike current
I = np.ndarray(shape=(len(time), 2))
I[0] = [0., 0.] 

# Exicitatory conductance
g_e = np.empty(len(time))
g_e[0] = 0
tau_e = 5e-3
rate_e = 1000. # 1000 hz max as multiple synapses  
w_e = 8e-9

# Inhibitatory conductance
g_i = np.empty(len(time))
g_i[0] = 0
tau_i = 5e-3
rate_i = 1000.
w_i = 8e-9

# Spike array


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



def LIF_ASC(we, wi, rateE, rateI):
    e_spike_train = np.asarray([homogeneous_poisson_process(rate=rateE*Hz, t_start=0.0*s, t_stop=1.0*s, as_array=True) for i in range(1)])#poisson_spike_trains(rate_e, T, 0)
    i_spike_train = np.asarray([homogeneous_poisson_process(rate=rateI*Hz, t_start=0.0*s, t_stop=1.0*s, as_array=True) for i in range(1)])#poisson_spike_trains(rate_i, T, 0)
    e_spike_train[:] = np.round_(e_spike_train[:] / dt)
    i_spike_train[:] = np.round_(i_spike_train[:] / dt)
    spikes = []
    for i in range(1, len(time)):

        dI = np.multiply(-1 * k, I[i -1])
        I[i,:] = I[i - 1,:] + dI * dt

        dg_e = -g_e[i -1]/tau_e
        g_e[i] = g_e[i - 1] + dg_e * dt

        dg_i = -g_i[i -1]/tau_i
        g_i[i] = g_i[i - 1] + dg_i * dt
        
        dV = 1/C * ( np.sum(I[i - 1,:]) - G*(V[i - 1] - e_L) - (g_e[i - 1] * (V[i - 1])) - (g_i[i - 1] * (V[i - 1] - e_L)))
        V[i] = V[i - 1] + dV * dt

        # Check for a spike and update voltage etc.  
        if V[i] > Theta_Inf + e_L:          
            I[i , :] = R * I[i - 1, :] + deltaI
            V[i] = V_r        
            spikes.append(i*dt)       
        
        # check if time has passed a synapse event time.
        
        if np.isin(i, e_spike_train):         
            g_e[i] = g_e[i - 1] + we
            
        if np.isin(i, i_spike_train):
            g_i[i] = g_i[i - 1] + wi
        
    return np.asarray(spikes)

# spike = LIF_ASC(w_e, w_i, rate_e, rate_i)
# plt.xlabel("Time (Seconds)")
# plt.ylabel("Membrane Potential (Volts)")
# plt.title("LIF-ASC model of a neuron")
# plt.plot(time, V)
# plt.show()
# print(spike)



# heatmap = []
# for i in np.arange(0, 20e-9, 1e-9):
#     heatmap.append([])
#     for j in np.arange(0, 20e-9, 1e-9):
#         heatmap[int(np.true_divide(i, 1e-9))].append(len(LIF_R_ASC_AT(i,j)))

# g_e_vec = np.arange(0, 20e-9, 1e-9)
# g_i_vec = np.arange(0, 20e-9, 1e-9)
# heatmap2 = np.zeros((len(g_e_vec),len(g_i_vec)))

# for i in range(len(g_e_vec)):
#     for j in range(len(g_i_vec)):        
#         heatmap2[i,j] = len(LIF_R_ASC_AT(g_e_vec[i], g_i_vec[j]))


# print(heatmap2)

rate_e_vec = np.arange(0, 1050, 50)
rate_i_vec = np.arange(0, 1050, 50)
heatmap2 = np.zeros((len(rate_e_vec),len(rate_i_vec)))

for i in range(len(rate_e_vec)):
    for j in range(len(rate_i_vec)):        
        heatmap2[i,j] = len(LIF_ASC(w_e, w_i, rate_e_vec[i], rate_i_vec[j]))

plt.plot(rate_i_vec, heatmap2[10, :])
plt.xlabel("Inhibtory synapse firing rate (Hz)")
plt.ylabel("Output firing rate (Hz)")
plt.show()

plt.plot(rate_e_vec, heatmap2[:, 10])
plt.xlabel("Excitatory synapse firing rate (Hz)")
plt.ylabel("Output firing rate (Hz)")
plt.show()

slope_i, intercept_i = np.polyfit(rate_i_vec, heatmap2[10, :], 1)
slope_e, intercept_e = np.polyfit(rate_e_vec, heatmap2[:, 10], 1)
print(slope_i)
print(slope_e)

sns.heatmap(pd.DataFrame(heatmap2, index=rate_i_vec, columns=rate_e_vec), cbar_kws={'label': 'Output firing rate (Hz)'}).invert_yaxis()
plt.xlabel("Inhibtory rate (Hz)")
plt.ylabel("Excitatory rate (Hz)")

plt.show()
