import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from quantities import Hz, s, ms
from elephant.spike_train_generation import homogeneous_poisson_process
from elephant.statistics import mean_firing_rate
import seaborn as sns
import json
import pandas as pd

files = ['neuronal_model_597605616/neuron_config.json',
'neuronal_model_603320017/neuron_config.json',
'neuronal_model_637930677/neuron_config.json',
'neuronal_model_591249612/neuron_config.json',
'neuronal_model_529894099/neuron_config.json']
out = ['synapse_597605616.eps','synapse_603320017.eps','synapse_637930677.eps','synapse_591249612.eps','synapse_529894099.eps']

with open(files[3]) as json_file:  
    data = json.load(json_file)
# Parameters
C = data['C'] #8.587009771685807e-11 # Capacitance
R = data['R_input'] #340045253.17350113 # Membrane Resistance
G = 1/R # Membrane conductance
e_L = data['El_reference'] #-0.07139990234375 # Resting Potential
Theta_Inf = data['th_inf'] #0.01724190991525249 # Instantaneous threshold
deltaTheta_s = data['threshold_reset_method']['params']['a_spike'] #  update to the spike-dependent component of the threshold - a_spike
deltaV = data['voltage_reset_method']['params']['b'] # 0.004055590135223554 # Voltage Addition following spike
b_s = data['threshold_dynamics_method']['params']['b_spike'] # 102.99935176846839 # Spike induced threshold time constant
deltaI = np.array(data['asc_amp_array']) # np.array([-2.2756493847005494e-11, -1.718122169268734e-10]) # Afterspike current amplitudes A_j
k = np.array(data['asc_tau_array']) # np.array([[1/0.3333333333333333, 1/0.03333333333333334]]) # Afterspike current time constants
f_v = data['voltage_reset_method']['params']['a'] # 0.3987646006290668 # voltage fraction following rest
a = data['threshold_dynamics_method']['params']['a_spike'] # 0.0024954940387932347 # Adaptation index of thereshold
b_v = data['threshold_dynamics_method']['params']['b_voltage'] # 34.344467877141945 # Voltage induced threshold time constant
f = np.array(data['AScurrent_reset_method']['params']['r']) # np.array([1., 1.]) #

# Time
T = 1.00 # Simulation length
dt = 1e-4 # Step size
time = np.arange(0, T+dt, dt) # Step values array

# Voltage
V = np.empty(len(time)) # array for saving Voltage history
V[0] = e_L # set initial to resting potential

# External current
Ie = np.zeros(len(time))
#Ie[10] = -0.0000000005 

# Afterspike current
I = np.ndarray(shape=(len(time), 2))
I[0] = [0., 0.] 

# Spike-dependent threshold component
Theta_s = np.empty(len(time))
Theta_s[0] = 0.

# Voltage-dependent threshold component
Theta_v = np.empty(len(time))
Theta_v[0] =  Theta_Inf

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

e_spike_train = np.asarray([homogeneous_poisson_process(rate=rate_e*Hz, t_start=0.0*s, t_stop=1.0*s, as_array=True) for i in range(1)])#poisson_spike_trains(rate_e, T, 0)
i_spike_train = np.asarray([homogeneous_poisson_process(rate=rate_i*Hz, t_start=0.0*s, t_stop=1.0*s, as_array=True) for i in range(1)])#poisson_spike_trains(rate_i, T, 0)
e_spike_train[:] = np.round_(e_spike_train[:] / dt)
i_spike_train[:] = np.round_(i_spike_train[:] / dt)

def LIF_R_ASC_AT(we, wi):
    spikes = []
    for i in range(1, len(time)):
        dTheta_s = -b_s * Theta_s[i - 1]
        Theta_s[i] = Theta_s[i - 1] + dTheta_s * dt

        dTheta_v = a * (V[i - 1] - e_L) - b_v * (Theta_v[i - 1] - Theta_Inf)
        Theta_v[i] = Theta_v[i - 1] + dTheta_v * dt

        dI = np.multiply(-1 * k, I[i -1])
        I[i,:] = I[i - 1,:] + dI * dt

        dg_e = -g_e[i - 1]/tau_e
        g_e[i] = g_e[i - 1] + dg_e * dt

        dg_i = -g_i[i -1]/tau_i
        g_i[i] = g_i[i - 1] + dg_i * dt
        
        dV = 1/C * ( np.sum(I[i - 1,:]) - G*(V[i - 1] - e_L) - (g_e[i - 1] * (V[i - 1])) - (g_i[i - 1] * (V[i - 1] - e_L)) + 0*Ie[i])
        # Why is the inhibitatory reversal potential the same as the resting potential 
        V[i] = V[i - 1] + dV * dt

        # Check for a spike and update voltage etc.  
        if V[i] > Theta_v[i] + Theta_s[i] + e_L:          
            I[i , :] = f * I[i - 1, :] + deltaI
            V[i] = e_L + f_v * (V[i - 1] - e_L) - deltaV
            Theta_s[i] = Theta_s[i - 1] + deltaTheta_s
            
            spikes.append(i*dt)       
        
        # check if time has passed a synapse event time.
        
        if np.isin(i, e_spike_train):         
            g_e[i] = g_e[i - 1] + we
            
        if np.isin(i, i_spike_train):
            g_i[i] = g_i[i - 1] + wi
        
    return np.asarray(spikes)

#spike = LIF_R_ASC_AT(w_e, w_i)

#plt.plot(time, V)
#plt.show()
#print(spike)



# heatmap = []
# for i in np.arange(0, 20e-9, 1e-9):
#     heatmap.append([])
#     for j in np.arange(0, 20e-9, 1e-9):
#         heatmap[int(np.true_divide(i, 1e-9))].append(len(LIF_R_ASC_AT(i,j)))

g_e_vec = np.linspace(0, 20, 21) * 1e-9#np.arange(0, 20e-9, 1e-9)
g_i_vec = np.linspace(0, 20, 21) * 1e-9#np.arange(0, 20e-9, 1e-9)

heatmap2 = np.zeros((len(g_e_vec),len(g_i_vec)))

for i in range(len(g_e_vec)):
    for j in range(len(g_i_vec)):        
        heatmap2[i,j] = len(LIF_R_ASC_AT(g_e_vec[i], g_i_vec[j]))
        if heatmap2[i ,j] > 100:
            heatmap2[i,j] = 100


# print(heatmap2)


sns.heatmap(pd.DataFrame(data = heatmap2,index = g_e_vec, columns = g_i_vec)).invert_yaxis()
# plt.pcolormesh(np.arange(0, 20e-9, 5e-9),np.arange(0, 20e-9, 5e-9),heatmap)
plt.xlabel("Inhibtory synaptic strength")
plt.ylabel("Excitatory synaptic strength")
plt.show()
