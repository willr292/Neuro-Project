import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from quantities import Hz, s, ms
from elephant.spike_train_generation import homogeneous_poisson_process, cpp
from elephant.spike_train_correlation import corrcoef, covariance
from elephant.conversion import BinnedSpikeTrain
from elephant.statistics import mean_firing_rate
import seaborn as sns
import pandas as pd
import json
from itertools import permutations
from numpy.polynomial.polynomial import polyfit

files = ['neuronal_model_597605616/neuron_config.json',
'neuronal_model_603320017/neuron_config.json',
'neuronal_model_637930677/neuron_config.json',
'neuronal_model_591249612/neuron_config.json',
'neuronal_model_529894099/neuron_config.json',
'neuronal_model_543890855/neuron_config.json',
'neuronal_model_523003560/neuron_config.json',
'neuronal_model_591249824/neuron_config.json']
out = ['synapse_597605616.eps','synapse_603320017.eps','synapse_637930677.eps','synapse_591249612.eps','synapse_529894099.eps']

with open(files[7]) as json_file:  
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


def LIF_R_ASC_AT(we, wi, e_spike_train, i_spike_train):
    np.asarray(e_spike_train)
    np.asarray(i_spike_train)
    # e_spike_train[:] = np.round_(e_spike_train[:] / dt)
    # i_spike_train[:] = np.round_(i_spike_train[:] / dt)
    e_result = np.round_((e_spike_train[:] * 0.001) / dt)
    i_result = np.round_((i_spike_train[:] * 0.001) / dt)


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
        
        if np.isin(i, e_result):         
            g_e[i] = g_e[i - 1] + we
            
        if np.isin(i, i_result):
            g_i[i] = g_i[i - 1] + wi
        
    return np.asarray(spikes)

# spiketrain_list = cpp(20 * Hz, (0, 0, 1), 1000 * ms)

# for i, spiketrain in enumerate(spiketrain_list):
#         t = spiketrain.rescale(ms)
#         plt.plot(t, i * np.ones_like(t), 'k.', markersize=2)
# plt.axis('tight')
# plt.xlim(0, 1000)
# plt.xlabel('Time (ms)', fontsize=16)
# plt.ylabel('Spike Train Index', fontsize=16)
# plt.gca().tick_params(axis='both', which='major', labelsize=14)
# #plt.show()
# cc_matrix = corrcoef(BinnedSpikeTrain(spiketrain_list, 1 * ms))
# print(cc_matrix[0][1])

rate_correlation = []
for x in permutations(np.divide(np.linspace(0, 100, 11),100), 3):
    if sum(x) == 1:
        spiketrain_list = cpp(500 * Hz, x, 1000 * ms)
        rate = len(LIF_R_ASC_AT(w_e, w_i, spiketrain_list[0], spiketrain_list[1]))
        cc_matrix = corrcoef(BinnedSpikeTrain(spiketrain_list, 5 * ms))
        rate_correlation.append([cc_matrix[0][1], rate])

print(rate_correlation)
x_val = [x[0] for x in rate_correlation]
y_val = [x[1] for x in rate_correlation]
#plt.scatter(x_val, y_val, marker="x")
sns.regplot(x_val, y_val, ci=None)
plt.ylim((0, 30))
plt.xlim((0, 1))
plt.xlabel("Pearson’s correlation coefficient")
plt.ylabel("Output firing rate (Hz)")
#sns.lmplot("Correlation", "Output firing rate (Hz)", pd.DataFrame((x_val, y_val), columns =['Correlation', 'Output firing rate (Hz)']))
plt.show()

slope, intecept = np.polyfit(x_val, y_val, 1)
print(slope)

