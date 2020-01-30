import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from quantities import Hz, s, ms
from corr_spike_trains import correlated_spikes
import seaborn as sns
import pandas as pd
import json

# total_e_rate = 500
# total_i_rate = 500
# probability = 0.3
# n_proc = 1000
# #P = np.random.randint(0, 2, (n_proc, n_proc))
# nu = np.full(n_proc, total_e_rate/probability) #np.random.random(n_proc) * 10 #rate
# rate = 100
# ee_cor = 0.5
# ei_cor = 0.5
# ii_cor = 0.5

# n_e = int(4/5 * n_proc)
# n_i = int(1/5 * n_proc)

# ee = np.full((n_e, n_e), ee_cor)
# ei = np.full((n_i, n_e), ei_cor)
# ie = np.transpose(ei)
# ii = np.full((n_i, n_i), ii_cor)

# C = np.concatenate((np.concatenate((ee, ie), axis=1),np.concatenate((ei, ii), axis=1)), axis = 0)
# np.fill_diagonal(C, (np.zeros(n_proc) * rate))
# P = C

# cor_spk = correlated_spikes(C, np.ones(n_proc) * rate, n_proc)
# res = cor_spk.offline_mixture(P, nu, n_src=n_proc, n_trg=n_proc, time=10000)

# spike_trains = []
# for i in range(n_proc):
#     spike_trains.append([])

# for x in range(len(res)):
#     spike_trains[res[x][1]].append(res[x][0])

# np.asarray(spike_trains)
# for i, spiketrain in enumerate(spike_trains):
#         t = spiketrain#.rescale(ms)
#         plt.plot(t, i * np.ones_like(t), 'k.', markersize=2)
# plt.axis('tight')
# plt.xlim(0, 10000)
# plt.xlabel('Time (ms)', fontsize=16)
# plt.ylabel('Spike Train Index', fontsize=16)
# plt.gca().tick_params(axis='both', which='major', labelsize=14)
# plt.show()

# ee_cor = 8
# ei_cor = 8
# ii_cor = 8

# n_proc = 100
# n_e = int(4/5 * n_proc)
# n_i = int(1/5 * n_proc)
# rate = 3

# ee = np.full((n_e, n_e), ee_cor)
# ei = np.full((n_i, n_e), ei_cor)
# ie = np.transpose(ei)
# ii = np.full((n_i, n_i), ii_cor)

# C = np.concatenate((np.concatenate((ee, ie), axis=1),np.concatenate((ei, ii), axis=1)), axis = 0)
# np.fill_diagonal(C, (np.ones(n_proc) * rate))

# print(C)

# cor_spk = correlated_spikes(C, np.ones(n_proc) * rate, n_proc)
# spikes = cor_spk.cox_process(time=200, tau_c=1, dt=0.1) # time(ms)
# cor_spk.raster_plot()
# plt.show()

# # print(spikes)
# # spike_times = np.empty()
# # for x in range(spikes):
# #     if spikes[x] = 1:
# #         spike_times.append(x * dt)

# n_proc = 5

# C = (np.ones((n_proc, n_proc)) +
#         np.random.uniform(0, 1, (n_proc, n_proc)) * 5.0)
# np.fill_diagonal(C, [5, 6, 7, 8, 9])
# C = np.maximum(C, C.T)
# rates = np.array([5, 15, 4, 6, 7])
# print(C)

# cor_spk = correlated_spikes(C, rates, n_proc)
# spikes = cor_spk.cox_process(time=20000)
# cor_spk.raster_plot()
# spk = cor_spk.extract_pyNCS_list()
# plt.show()

# def corr_spikes(r, t_c, time, dt, n_syn, c):
#     spike_id = []
#     spike_train = []
#     x_save = []
#     x = r
#     sigma = np.sqrt(c * r * t_c)
#     for t in np.arange(0, time, dt):
#         x = x + dt*((r - x)/t_c) + np.sqrt(dt) * (np.random.normal(0, sigma) / t_c)
#         x_save.append(x)
#         p = x * dt
#         z = np.random.uniform(size=n_syn)
#         n_spikes = np.where( z < p, 1, 0)
#         if np.sum(n_spikes) > 0:
#             spike_id.append(np.where(n_spikes == 1)[0])
#             spike_train.append(t)

#     return np.asarray(spike_train), np.asarray(spike_id), np.asarray(x_save)

# (trains, ids, x) = corr_spikes(30, 0.01, 1, 1e-6, 30, 30)
# print(trains)
# print(ids)
# plt.plot(x)
# plt.show()
# plt.plot(trains, ids, '.')
# plt.show()

total_e_rate = 500
total_i_rate = 500
e_probability = 0.3
i_probability = 0.3
n_proc = 2
P = np.full((n_proc, n_proc), e_probability)
C = np.full((n_proc, n_proc), 0)
nu = np.full(n_proc, total_e_rate)

cor_spk = correlated_spikes(C, nu, n_proc)
res = cor_spk.offline_mixture(P, nu, n_src=n_proc, n_trg=1, time=1000)

spike_trains = []
for i in range(n_proc):
    spike_trains.append([])

for x in range(len(res)):
    spike_trains[res[x][1]].append(res[x][0])

np.asarray(spike_trains)
for i, spiketrain in enumerate(spike_trains):
        t = spiketrain#.rescale(ms)
        plt.plot(t, i * np.ones_like(t), 'k.', markersize=2)
plt.axis('tight')
plt.xlim(0, 1000)
plt.xlabel('Time (ms)', fontsize=16)
plt.ylabel('Spike Train Index', fontsize=16)
plt.gca().tick_params(axis='both', which='major', labelsize=14)
#plt.show()


with open('neuronal_model_597605616/neuron_config.json') as json_file:  
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
print(V)

def LIF_R_ASC_AT(we, wi, rateE, rateI, e_prob, i_prob, n_proc, n_trg):
    print(V)
    e_P = np.full((n_proc, n_proc), e_prob)
    i_P = np.full((n_proc, n_proc), i_prob)
    C = np.full((n_proc, n_proc), 0)
    e_nu = np.full(n_proc, rateE)
    i_nu = np.full(n_proc, rateI)
    cor_spk_e = correlated_spikes(C, e_nu, n_proc)
    e_res = cor_spk_e.offline_mixture(e_P, e_nu, n_src=n_proc, n_trg=1, time=1000)

    e_spike_train = []
    for x in range(len(e_res)):
        e_spike_train.append(e_res[x][0])

    cor_spk_i = correlated_spikes(C, i_nu, n_proc)
    i_res = cor_spk_i.offline_mixture(i_P, i_nu, n_src=n_proc, n_trg=1, time=1000)

    i_spike_train = []
    for x in range(len(i_res)):
        i_spike_train.append(i_res[x][0])
    
    np.asarray(e_spike_train)
    np.asarray(i_spike_train)
    print(e_spike_train)
    e_spike_train = np.round_(np.divide(e_spike_train, dt))
    i_spike_train = np.round_(np.divide(i_spike_train, dt))
    print(e_spike_train)
    spikes = []
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
        
        dV = 1/C * ( np.sum(I[i - 1,:]) - G*(V[i - 1] - e_L) - (g_e[i - 1] * (V[i - 1])) - (g_i[i - 1] * (V[i - 1] - e_L)) + 0*Ie[i])
        
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

LIF_R_ASC_AT(w_e, w_i, rate_e, rate_i, 0.3, 0.3, 2, 1)