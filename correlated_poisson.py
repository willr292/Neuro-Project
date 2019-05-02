import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from quantities import Hz, s, ms
from corr_spike_trains import correlated_spikes
import seaborn as sns

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
probability = 0.3
n_proc = 1000
P = np.full((n_proc, n_proc), probability)
C = np.full((n_proc, n_proc), 0)
nu = np.full(n_proc, total_e_rate/probability)

cor_spk = correlated_spikes(C, nu, n_proc)
res = cor_spk.offline_mixture(P, nu, n_src=n_proc, n_trg=1, time=10000)

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
plt.xlim(0, 10000)
plt.xlabel('Time (ms)', fontsize=16)
plt.ylabel('Spike Train Index', fontsize=16)
plt.gca().tick_params(axis='both', which='major', labelsize=14)
plt.show()