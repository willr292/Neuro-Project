import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from quantities import Hz, s, ms
from corr_spike_trains import correlated_spikes

ee_cor = 0.1
ei_cor = 0.1
ii_cor = 0.1



total_n = 10
n_e = int(4/5 * total_n)
n_i = int(1/5 * total_n)
rate = 1

ee = np.full((n_e, n_e), ee_cor)
ei = np.full((n_i, n_e), ei_cor)
ie = np.transpose(ei)
ii = np.full((n_i, n_i), ii_cor)

A = np.concatenate((ee, ie), axis=1)
B = np.concatenate((ei, ii), axis=1)
print(A)
print('------------------')
print(B)
C = np.concatenate((np.concatenate((ee, ie), axis=0),np.concatenate((ei, ii), axis=0)), axis = 1)
np.fill_diagonal(C, (np.ones(total_n) * rate))

print(C)