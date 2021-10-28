import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import interpolate
import sys

q_1_fine = pd.read_csv('src/output/M_1_not_yaw.csv',header=None).to_numpy()
q_1_fine = np.reshape(q_1_fine,[29,1001],order='F')
time = np.linspace(0,100,1001)
print(q_1_fine)
plt.plot(time[:-1],q_1_fine[0,:-1])
# plt.xlim([90,100])
plt.savefig('results.png')

plt.figure()
q_init = pd.read_csv('src/input/M_init.csv',header=None).to_numpy()
q_init = np.reshape(q_init,[2,31],order='F')
print(q_init.shape)
plt.plot(q_init[1,1:-1],'.')
plt.savefig('shape.png')