

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate

df = pd.read_csv('data/turbine_structure_UAE.csv')
radius = df['r'].to_numpy()
rho = df['rho'].to_numpy()
EI_edge = df['ES'].to_numpy()
EI_flap = df['FS'].to_numpy()

plt.plot(radius,rho)

N_total = 32#int(input("Please enter the ALM node number:\n"))
r_nacelle = 0.508#float(input("Please enter the nacelle radius:\n"))

def find_nearest(a, a0):
    idx = np.abs(a - a0).argmin()
    return a.flat[idx]

dr = 5.5/N_total
alm_points = np.linspace(0,5.5,N_total+1)[1:] - r_nacelle

location = np.where(alm_points == find_nearest(alm_points,0))

print('Starting point:',location[0])

f = interpolate.interp1d(radius, rho)
rho_alm = f(alm_points)
f = interpolate.interp1d(radius, EI_edge)
EI_e = f(alm_points)
f = interpolate.interp1d(radius, EI_flap)
EI_f = f(alm_points)

df = pd.read_csv('data/turbine_geo_uae.csv')
radius = df['radius'].to_numpy()
chord_length = df['chord length'].to_numpy()
twist_angle = df['twist angle'].to_numpy()

f = interpolate.interp1d(radius, chord_length)
chord_length_alm = f(alm_points)
f = interpolate.interp1d(radius, twist_angle)
twist_angle_alm = f(alm_points)
twist = twist_angle_alm/180*np.pi

np.savetxt("interp_data/chord_length.csv", chord_length_alm, delimiter=",",fmt='%1.8f')
np.savetxt("interp_data/twist_angle.csv", twist, delimiter=",",fmt='%1.8f')
np.savetxt("interp_data/edgewise_stiffness.csv", EI_e, delimiter=",",fmt='%1.8f')
np.savetxt("interp_data/flapwise_stiffness.csv", EI_f, delimiter=",",fmt='%1.8f')
np.savetxt("interp_data/density.csv", rho_alm, delimiter=",",fmt='%1.8f')
np.savetxt("interp_data/alm_node.csv", alm_points+r_nacelle, delimiter=",",fmt='%1.8f')