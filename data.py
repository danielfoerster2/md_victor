import os

import matplotlib.pyplot    as plt
import numpy                as np
import pandas               as pd
import scipy.constants      as cst


path = 'Data/data.csv'
df = pd.read_csv(path)

#df.info()
#print(df.describe())
#print(df.head())

n_atoms = df['n_atoms'][0]
t = df['time']
E_pot = df["E_pot"]
E_kin = df["E_kin"]
E_th = df["E_th"]
E = df["E_tot"]
T = 2 * E_kin / (3 * n_atoms * cst.k/cst.e)


fig, ax = plt.subplots(4, 1, figsize=(10,8))
ax[0].plot(t, E_pot, lw=.7)
ax[1].plot(t, E_kin, lw=.7)
ax[2].plot(t, E_th, lw=.7)
ax[3].plot(t, E, lw=.7)
[ax[i].grid(alpha=.3) for i in range(len(ax))]
[ax[i].tick_params(labelbottom=False) for i in range(len(ax)-1)]

ax[-1].set_xlabel("Time [fs]")
ax[0].set_ylabel("Potential Eg [eV]")
ax[1].set_ylabel("Kinetic Eg [eV]")
ax[2].set_ylabel("Thermostat Eg [eV]")
ax[3].set_ylabel("Total Eg [eV]")
ax[0].set_title("Energy evolution during MD simulation", fontsize=15)

[ax[i].set_xlim(t.min(), t.max()) for i in range(len(ax))]
plt.savefig("Curve/energy.png")
plt.show()


plt.plot(t, T, lw=.7)
plt.xlabel("Time [fs]")
plt.ylabel("Temperature [K]")
plt.grid(alpha=.3)
plt.xlim(t.min(), t.max())
plt.title("Temperature as a function of time", fontsize=15)
plt.savefig("Curve/temperature.png")
plt.show()
