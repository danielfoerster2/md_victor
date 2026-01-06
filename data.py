#!/usr/bin/python3

import os

import matplotlib.pyplot    as plt
import numpy                as np
import pandas               as pd
import scipy.constants      as cst

from scipy.stats        import norm

study_energy = os.environ.get("study_energy") == "1"
study_temperature = os.environ.get("study_temperature") == "1"
study_vel_distribution = os.environ.get("study_vel_distribution") == "1"

initial_temperature = float(os.environ.get("initial_temperature"))
final_temperature = float(os.environ.get("final_temperature"))

n_atoms = int(os.environ.get("n_atoms"))

path = 'data.dat'
path_init = "vel_init.dat"
path_final = "vel_final.dat"

try:
    df = pd.read_csv(path, delim_whitespace=True)
    df.columns = ['time', 'E_pot', 'E_kin', 'E_th']
except FileNotFoundError:
    exit(f"Error : file {path} not exists")
except Exception as e:
    exit(f"Error : {e}")

time = df['time'].values
epot = df["E_pot"].values
ekin = df["E_kin"].values
eth = df["E_th"].values
etot = epot + ekin
extend_hamiltonian = epot + ekin + eth

print(f"Final time : {time[-1]} fs")
print(f"Simulation time : {time[-1] - time[0]} fs")
temperature = 2 * ekin / ((3 * n_atoms - 6) * cst.k / cst.e)
target_temperature = (final_temperature - initial_temperature) / time[-1] * time + initial_temperature


if study_vel_distribution:
    try:
        df_init_velocity = pd.read_csv("vel_init.dat", delim_whitespace=True)
        df_init_velocity.columns = ["vx", "vy", "vz"]
    except FileNotFoundError:
        exit(f"Error : file {path_init} not exists")
    except Exception as e:
        exit(f"Error : {e}")
    try:
        df_final_velocity = pd.read_csv("vel_final.dat", delim_whitespace=True)
        df_final_velocity.columns = ["vx", "vy", "vz"]
    except FileNotFoundError:
        exit(f"Error : file {path_final} not exists")
    except Exception as e:
        exit(f"Error : {e}")

    vx_i = df_init_velocity['vx'].values
    vy_i = df_init_velocity['vy'].values
    vz_i = df_init_velocity['vz'].values
    vx_f = df_final_velocity['vx'].values
    vy_f = df_final_velocity['vy'].values
    vz_f = df_final_velocity['vz'].values

    x_min_i = min(vx_i.min(), vy_i.min(), vz_i.min())
    x_max_i = max(vx_i.max(), vy_i.max(), vz_i.max())
    margin = 0.1 * (x_max_i - x_min_i)
    x_min_i -= margin
    x_max_i += margin

    x_min_f = min(vx_f.min(), vy_f.min(), vz_f.min())
    x_max_f = max(vx_f.max(), vy_f.max(), vz_f.max())
    margin = 0.1 * (x_max_f - x_min_f)
    x_min_f -= margin
    x_max_f += margin

    mu_x, std_x = norm.fit(vx_f)
    mu_y, std_y = norm.fit(vy_f)
    mu_z, std_z = norm.fit(vz_f)
    x = np.linspace(x_min_f, x_max_f)
    gauss_curve_x = norm.pdf(x, mu_x, std_x)
    gauss_curve_y = norm.pdf(x, mu_y, std_y)
    gauss_curve_z = norm.pdf(x, mu_z, std_z)

    fig,ax = plt.subplots(3,2, figsize=(14,7))
    ax[0][0].hist(vx_i, bins=50, edgecolor='black', alpha=.5)
    ax[1][0].hist(vy_i, bins=50, edgecolor='black', alpha=.5)
    ax[2][0].hist(vz_i, bins=50, edgecolor='black', alpha=.5)
    ax[0][1].hist(vx_f, bins=50, edgecolor='black', alpha=.5)
    ax[0][1].plot(x, gauss_curve_x, c='r', lw=.7)
    ax[1][1].hist(vy_f, bins=50, edgecolor='black', alpha=.5)
    ax[1][1].plot(x, gauss_curve_y, c='r', lw=.7)
    ax[2][1].hist(vz_f, bins=50, edgecolor='black', alpha=.5)
    ax[2][1].plot(x, gauss_curve_z, c='r', lw=.7)
    ax[0][0].set_ylabel(r"Velocity $\mathrm{v}_x$ [$\AA\cdotp fs^{-1}$]")
    ax[1][0].set_ylabel(r"Velocity $\mathrm{v}_y$ [$\AA\cdotp fs^{-1}$]")
    ax[2][0].set_ylabel(r"Velocity $\mathrm{v}_z$ [$\AA\cdotp fs^{-1}$]")
    ax[0][0].set_title("Distribution of initial velocities", fontsize=12)
    ax[0][1].set_title("Distribution of final velocities", fontsize=12)
    [ax[i][j].grid(alpha=.3) for i in range(len(ax)) for j in range(len(ax[i]))]
    [ax[i][0].set_xlim(x_min_i, x_max_i) for i in range(len(ax))]
    [ax[i][1].set_xlim(x_min_f, x_max_f) for i in range(len(ax))]
    [ax[i][j].tick_params(labelbottom=False) for i in range(len(ax)-1) for j in range(len(ax[i]))]
    plt.show()


if study_energy:
    fig, ax = plt.subplots(5, 1, figsize=(12,10))
    ax[0].plot(time, epot, lw=.7)
    ax[1].plot(time, ekin, lw=.7)
    ax[2].plot(time, etot, lw=.7)
    ax[3].plot(time, eth, lw=.7)
    ax[4].plot(time, extend_hamiltonian, lw=.7)
    [ax[i].grid(alpha=.3) for i in range(len(ax))]
    [ax[i].tick_params(labelbottom=False) for i in range(len(ax)-1)]

    ax[-1].set_xlabel(r"Time [$\mathrm{fs}$]")
    ax[0].set_ylabel(r"Potential Eg [$\mathrm{eV}$]")
    ax[1].set_ylabel(r"Kinetic Eg [$\mathrm{eV}$]")
    ax[2].set_ylabel(r"Total Eg [$\mathrm{eV}$]")
    ax[3].set_ylabel(r"Thermal Eg [$\mathrm{eV}$]")
    ax[4].set_ylabel(r"Extend Hamiltonian [$\mathrm{eV}$]")
    ax[0].set_title("Energy evolution during MD simulation", fontsize=15)

    [ax[i].set_xlim(time.min(), time.max()) for i in range(len(ax))]
    plt.show()


if study_temperature:
    plt.plot(time, temperature, lw=.7, label="Temperature")
    plt.plot(target_temperature, c='r', ls='--', lw=.7, label='Target temperature')
    plt.xlabel(r"Time [$\mathrm{fs}$]")
    plt.ylabel(r"Temperature [$\mathrm{K}$]")
    plt.grid(alpha=.3)
    plt.xlim(time.min(), time.max())
    plt.title("Temperature as a function of time", fontsize=15)
    plt.legend()
    plt.show()
