#constants and initial states
import qutip as qt
import numpy as np

#defining constants
hbar = 1
c = 10
k_B = 1
cos = np.cos
exp = np.exp

#in this first part of the program we define the integration time and also the initial states
tspan = np.linspace(0, 400, 100000)
psi0 = qt.tensor(qt.basis(2,0), qt.basis(2,0))
rho0 = psi0*psi0.dag()

#temperature
T = 1
beta = 1.0/(k_B*T)

#Rabi frequency and Rydberg interaction
bigomega = 1
V = 1
delta = 1

#driving frequency and laser frequency
omega_d = 40
omega_0 = 100
tau = 2.0*np.pi/omega_d
#tau = 2.0*np.pi/omega_d

N_tot = 25


