#constants and initial states
import qutip as qt
import numpy as np

#defining constants
#hbar = 1
#c = 10
#k_B = 1
cos = np.cos
exp = np.exp

#in this first part of the program we define the integration time and also the initial states
tspan = np.linspace(0, 500, 100000)

#single ion case: initial states
psi0 = qt.basis(2,0)
rho0 = psi0*psi0.dag()

#two ions case: initial states
#psi0 = qt.tensor(qt.basis(2,0), qt.basis(2,0))
#rho0 = psi0*psi0.dag()

#temperature
T = 1
#beta = 1.0/(k_B*T)
beta = 1.0/T

#Rabi frequency and Rydberg interaction
#Everything is in units of delta
bigomega = 1
V = 1
delta = 1

#decay rate and dissipation constant
gamma = 0.1
kappa0 = gamma/(np.sqrt(1 + 4.0*bigomega)**3)

#driving frequency and laser frequency
omega_d = 8
omega_0 = 100
tau = 2.0*np.pi/omega_d
#tau = 2.0*np.pi/omega_d

N_tot = 21


