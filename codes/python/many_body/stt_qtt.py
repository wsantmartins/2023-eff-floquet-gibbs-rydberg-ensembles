#constants and initial states
import qutip as qt
import numpy as np

#defining constants
#hbar = 1
#c = 10
#k_B = 1
cos = np.cos
exp = np.exp


#temperature
T = 1
#beta = 1.0/(k_B*T)
beta = 1.0/T

#Rabi frequency and Rydberg interaction
#Everything is in units of delta
bigomega = 1
V = 1
delta = -0.5

#decay rate and dissipation constant
gamma = 0.1
kappa0 = gamma/(np.sqrt(1 + 4.0*bigomega)**3)

#driving frequency and laser frequency
omega_d = 4
omega_0 = 100
tau = 2.0*np.pi/omega_d
#tau = 2.0*np.pi/omega_d

N = 4
N_tot = N

#single ion case: initial states
state_list = [qt.basis(2, 0)]*N
psi0 = qt.tensor(state_list)
rho0 = psi0*psi0.dag()



