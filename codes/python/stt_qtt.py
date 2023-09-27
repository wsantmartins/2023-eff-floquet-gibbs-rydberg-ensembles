#constants and initial states
import qutip as qt
import numpy as np

#defining constants
hbar = 1
c = 1
k_B = 1
cos = np.cos
exp = np.exp

#in this first part of the program we define the integration time and also the initial states
time = np.linspace(0, 20, 1000)
psi0 = qt.tensor(qt.basis(2,0), qt.basis(2,0))

#temperature
T = 1
beta = 1.0/(k_B*T)

#Rabi frequency and Rydberg interaction
bigomega = 1
V = 1
delta = 1

#driving frequency and laser frequency
omega0_d = 100
#omega_0 = 100
#tau = 2.0*np.pi/omega_d

#Rabi frequency and Rydberg interaction
bigomega = 1
V = 1
delta = 1

