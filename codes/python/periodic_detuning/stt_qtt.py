#constants and initial states
import qutip as qt
import numpy as np

#defining constants
cos = np.cos
exp = np.exp

#setting initial conditions 
# initial state

tmin = 0.0
tmax = 5e5
nsteps = 1000000

tspan = np.linspace(tmin, tmax, nsteps)

#temperature 
beta = 20

#Rabi frequency and Rydberg interaction
#Everything is in units of delta
bigomega = 0.1
V = 0.1
delta = 0.1

#decay rate and dissipation constant

kappa0 = 0.02

#driving frequency and laser frequency
omega_d = 1
omega_0 = 10
tau = 2.0*np.pi/omega_d

N = 2
Nf = 2

#single ion case: initial states
state_list = [qt.basis(2, 0)]*N
psi0 = qt.tensor(state_list)
rho0 = psi0*psi0.dag()



