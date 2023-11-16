import qutip as qt
import matplotlib.pyplot as plt 
import numpy as np
import time
from oprts import*
from stt_qtt import*

start_time = time.time()

#setting initial conditions 
# initial state

tmin = 0.0
tmax = 200.0
nsteps = 100000

tspan = np.linspace(tmin, tmax, nsteps)

#periodic detuning
sw_itm = '1 + np.sign(np.sin({}*t))'
sw = sw_itm.format(omega_d)

#two particle case: Hamiltonian
Ht = [H0, [Hv, sw]]

#operators for calculating the expectation value
e_ops = [magx, magy, magz]

#jump operators 
j_ops = [np.sqrt(gamma)*magm]

#integrator Lindblad equation
res_lind = qt.mesolve(Ht, rho0, tspan, j_ops, e_ops)

magx_t = res_lind.expect[0]/N
magy_t = res_lind.expect[1]/N
magz_t = res_lind.expect[2]/N

plt.plot(tspan, magx_t)
plt.plot(tspan, magy_t)
plt.plot(tspan, magz_t)

print("--- %s seconds ---" % (time.time() - start_time))

plt.show()