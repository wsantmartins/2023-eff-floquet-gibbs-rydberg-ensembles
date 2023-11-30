import qutip as qt
import matplotlib.pyplot as plt 
import numpy as np
import time
from oprts import*
from stt_qtt import*
from floquet_magnus import*

start_time = time.time()

#periodic detuning
sw_itm = '0.5 + 0.5*np.sign(np.sin({}*t))'
sw = sw_itm.format(omega_d)

#two particle case: Hamiltonian
Ht = [H0, [Hv, sw]]

#operators for calculating the expectation value
e_ops = [magx, magy, magz]

#jump operators 
j_ops = [np.sqrt(gamma)*magm]

#store the states
opts = qt.Options(store_states=True)

#integrator Lindblad equation
res_lind = qt.mesolve(Ht, rho0, tspan, j_ops, e_ops, options=opts)

Rt = 0
rho_FG = 0
for k in range(0, Nf):
    Rt += R(k+1)

H_floquet = 1j*Rt/tau 

rho_FG = -beta*H_floquet
rho_FG = rho_FG.expm()
rho_FG = rho_FG/np.trace(rho_FG) 

states_at_times = res_lind.states

magx_t = res_lind.expect[0]/N
magy_t = res_lind.expect[1]/N
magz_t = res_lind.expect[2]/N

#mean value for effective operators 
magx_FG = [np.trace(rho_FG*magx)/N]*len(tspan) 
magy_FG = [np.trace(rho_FG*magy)/N]*len(tspan) 
magz_FG = [np.trace(rho_FG*magz)/N]*len(tspan) 

tr_dist = [qt.tracedist(rho_FG, rho_t) for rho_t in states_at_times]
fidlty = [qt.fidelity(rho_FG, rho_t) for rho_t in states_at_times]

#plt.plot(tspan, magx_t)
#plt.plot(tspan, magy_t)
#plt.plot(tspan, magz_t)

#plt.plot(tspan, tr_dist, label = r'$\mathrm{tr}|\varrho_t - \varrho^{(m)}_\mathrm{FG}|$', color ='#b2182b') #red
#plt.plot(tspan, fidlty, label = r'$\mathcal{F}|\varrho_t - \varrho^{(m)}_\mathrm{FG}|$', color ='#2166ac') #blue
plt.plot(tspan, magx_FG)
plt.plot(tspan, magx_t)
plt.plot(tspan, magy_FG)
plt.plot(tspan, magy_t)
plt.plot(tspan, magz_FG)
plt.plot(tspan, magz_t)

print("--- %s seconds ---" % (time.time() - start_time))

plt.show()