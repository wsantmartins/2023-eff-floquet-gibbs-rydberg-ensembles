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


#plt.plot(tspan, tr_dist, label = r'$\mathrm{tr}|\varrho_t - \varrho^{(m)}_\mathrm{FG}|$', color ='#b2182b') #red
#plt.plot(tspan, fidlty, label = r'$\mathcal{F}|\varrho_t - \varrho^{(m)}_\mathrm{FG}|$', color ='#2166ac') #blue
plt.plot(tspan, magx_FG)
plt.plot(tspan, magx_t)
plt.plot(tspan, magy_FG)
plt.plot(tspan, magy_t)
plt.plot(tspan, magz_FG)
plt.plot(tspan, magz_t)

fid_file = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\ll_f_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
td_file = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\ll_td_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magxt = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\ll_magxt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magxfg = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\ll_magxfg_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magyt = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\ll_magyt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magyfg = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\ll_magyfg_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magzt = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\ll_magzt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magzfg = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\ll_magzfg_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'

np.savetxt(fid_file, np.c_[tspan, fidlty])
np.savetxt(td_file, np.c_[tspan, tr_dist])
np.savetxt(magxt, np.c_[tspan, magx_t])
np.savetxt(magyt, np.c_[tspan, magy_t])
np.savetxt(magzt, np.c_[tspan, magz_t])
np.savetxt(magxfg, np.c_[tspan, magx_FG])
np.savetxt(magyfg, np.c_[tspan, magy_FG])
np.savetxt(magzfg, np.c_[tspan, magz_FG])

print("--- %s seconds ---" % (time.time() - start_time))

plt.show()