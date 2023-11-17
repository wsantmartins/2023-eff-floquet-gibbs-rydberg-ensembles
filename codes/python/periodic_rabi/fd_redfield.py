import qutip as qt
import matplotlib.pyplot as plt 
import numpy as np
import time
from oprts import*
from stt_qtt import*
from floquet_magnus import*

start_time = time.time()

#This program implements the Bloch-Redfield equation for a many-body Rydberg ions interacting with a thermal reservoir. In each step we have defined operators for a single and two particles 

#General J factor
gamma = (
    "kappa0*w**3"
)

#Term for positive frequencies
gamma_p = (
    gamma
    + "* (1 + exp(-beta*w)) / \
          (1-exp(-beta*w))"
)

#Term for negative frequencies
gamma_m = (
    gamma
    + "* exp(beta*w) / \
            (1-exp(beta*w))"
)


# define spectra with variable names
spectra_cb = "(w > 0) * " + gamma_p + "+ (w < 0) * " + gamma_m

# add a check for min size of w to avoid numerical problems
spectra_cb = "0 if (w > -1e-4 and w < 1e-4) else " + spectra_cb

# define string with numerical values expect for w
#constants = ["hbar", "c", "beta"]
constants = ["kappa0", "beta"]
spectra_cb_numerical = spectra_cb
for ct in constants:
    # replace constants with numerical value
    spectra_cb_numerical = spectra_cb_numerical.replace(ct, str(eval(ct)))

em_itm = 'exp(-{}*1j*t)'
ep_itm = 'exp({}*1j*t)'

em = em_itm.format(omega_0)
ep = ep_itm.format(omega_0)

#single particle case: operators and spectral density
a_ops = [[(magp, magm), (spectra_cb_numerical, ep, em)]]


#operators for calculating the expectation value
e_ops = [magx, magy, magz]


#periodic detuning
sw_itm = "0.5 + 0.5*np.sign(np.sin({}*t))"

sw = sw_itm.format(omega_d)

#two particle case: Hamiltonian
Ht = [H0, [Hv, sw]]

#store the states of the system
opts = qt.Options(store_states=True)

#Here we calculate the Bloch-Redfield equation for this system
res_brme = qt.brmesolve(Ht, rho0, tspan, a_ops, e_ops, use_secular=True, options=opts)

#Floquet-Gibbs states
Rt = 0
rho_FG = 0
for k in range(0, Nf):
    Rt += R(k+1)

H_floquet = 1j*Rt/tau 

rho_FG = -beta*H_floquet
rho_FG = rho_FG.expm()
rho_FG = rho_FG/np.trace(rho_FG) 

states_at_times = res_brme.states

#mean value for observables 
magx_t = res_brme.expect[0]/N
magy_t = res_brme.expect[1]/N
magz_t = res_brme.expect[2]/N

#mean value for effective operators 
magx_FG = [np.trace(rho_FG*magx)/N]*len(tspan) 
magy_FG = [np.trace(rho_FG*magy)/N]*len(tspan) 
magz_FG = [np.trace(rho_FG*magz)/N]*len(tspan) 

tr_dist = [qt.tracedist(rho_FG, rho_t) for rho_t in states_at_times]
fidlty = [qt.fidelity(rho_FG, rho_t) for rho_t in states_at_times]

print(len(tr_dist))
print(rho_FG)

plt.figure()

# Add axis labels with LaTeX formatting
plt.xlabel(r'$t$', fontsize=14)
#plt.ylabel(r'$\mathrm{tr}|\varrho_t - \varrho^{(m)}_\mathrm{FG}|$', fontsize=14)

plt.plot(tspan, tr_dist, label = r'$\mathrm{tr}|\varrho_t - \varrho^{(m)}_\mathrm{FG}|$', color ='#b2182b')
plt.plot(tspan, fidlty, label = r'$\mathcal{F}|\varrho_t - \varrho^{(m)}_\mathrm{FG}|$', color ='#2166ac')
#plt.plot(tspan, magx_FG)
#plt.plot(tspan, magx_t)
#plt.plot(tspan, magy_FG)
#plt.plot(tspan, magy_t)
#plt.plot(tspan, magz_FG)
#plt.plot(tspan, magz_t)

print("--- %s seconds ---" % (time.time() - start_time))


fid_file = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_rabi\\data\\f_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
td_file = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_rabi\\data\\td_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magxt = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_rabi\\data\\magxt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magxfg = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_rabi\\data\\magxfg_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magyt = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_rabi\\data\\magyt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magyfg = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_rabi\\data\\magyfg_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magzt = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_rabi\\data\\magzt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'
magzfg = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_rabi\\data\\magzfg_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}'

np.savetxt(fid_file, np.c_[tspan, fidlty])
np.savetxt(td_file, np.c_[tspan, tr_dist])
np.savetxt(magxt, np.c_[tspan, magx_t])
np.savetxt(magyt, np.c_[tspan, magy_t])
np.savetxt(magzt, np.c_[tspan, magz_t])
np.savetxt(magxfg, np.c_[tspan, magx_FG])
np.savetxt(magyfg, np.c_[tspan, magy_FG])
np.savetxt(magzfg, np.c_[tspan, magz_FG])

plt.show()

