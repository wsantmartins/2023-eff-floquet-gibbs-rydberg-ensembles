import trace
from stt_qtt import*
from oprts import*
from floquet_magnus import*
import qutip as qt 
import matplotlib.pyplot as plt
import numpy as np

#this program implements the Bloch-Redfield equation for two Rydberg ions interacting with a thermal reservoir

#General J factor
gamma = (
    "(4* w**3) / (3*hbar*c**3)"
)

#Term for positive frequencies
gamma_p = (
    gamma
    + "* (1 + exp(-beta*hbar*w)) / \
          (1-exp(-beta*hbar*w))"
)

#Term for negative frequencies
gamma_m = (
    gamma
    + "* exp(beta*hbar*w) / \
            (1-exp(beta*hbar*w))"
)


# define spectra with variable names
spectra_cb = "(w > 0) * " + gamma_p + "+ (w < 0) * " + gamma_m

# add a check for min size of w to avoid numerical problems
spectra_cb = "0 if (w > -1e-4 and w < 1e-4) else " + spectra_cb

# define string with numerical values expect for w
constants = ["hbar", "c", "beta"]
spectra_cb_numerical = spectra_cb
for ct in constants:
    # replace constants with numerical value
    spectra_cb_numerical = spectra_cb_numerical.replace(ct, str(eval(ct)))

em_itm = 'exp(-{}*1j*t)'
ep_itm = 'exp({}*1j*t)'

em = em_itm.format(omega_0)
ep = ep_itm.format(omega_0)

a_ops = [[(sp1 + sp2, sm1 + sm2), (spectra_cb_numerical, ep, em)]]

#operators for calculating the expectation value
e_ops = [sx1, sx1 + sx2, sz1, sz1 + sz2]


#delta = eval(delta)
delta_itm = "cos({}*t)"
delta = delta_itm.format(omega_d)

#This part of the code initialise the Hamiltonian of the system
Hz = sz1 + sz2
H_0 = V*(n1*n2) + bigomega*(sx1 + sx2)

H_t = [H_0, [Hz, delta]]

#store the states of the system
opts = qt.Options(store_states=True)

#Here we calculate the Bloch-Redfield equation for this system
res_brme = qt.brmesolve(H_t, rho0, tspan, a_ops, e_ops, use_secular=False, options=opts)

#Here are the Floquet quasi-energies
args = {'wd': omega_d}
Hf_t = [H_0, [Hz, 'cos(wd*t)']]

#f_modes_0, f_energies = qt.floquet_modes(Hf_t, T, args)

#print(f_energies)

#Floquet-Gibbs states
Rt = 0
rho_FG = 0
for k in range(0, N_tot):
    Rt += R(k+1)

H_floquet = 1j*Rt/tau 

rho_FG = -beta*H_floquet
rho_FG = rho_FG.expm()
rho_FG = rho_FG/np.trace(rho_FG) 

states_at_times = res_brme.states

tr_dist = [qt.tracedist(rho_FG, rho_t) for rho_t in states_at_times]

print(len(tr_dist))
print(rho_FG)

plt.figure()

# Add axis labels with LaTeX formatting
plt.xlabel(r'$t$', fontsize=14)
plt.ylabel(r'$\mathrm{tr}|\varrho_t - \varrho^{(m)}_\mathrm{FG}|$', fontsize=14)

plt.plot(tspan, tr_dist)

plt.show()