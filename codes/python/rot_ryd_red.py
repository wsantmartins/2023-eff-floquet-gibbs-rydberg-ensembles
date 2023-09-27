#from qutip import*
from stt_qtt import*
from oprts import*
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

em = 'exp(-100j*t)'
ep = 'exp(100j*t)'

#expm = eval(expm)
#expp = eval(expp)

a_ops = [[(sp1 + sp2, sm1 + sm2), (spectra_cb_numerical, ep, em) ] ]

#operators for calculating the expectation value
e_ops = [sx1, sx1 + sx2, sz1, sz1 + sz2]


#delta = eval(delta)
#args = {'w': omega_d}
delta = "cos(2*t)"

#This part of the code initialise the Hamiltonian of the system
Hz = sz1 + sz2
H_0 = V*(n1*n2) + bigomega*(sx1 + sx2)

H_t = [H_0, [Hz, delta]]

#Here we calculate the Bloch-Redfield equation for this system
res_brme = qt.brmesolve(H_t, psi0, time, a_ops, e_ops)

#Here are the Floquet quasi-energies
args = {'w': omega0_d}
Hf_t = [H_0, [Hz, 'cos(w*t)']]

T = 2*np.pi/omega0_d
f_modes_0, f_energies = qt.floquet_modes(Hf_t, T, args)

#print(f_energies)

#Floquet-Gibbs states
rho_fg = qt.qdiags(f_energies, 0)
rho_fg = -beta*rho_fg
rho_fg = rho_fg.expm() 

print(rho_fg)

plt.figure()

plt.plot(time, res_brme.expect[0], label=r'$\sigma_x$')
plt.plot(time, res_brme.expect[1], label=r'$\sigma^1_x + \sigma^2_x$')
plt.plot(time, res_brme.expect[2], label=r'$\sigma_z$')
plt.plot(time, res_brme.expect[3], label=r'$\sigma^1_z + \sigma^2_z$')

plt.legend()

plt.show()