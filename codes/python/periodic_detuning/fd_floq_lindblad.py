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


#operators and spectral density
a_ops = [[magx, spectra_cb_numerical]]


#operators for calculating the expectation value
e_ops = [magx, magy, magz]

#Floquet Hamiltonian
Rt = 0
for k in range(0, Nf):
    Rt += R(k+1)

H_floquet = 1j*Rt/tau 

#store the states of the system
opts = qt.Options(store_states=True)

#Here we calculate the Bloch-Redfield equation for this system
res_red_floq = qt.brmesolve(H_floquet, rho0, tspan, a_ops, e_ops, use_secular=True, options=opts)

states_at_times = res_red_floq.states

#mean value for observables 
magx_t = res_red_floq.expect[0]/N
magy_t = res_red_floq.expect[1]/N
magz_t = res_red_floq.expect[2]/N

#mean value for effective operators 

plt.plot(tspan, magx_t)
plt.plot(tspan, magy_t)
plt.plot(tspan, magz_t)


print("--- %s seconds ---" % (time.time() - start_time))

magxt = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_floq_magxt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}_V={V}'
magyt = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_floq_magyt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}_V={V}'
magzt = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_floq_magzt_N={N}_w={omega_d}_w0={omega_0}_Nf={Nf}_V={V}'

np.savetxt(magxt, np.c_[tspan, magx_t])
np.savetxt(magyt, np.c_[tspan, magy_t])
np.savetxt(magzt, np.c_[tspan, magz_t])

plt.show()