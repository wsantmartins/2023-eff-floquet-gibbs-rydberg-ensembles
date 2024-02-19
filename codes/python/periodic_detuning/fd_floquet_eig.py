import qutip as qt
import matplotlib.pyplot as plt 
import numpy as np
import time
from oprts import*
from stt_qtt import*
from floquet_magnus import*

#first we call the calculation of the Floquet Hamiltonian
Rt = 0
for k in range(0, Nf):
    Rt += R(k+1)

H_floquet = 1j*Rt/tau 

#spectral decompisition of the Floquet Hamiltonian; spec_dec[0]([1]) corresponds to eigenvalues(states)
spec_dec = H_floquet.eigenstates()

#eigenenergies and eigenstates
floquet_ee = spec_dec[0]
floquet_es = spec_dec[1] 

#calculating eigenoperator and its respective eigenfrequency 
EG_list = []
EF_list = []

#building the jump operators
for i in range(N):
    for j in range(N):
        
        X_ij = floquet_es[i] * floquet_es[j].dag()
        EG_list.append(X_ij)
        
        mu_ij = (floquet_ee[i] - floquet_ee[j])
        print(mu_ij)
        if i == j:
            gamma_ij = 0  # Skip the current iteration if i is equal to j
        elif mu_ij < 0:
            mu_ij = np.abs(mu_ij)
            gamma_ij = kappa0 * (mu_ij**3.0) * (np.exp(beta * mu_ij)/(1.0 - np.exp(beta * mu_ij)))
            gamma_ij = np.sqrt(gamma_ij)
        else:
            gamma_ij = kappa0 * (mu_ij**3.0) * (1.0 + np.exp(-beta * mu_ij)/(1 - np.exp(-beta * mu_ij)))                
            gamma_ij = np.sqrt(gamma_ij)
        
        EF_list.append(gamma_ij)


#operators for calculating the expectation value
e_ops = []        
e_ops = [magx, magy, magz]

#print(EF_list)

#print(floquet_ee)
#jump operators
j_ops = [] 
j_ops = [EF_list[i]*EG_list[i] for i in range(N**2)]

#store the states
opts = qt.Options(store_states=True)

res_sch = qt.sesolve(H_floquet, psi0, tspan, e_ops)
res_floq_lind = qt.mesolve(H_floquet, rho0, tspan, j_ops, e_ops, options=opts)

#magx_t = res_sch.expect[0]/N
#magy_t = res_sch.expect[1]/N
#magz_t = res_sch.expect[2]/N
print(mu_ij)

magx_t = res_floq_lind.expect[0]/N
magy_t = res_floq_lind.expect[1]/N
magz_t = res_floq_lind.expect[2]/N

plt.plot(tspan, magx_t)
plt.plot(tspan, magy_t)
plt.plot(tspan, magz_t)

plt.show()