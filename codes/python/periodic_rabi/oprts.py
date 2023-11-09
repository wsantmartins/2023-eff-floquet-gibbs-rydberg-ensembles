#script for the operator's definition
import qutip as qt 
from stt_qtt import*

#defining the single-body operators
sx = qt.sigmax()
sy = qt.sigmay()
sz = qt.sigmaz()
sp = qt.sigmap()
sm = qt.sigmam()
sn = sp*sm

# Setup operators for individual qubits
sx_list, sy_list, sz_list, sp_list, sm_list, sn_list = [], [], [], [], [], []

for i in range(N):
    op_list = [qt.qeye(2)]*N
    op_list[i] = sx
    sx_list.append(qt.tensor(op_list))
    op_list[i] = sy
    sy_list.append(qt.tensor(op_list))
    op_list[i] = sz
    sz_list.append(qt.tensor(op_list))
    op_list[i] = sp
    sp_list.append(qt.tensor(op_list))
    op_list[i] = sm
    sm_list.append(qt.tensor(op_list))
    op_list[i] = sn
    sn_list.append(qt.tensor(op_list))

#magnetization terms 
H0 = 0
mag = 0
magx = 0
magy = 0
magz = 0
magp = 0
magm = 0
for i in range(N):
    mag += sn_list[i]
    magx += sx_list[i]
    magy += sy_list[i]
    magz += sz_list[i]
    magp += sp_list[i]
    magm += sm_list[i]

#interaction terms: here snn is the two-body interaction
snn = 0
if N == 1:
    snn = 0.0
else:
    for i in range(N - 1):
        snn += sn_list[i]*sn_list[i+1]

#time-independent part of the Hamiltonian for time dependent Rabi freq.
H0 = delta*mag + V*snn
Hv = magx