# 2023ryd_eng_repo
The codes contained in this folder refer to the work Quasiperiodic Floquet-Gibbs states in Rydberg atomic systems, Wilson S. Martins, Federico Carollo, Kay Brandner and Igor Lesanovsky.

All codes are written in Python 3.11.5, and the libraries required to run the codes are Matplotlib, NumPy and QuTip, in the following versions:

QuTiP Version: 4.7.3
Numpy Version: 1.25.2
Matplotlib Version: 3.7.2

All the files containing "fd_" generate and/or display data, and they are: 
- fd_redfield.py, that simulates the Redfield equation for a chain of spins using pre-defined noise-power spectrum, time-dependent Hamiltonian and time-dependent interaction terms;
- fd_plt_script_strob_fid and fd_plt_script_strob_pol, that generate the plots for the stroboscopic fidelities and polarizations, contained in the main text.

The script floquet_magnus.py generates a Floquet-Magnus expantion for a given Hamiltonian at arbitrary order, exploring the recursive structure of the BCH formula expansion. Also, stt_qtts and oprts contain variable entries and operators used to construct the Hamiltonian. 

Data and figures are contained in the namesake folders. 