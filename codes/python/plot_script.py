import matplotlib.pyplot as plt
from floquet_magnus import*

number = np.arange(1, N_tot + 1, 1)

op_norm = np.zeros((N_tot, N_tot))

for i in range(0, N_tot):
    N = number[i]
    for j in range(0, N):
        omega_d = 5.0 + (j+1)*1.0
        tau = 2.0*np.pi/omega_d
        #print(omega_d)
        #print(tau)
        P1 = tau*(X + Y)
        Rt = 0.0
        for k in range(0, N):
            Rt += R(k+1)
        U_t = ((tau*X).expm())*((tau*Y).expm())        
        op_norm[i,j] = (Rt.expm() - U_t).norm()    


# Create a sample matrix (replace this with your actual data)
matrix = op_norm

# Create a figure and axis
fig, ax = plt.subplots()

# Set aspect ratio to 'equal' before adding the colorbar
ax.set_aspect('equal')

# Define the x and y axis limits
x_min, x_max, y_min, y_max = 0, 25, 5, 55

# Plot the matrix as an image
cax = ax.imshow(matrix, origin='lower', extent=[x_min, x_max, y_min, y_max], cmap='bone')  # You can choose a different colormap

# Adjust the aspect ratio to change plot dimensions
ax.set_aspect(0.4)  # Modify this value as needed

# Add axis labels with LaTeX formatting
plt.xlabel(r'$n$', fontsize=14)
plt.ylabel(r'$\omega_\mathrm{d}$', fontsize=14)

# Add a colorbar
cbar = plt.colorbar(cax)

# Optionally, you can set a label for the colorbar
cbar.set_label(r'$\left\|U_t - U^{[n]}_t\right\|$', fontsize=14)

# Save the plot as a PDF file
plt.savefig('floquet.pdf', format='pdf')

# Show the plot
plt.show()

#print(Rt.expm() - U_t)
#print((Rt.expm() - U_t).norm())
#print(op_norm)