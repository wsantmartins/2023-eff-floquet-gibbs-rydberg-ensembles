import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

#Reading the files from my local computer
data_paths = [
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=3_w=1_w0=10_Nf=2_V=0.1_beta=20',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=3_w=1_w0=10_Nf=2_V=0.15_beta=20',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=3_w=1_w0=10_Nf=2_V=0.2_beta=20',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=3_w=1_w0=10_Nf=2_V=0.25_beta=20',
]

# Define the custom color
custom_colors = ['#c6dbef', '#6baed6', '#2171b5', '#08306b']


# Set the font family and size
font_style = {'font.family': 'Times New Roman', 'font.serif': ['Times New Roman'], 'font.size': 24}
plt.rcParams.update(font_style)
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

# Initialize lists to store data for each dataset
data_sets = []

for data_path in data_paths:
    try:
        with open(data_path, 'r') as file:
            # Read the entire contents of the file into a variable
            data = file.readlines()
    except FileNotFoundError:
        print(f"File '{data_path}' not found.")
        continue

    # Initialize lists to store t and y values for this dataset
    t_values = []
    y_values = []

    # Process each line in the data
    for line in data:
        # Split each line into two values based on whitespace
        values = line.strip().split()
        if len(values) == 2:  # Ensure there are two values in the line
            t, y = map(complex, values)
            t_values.append(t.real)
            y_values.append(y.real)

    # Convert lists to NumPy arrays and store them as a tuple in the data_sets list
    data_sets.append((np.array(t_values), np.array(y_values)))

# Define the frequency omega and calculate the period T
omega = 1  # Example value, set it to your desired frequency
T = 2 * np.pi / omega

# Debug print the period T
print(f"Period T: {T}")


#plot fidelities
plt.figure(figsize=(5,4))

# Set the thickness of the axes lines
plt.gca().spines['top'].set_linewidth(1)
plt.gca().spines['bottom'].set_linewidth(1)
plt.gca().spines['left'].set_linewidth(1)
plt.gca().spines['right'].set_linewidth(1)

# Set the thickness of tick marks
plt.tick_params(axis='both', which='both', width=1, labelsize=24)

# Set y-axis limits
plt.xlim(0, 5e4)
plt.ylim(0.5, 1.1)

# Plot label
plt.xlabel(r'$\omega t$')
#plt.ylabel(r'$F(\varrho^{\prime}_{\mathrm{S}}, \varrho^{\prime 2}_\mathrm{F}) $')
plt.ylabel('fidelities')

for i, (t_values, y_values) in enumerate(data_sets):
    color = custom_colors[i]  # Get the custom color for this dataset

    # Select data points stroboscopically
    stroboscopic_indices = np.isclose((t_values / T) % 1, 0, atol=1e-3)
    stroboscopic_t_values = t_values[stroboscopic_indices]
    stroboscopic_y_values = y_values[stroboscopic_indices]

    # Debug print the stroboscopic values
    print(f"Dataset {i}:")
    print(f"Stroboscopic t values: {stroboscopic_t_values}")
    print(f"Stroboscopic y values: {stroboscopic_y_values}")

    # Plot the stroboscopic data
    plt.plot(stroboscopic_t_values, stroboscopic_y_values, linewidth=4, color=color)


plt.tight_layout()

#Add a legend inside the plot for the specified labels
plt.text(40000, 0.8, '(e)', fontsize = 24)

# Save the figure as a PDF file
plt.savefig('fide.pdf',format='pdf')

# Display the plot
plt.show()
