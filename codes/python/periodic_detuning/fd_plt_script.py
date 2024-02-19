import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert

# Specify the paths to your data files with double backslashes
#data_paths = [
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magxt_N=2_w=8_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magxfg_N=2_w=8_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magyt_N=2_w=8_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magyfg_N=2_w=8_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magzt_N=2_w=8_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magzfg_N=2_w=8_w0=100_Nf=2',
#]

#data_paths = [
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=2_w=6_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=3_w=6_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=4_w=6_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=5_w=6_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=6_w=6_w0=100_Nf=2',
#]

#data_paths = [
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=2_w=4_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=2_w=6_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=2_w=8_w0=100_Nf=2',
#    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_f_N=2_w=12_w0=100_Nf=2',
#]

data_paths = [
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magxt_N=6_w=8_w0=100_Nf=2',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magxfg_N=6_w=8_w0=100_Nf=2',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magyt_N=6_w=8_w0=100_Nf=2',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magyfg_N=6_w=8_w0=100_Nf=2',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magzt_N=6_w=8_w0=100_Nf=2',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\red_magzfg_N=6_w=8_w0=100_Nf=2',
]



# Define the custom color
#custom_colors = ['#252525', '#252525','#252525', '#252525','#252525']
custom_colors = ['#d9d9d9', '#252525','#a6bddb', '#034e7b', '#fcbba1', '#cb181d']

# Set the font family and size
font_style = {'font.family': 'Times New Roman', 'font.serif': ['Times New Roman'], 'font.size': 20}

# Apply the custom font style to the Matplotlib configuration
plt.rcParams.update(font_style)

# Specify LaTeX preamble to set font family for LaTeX-rendered text
#latex_preamble = r'\usepackage{mathptmx}'  # Change to mathptmx for Times New Roman
#plt.rcParams['text.latex.preamble'] = latex_preamble

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

# Set the figure size to be smaller
#plt.figure(figsize=(5, 5))  # Adjust the width and height as needed

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

    # Initialize lists to store x and y values for this dataset
    t_values = []
    x_values = []

    # Process each line in the data
    for line in data:
        # Split each line into two values based on whitespace
        values = line.strip().split()
        if len(values) == 2:  # Ensure there are two values in the line
            x, y = map(complex, values)
            t_values.append(x.real)
            x_values.append(y.real)

    # Convert lists to NumPy arrays and store them as a tuple in the data_sets list
    data_sets.append((t_values, x_values))

# Set the thickness of the axes lines
plt.gca().spines['top'].set_linewidth(1.5)    # Top border
plt.gca().spines['bottom'].set_linewidth(1.5) # Bottom border
plt.gca().spines['left'].set_linewidth(1.5)   # Left border
plt.gca().spines['right'].set_linewidth(1.5)  # Right border

# Set the thickness of tick marks
plt.tick_params(axis='both', which='both', width=1.5, labelsize=20)

# Plot the data for each dataset
plt.xlim(200, 204)  # Set x-axis limits
plt.ylim(-0.38, -0.44)   # Set y-axis limits

# Plot label
plt.xlabel(r'$\Omega t$')
plt.ylabel(r'$f(\varrho_\mathrm{S}, \varrho^{(2)}_\mathrm{FG})$')

for i, (t_values, x_values) in enumerate(data_sets):
    color = custom_colors[i]  # Get the custom color for this dataset
    window_size = 20000  # You can adjust the window size
    moving_avg = np.convolve(x_values, np.ones(window_size)/window_size, mode='valid')
    t_values_avg = t_values[window_size-1:]  # Adjust t_values for the moving average
    analytic_signal = hilbert(x_values)
    envelope = np.abs(analytic_signal)
    # Plot the data with a label for the legend only for the first three datasets
    
    plt.plot(t_values, x_values, linewidth=1.5, color=color)
    #plt.loglog(t_values_avg, moving_avg, linewidth=1.5, color=color)


plt.tight_layout()

# Add a legend inside the plot for the specified labels
#plt.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98), fancybox=True, ncol=1)

# Save the figure as a PDF file
plt.savefig('mag_xyz_high.pdf', format='pdf')
#plt.savefig('fid_omega_ma.pdf',format='pdf')

# Display the plot
plt.show()


