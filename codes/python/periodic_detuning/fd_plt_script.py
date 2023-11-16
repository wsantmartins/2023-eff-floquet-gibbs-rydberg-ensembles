import matplotlib.pyplot as plt
import numpy as np

# Specify the paths to your data files with double backslashes
data_paths = [
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\f_N=2_w=6_w0=100_Nf=2',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\f_N=3_w=6_w0=100_Nf=2',
    'C:\\Users\\wsant\\OneDrive\\Dokumente\\2023ryd_eng_repo\\codes\\python\\periodic_detuning\\data\\f_N=4_w=6_w0=100_Nf=2'
]

# Define the custom color

#colors without transparency
#custom_colors = ['#d9d9d9','#969696','#525252','#000000']
custom_colors = ['#000000','#67001f','#053061']

#colors with transparency
#custom_colors = ['#92c5de', '#c2a5cf', '#f4a582']

# Define the custom font family
custom_font = {'font.family': 'STIXGeneral'}

# Specify LaTeX preamble to set font family for LaTeX-rendered text
plt.rcParams['text.latex.preamble'] = r'\usepackage{stixgeneral}'  # Use STIX font family in LaTeX

# Apply the custom font family to the Matplotlib configuration
plt.rcParams.update(custom_font)

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
plt.tick_params(axis='both', which='both', width=1.5, labelsize=24)

# Plot the data for each dataset
plt.xlim(0, 100)  # Set x-axis limits
plt.ylim(0, 1)   # Set y-axis limits

#plt.yscale('log')

for i, (t_values, x_values) in enumerate(data_sets):
    #label = f'Dataset {i+1}'  # Label for the legend
    color = custom_colors[i]  # Get the custom color for this dataset
    plt.plot(t_values, x_values, linewidth=1.5, color = color) 

# Save the figure as a PDF file
plt.savefig('fidelity.pdf', format='pdf')

# Display the plot
plt.show()

