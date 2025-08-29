import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import csv
import tkinter as tk
from tkinter import filedialog as fd

fontsize = 20

# Global variables to store selected directories
selected_directory = None
save_directory = None

def folderDataset(initDir, heading):
    root = tk.Tk()
    root.withdraw()
    wd = fd.askdirectory(title=heading, initialdir=initDir)
    return wd

def read_stress_data(file_prefix, case_number, layer_number):
    global selected_directory
    if selected_directory is None:
        selected_directory = folderDataset(os.getcwd(), "Select the directory containing stress data files")
    
    filename = os.path.join(selected_directory, f"{file_prefix}{case_number}_{layer_number}")
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        return None, None
    
    data = np.genfromtxt(filename, delimiter="\t", invalid_raise=False, filling_values=np.nan)
    x = data[:, 0]
    stress_values = data[:, 1]
    #print(f"[DEBUG] Read file: {filename}, x: {x}, stress_values: {stress_values}")
    
    valid_indices = ~np.isnan(x) & ~np.isnan(stress_values)
    x = x[valid_indices]
    stress_values = stress_values[valid_indices]
    
    #print(f"[DEBUG] After removing NaN values: x: {x}, stress_values: {stress_values}")
    
    return x, stress_values

def save_plot_data_to_csv(x_data, stress_data, file_prefix, case_number):
    global save_directory
    if save_directory is None:
        save_directory = selected_directory#folderDataset(os.getcwd(), "Select the directory to save output files")
    
    stress_component = file_prefix.split('AnaInner')[0]
    csv_filename = os.path.join(save_directory, f"{stress_component}{case_number}_ANA.csv")
    
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['x', stress_component])
        writer.writerows(zip(x_data, stress_data))
    
    #print(f"[DEBUG] Saved data to {csv_filename}")

def plot_inner():
    global save_directory
    
    # Explicitly set fonts and font properties for both math and text
    plt.rcParams.update({
        'font.family': 'Arial',           # Set global font to Arial
        'mathtext.default': 'it',    # Use regular mathtext instead of LaTeX
        'mathtext.fontset': 'custom',     # Custom fontset for math
        'mathtext.rm': 'Arial',           # Use Arial for regular math
        'mathtext.it': 'Arial:italic',    # Use italic Arial for italic math
        'mathtext.bf': 'Arial:bold',      # Use bold Arial for bold math
    })
    # Set global font size for all plot elements
    plt.rcParams.update({
        'font.size': fontsize,               # Set the default font size
        'axes.titlesize': fontsize + 2,      # Title font size
        'axes.labelsize': fontsize + 2,      # Axis labels font size
        'xtick.labelsize': fontsize,         # X-axis tick labels font size
        'ytick.labelsize': fontsize,         # Y-axis tick labels font size
        'legend.fontsize': fontsize + 2,     # Legend font size
    })
    ylabelNorm = [
        r'$\sigma_{xx}$ in MPa', 
        r'$\sigma_{yy}$ in MPa',r'$\sigma_{zz}$ in MPa', r'$\tau_{xy}$ in MPa', r'$\epsilon_{zz}$', r'$\epsilon_{xx}$']#, r'$U_zz$ in mm']
    
    outputStressComp = ['Sigxx','Sigyy','Sigzz','Tauxy','Epsxx','Epsyy']#,'Uxx']
    
    case_numbers = [1]  #0,
    model_names = ['Present']#'CLPT', 
    N_layers = 4
    file_prefixes = ['SigxxAnaInner','SigyyAnaInner','SigzzAnaInner','TauxyAnaInner','EpsxxAnaInner','EpsyyAnaInner']#,'UxxAnaInner']
    
    line_styles = ['-', '--']  
    markers = ['o', 's']  
    colors = ['#004E73','#AFCC50']  
        
    for ii in range(len(ylabelNorm)):
        plt.figure(figsize=(13.33, 5))

        for idx, (case_number, model_name) in enumerate(zip(case_numbers, model_names)):
            x_combined = []
            stress_combined = []

            for layer_number in range(1, N_layers + 1):
                x, stress_values = read_stress_data(file_prefixes[ii], case_number, layer_number)

                if x is not None and stress_values is not None:
                    x_min, x_max = np.min(x), np.max(x)
                    x_range = x_max - x_min
                    #print(f"[DEBUG] Layer {layer_number}, Case {case_number}: x_min = {x_min}, x_max = {x_max}, x_range = {x_range}")

                    x_start = -0.5 + (layer_number - 1) / N_layers
                    x_end = x_start + 1/N_layers
                    #print(f"[DEBUG] Layer {layer_number}, Case {case_number}: x_start = {x_start}, x_end = {x_end}")

                    x_scaled = ((x - x_min) / x_range) * (x_end - x_start) + x_start
                    #print(f"[DEBUG] Layer {layer_number}, Case {case_number}: x_scaled = {x_scaled}")
                    x_combined.extend(x_scaled)
                    stress_combined.extend(stress_values)

                    plt.plot(x_scaled, stress_values, 
                             linestyle='-',#line_styles[idx], 
                             marker=markers[idx], 
                             markersize=5, 
                             color=colors[idx],
                             alpha = 0.8,
                             lw=5,  
                             label=model_name if layer_number == 1 else None)

            save_plot_data_to_csv(x_combined, stress_combined, outputStressComp[ii], case_number)

        plt.grid(color='grey', linestyle='solid', linewidth=0.5, alpha=0.2)
        maxXValue = 0.51
        minYValue, maxYValue = plt.ylim()
        plt.xlim([-maxXValue, maxXValue])
        plt.xticks([-0.5, -0.25, 0, 0.25, 0.5])
        plt.locator_params(axis='y', nbins=6)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.tick_params(axis='both', which='major', color='k', direction='in', width=1,
                        bottom=True, top=True, left=True, right=True)
        plt.tick_params(axis='both', which='minor', color='k', direction='in',
                        bottom=True, top=True, left=True, right=True)
        plt.minorticks_on()
        plt.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        plt.gca().yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles, labels, loc='upper right',fontsize=fontsize+2)

        plt.xlabel(r'$\bar{z}$', fontsize=fontsize+4)
        plt.ylabel(ylabelNorm[ii], fontsize=fontsize+4)
        plt.tight_layout(pad=2.0)
        # Save as both .pdf and .svg formats
        pdf_filename = os.path.join(save_directory, f"{outputStressComp[ii]}.pdf")
        svg_filename = os.path.join(save_directory, f"{outputStressComp[ii]}.svg")
        
        plt.savefig(pdf_filename, format='pdf', bbox_inches='tight', dpi=400)
        plt.savefig(svg_filename, format='svg', bbox_inches='tight', dpi=400)
        plt.close()

# Call the plot function to execute
plot_inner()

