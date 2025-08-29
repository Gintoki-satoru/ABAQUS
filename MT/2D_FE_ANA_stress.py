import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import tkinter as tk
from tkinter import filedialog
from scipy.interpolate import interp1d
import numpy as np

def plotInnerAndCalculateError():
    ylabelNorm = [ r'$\sigma_{xx}$ in MPa', r'$\sigma_{yy}$ in MPa',r'$\sigma_{zz}$ in MPa',r'$\tau_{xy}$ in MPa']
    outputStressComp = ['Sigxx','Sigyy','Sigzz','Tauxy']
    outputStressComp2 = ['Sigxx','Sigyy','Sigzz','Tauxy']
    fontsize = 20

    # Initialize Tkinter root and hide the main window
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    # Request the user to select multiple ANA directories
    ana_directories = filedialog.askdirectory(title="Select First ANA Directory")
    ana_directories_list = []
    while ana_directories:
        ana_directories_list.append(ana_directories)
        ana_directories = filedialog.askdirectory(title="Select Next ANA Directory (Cancel to Finish)")
    
    if not ana_directories_list:
        print("No ANA directories selected. Exiting.")
        return

    # Request the user to select corresponding FEM directories
    fem_directories_list = []
    for ana_dir in ana_directories_list:
        fem_directory = filedialog.askdirectory(
            title=f"Select FEM Directory Corresponding to {os.path.basename(ana_dir)}"
        )
        if not fem_directory:
            print(f"No FEM directory selected for {os.path.basename(ana_dir)}. Exiting.")
            return
        fem_directories_list.append(fem_directory)

    # Ensure that each ANA directory has a corresponding FEM directory
    if len(ana_directories_list) != len(fem_directories_list):
        print("Mismatch between the number of ANA and FEM directories. Exiting.")
        return

    # Matplotlib settings for fonts and LaTeX rendering
    plt.rcParams.update({
        'font.family': 'Arial',
        'mathtext.default': 'it',
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',
    })
    plt.rcParams.update({
        'font.size': fontsize,
        'axes.titlesize': fontsize + 2,
        'axes.labelsize': fontsize + 2,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'legend.fontsize': fontsize + 2,
    })

    # Case names and numbers
    case_names = {1: 'Present'}#0: 'CLPT', 
    case_numbers = [1]

    # Loop over each pair of ANA and FEM directories
    for ana_directory, fem_directory in zip(ana_directories_list, fem_directories_list):
        error_data = []  # Store errors for each component and case
        pdfFiles = []  # Store paths to saved plot files

        # Loop over each stress component
        for ii in range(len(ylabelNorm)):
            fig, axs = plt.subplots(1, 2, figsize=(13.33, 7.5), sharex=True, sharey=False)

            # Process each case
            for i, case in enumerate(case_numbers):
                # Read analytical and FEM data
                ana_file = os.path.join(ana_directory, f"{outputStressComp2[ii]}{case}_ANA.csv")
                fem_file = os.path.join(fem_directory, f"{outputStressComp[ii]}_FEM.csv")
                ana_data = pd.read_csv(ana_file)
                fem_data = pd.read_csv(fem_file)

                x_ana = ana_data['x']
                y_ana = ana_data[outputStressComp2[ii]]
                x_fem = fem_data['x']
                y_fem = fem_data[outputStressComp[ii]]

                # Interpolate ANA data to FEM x-axis if point counts differ
                if len(x_ana) != len(x_fem):
                    ana_interp_func = interp1d(x_ana, y_ana, kind='linear', fill_value="extrapolate")
                    y_ana_interp = ana_interp_func(x_fem)
                else:
                    y_ana_interp = y_ana

                # Calculate relative errors
                relative_errors = np.where(y_fem != 0, np.abs((y_fem - y_ana_interp) / y_fem) * 100, np.nan)
                relative_errors = relative_errors[~np.isnan(relative_errors)]  # Remove NaN values
                median_relative_error = np.median(relative_errors) if len(relative_errors) > 0 else np.nan

                # Append error data
                error_data.append(
                    f"Median Relative Error for {outputStressComp[ii]} (Case {case_names[case]}): {median_relative_error:.4f}%"
                )

                # Plot data
                ax = axs[i]
                ax.plot(x_ana, y_ana, linestyle='-', lw=5, c='#004E73', label=f'{case_names[case]}')
                ax.plot(x_fem, y_fem, linestyle='-', markersize=3, lw=5, c='#FDCA00', alpha=0.8, label='FEM')
                ax.set_xlabel(r'$\bar{z}$', fontsize=fontsize + 4)
                ax.set_ylabel(ylabelNorm[ii], fontsize=fontsize + 4, labelpad=10)
                ax.grid(color='grey', linestyle='solid', linewidth=0.5, alpha=0.2)
                ax.set_xlim([-0.49, 0.51])
                ax.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                ax.tick_params(axis='both', which='major', direction='in', width=1, bottom=True, top=True, left=True, right=True)
                ax.minorticks_on()
                ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
                ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
                y_min = round(1.1 * min(y_ana.min(), y_fem.min()), 12)
                y_max = round(1.1 * max(y_ana.max(), y_fem.max()), 12)
                ax.set_ylim(y_min, y_max)
                ax.legend(loc='best', fontsize=fontsize + 4)

            # Adjust layout and save plots
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            file_prefix = 'verification_stress_'
            pdf_file_name = os.path.join(ana_directory, f"{file_prefix}{outputStressComp[ii]}_FEM_ANA.pdf")
            pdfFiles.append(pdf_file_name)
            svg_file_name = os.path.join(ana_directory, f"{file_prefix}{outputStressComp[ii]}_FEM_ANA.svg")
            plt.savefig(pdf_file_name, format='pdf', bbox_inches='tight')
            plt.savefig(svg_file_name, format='svg', bbox_inches='tight')
            plt.close()

        # Save errors to a text file in the ANA directory
        error_file_name = os.path.join(ana_directory, "median_relative_errors_stress.txt")
        with open(error_file_name, "w") as file:
            file.write("\n".join(error_data))

        print(f"Saved plots in {ana_directory}: {pdfFiles}")
        print(f"Median relative errors saved to {error_file_name}")

# Call the function to generate plots and calculate errors
plotInnerAndCalculateError()
