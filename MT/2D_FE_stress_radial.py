import os
import matplotlib.pyplot as plt
from matplotlib import rcParams, ticker
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog as fd

# Global variables
pdfFiles = []
csvFiles = []

# Tkinter-Window for folder selection
win = tk.Tk()
win.withdraw()  # Hide the main tkinter window

# Function to select multiple folders for datasets
def select_multiple_folders(initDir, heading):
    print(f"\n--- {heading} ---")
    folders = []
    while True:
        folder = fd.askdirectory(title="Select a folder (or cancel to finish)", initialdir=initDir)
        if not folder:
            break
        folders.append(folder)
        print(f"Selected: {folder}")
    return folders

# Function to import dataset
def importDataset(nameDataset, directory):
    file = [directory, '/', nameDataset]
    data = pd.read_csv(''.join(file), header=None, sep="\t")
    return data

# Function to import 2D FEM dataset
def importDataset2DFE(nameDataset, directory):
    file = [directory, '/', nameDataset, '.txt']
    data = pd.read_csv(''.join(file), header=None, sep=" ")
    return data

# Function to normalize radii
def normRadii(dataset):
    modDataset = dataset.copy()
    r_i = dataset[0, 0]
    r_n = dataset[-1, 0]
    for ii in range(0, np.shape(dataset)[1], 2):
        for jj in range(np.shape(dataset)[0]):
            modDataset[jj, ii] = (dataset[jj, ii] - ((r_i + r_n) / 2)) / (r_n - r_i)
    return modDataset

fontsize = 24

# Function to import stress data from 2D FEM
def importStressInner2DFiniteElement(directory, N):
    stressComponent = [1, 2, 3, 4]  # Stress components to be processed
    sumStresses = []

    for ii in range(N):
        dataset = []
        radii = []
        noStressComponents = 4
        circPos = ['TS', 'TG']

        for jj in range(len(circPos)):
            file = ['FEInner', str(ii + 1), '_', circPos[jj]]
            dataset.append(importDataset2DFE(''.join(file), directory))
            if len(radii) == 0:
                radii.append(list(dataset[jj][6]))

        stressMean = np.zeros((len(dataset[0].index), noStressComponents))
        modData = np.zeros((len(dataset[0].index), 2 * len(stressComponent)))

        for mm in range(len(dataset[0].index)):
            for nn in range(noStressComponents):
                sum = 0.0
                for jk in range(len(circPos)):
                    sum += dataset[jk][nn][mm]
                stressMean[mm][nn] = sum / len(circPos)

        for mm in range(len(stressComponent)):
            modData[:, 2 * mm] = radii[0]
            modData[:, 2 * mm + 1] = stressMean[:, stressComponent[mm] - 1]

        sumStresses.append(modData)

    sumStressesNormRadii = normRadii(np.concatenate(sumStresses, axis=0))
    return sumStressesNormRadii

# Function to save data to CSV
def saveToCSV(directory, save_dir, N):
    ylabelNorm = [r'$\sigma_{xx}$ in MPa', r'${\sigma}_{yy}$ in MPa', r'${\sigma}_{zz}$ in MPa', r'${\tau}_{xy}$ in MPa']
    outputStressComp = ['Sigxx', 'Sigyy', 'Sigzz', 'Tauxy']

    for ii in range(len(ylabelNorm)):
        resultsInner2DFE = importStressInner2DFiniteElement(directory, N)
        fileName = [outputStressComp[ii], '_FEM.csv']
        csvFiles.append(''.join(fileName))

        # Save data to CSV
        csv_data = pd.DataFrame(resultsInner2DFE[:, [2 * ii, 2 * ii + 1]], columns=['x', outputStressComp[ii]])
        csv_data.to_csv(os.path.join(save_dir, ''.join(fileName)), index=False)

# Function to plot stress data
def plotInner(directory, save_dir, N):
    ylabelNorm = [r'$\sigma_{xx}$ in MPa', r'${\sigma}_{yy}$ in MPa', r'${\sigma}_{zz}$ in MPa', r'${\tau}_{xy}$ in MPa']
    outputStressComp = ['Sigxx', 'Sigyy', 'Sigzz', 'Tauxy']
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

    for ii in range(len(ylabelNorm)):
        # Start a new figure for each plot
        plt.figure(figsize=(16,9))

        resultsInner2DFE = importStressInner2DFiniteElement(directory, N)
        plt.plot(resultsInner2DFE[:, 2 * ii], resultsInner2DFE[:, 2 * ii + 1], linestyle='-', lw=4,
                 c='#004E73', label=r'2D FEA')

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
        plt.legend(loc='lower center')
        plt.xlabel(r'$\bar{x}$', fontsize=fontsize + 4)
        plt.ylabel(ylabelNorm[ii], fontsize=fontsize + 4, labelpad=10)
        plt.ylim((round(minYValue - 0.5 * abs(minYValue), 4), round(maxYValue + 0.1 * abs(maxYValue), 4)))

        # Save plots with updated file names
        file_prefix = '0_90S_'
        pdf_file_name = f"{file_prefix}{outputStressComp[ii]}_FE.pdf"
        pdfFiles.append(pdf_file_name)

        svg_file_name = f"{file_prefix}{outputStressComp[ii]}_FE.svg"
        plt.savefig(os.path.join(save_dir, pdf_file_name), bbox_inches='tight', format='pdf')
        plt.savefig(os.path.join(save_dir, svg_file_name), bbox_inches='tight', format='svg')

        # Clear the figure after saving
        plt.close()

    print(f"Saved plots: {pdfFiles}")

# Main function to handle multiple folders
if __name__ == "__main__":
    initDir = 'C:\\Users\\Pavit\\HESSENBOX-DA\\austausch\\plot\\L-CFRP\\0-90-90-0'
    folders = select_multiple_folders(initDir, "Select folders for processing")

    for folder in folders:
        print(f"\nProcessing folder: {folder}")
        save_dir = os.path.join(folder, "output")
        os.makedirs(save_dir, exist_ok=True)

        # Load required data
        compositeLayup = importDataset('compositeLayup', folder)
        N = len(compositeLayup.columns)

        plotInner(folder, save_dir, N)
        saveToCSV(folder, save_dir, N)
        print(f"Saved output to {save_dir}")
