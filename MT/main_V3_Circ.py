import datetime
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PyPDF2 import PdfFileMerger
import tkinter as tk
from tkinter import filedialog as fd

# Global variables
pdfFiles = []

# Tkinter-Window:
win = tk.Tk()

# Working directory:
mainFolder = 'C:\\00_Promotion\\21_Auswertung_Verifikation\\00_Post_Processing_Python_Skript'

# Define directory of the data:
def folderDataset(initDir, heading):
    wd = fd.askdirectory(title=heading, initialdir=initDir)
    return wd


# Data of finite element method:
initDirFE = 'C:\\00_Promotion\\12_Numerik'
wd_fe = folderDataset(initDirFE, 'Dataset of the finite element computations')

# Data of semi-analytical method:
initDirSA = 'C:\\00_Promotion\\11_Analytik'
wd_sa = folderDataset(initDirSA, 'Dataset of the semi-analytical analyses')

# Close Tkinter-Window:
win.destroy()

def importDatasetFE(nameDataset):
    file = [wd_fe, '/', nameDataset, '.txt']
    data = pd.read_csv(''.join(file), header=None,  sep=" ")
    return data


def importDatasetSA(nameDataset):
    file = [wd_sa, '/', nameDataset]
    data = pd.read_csv(''.join(file), header=None,  sep="\t")
    return data


# Import layup of the analysed composite:
compositeLayup = importDatasetSA('compositeLayup')
compositeLayupTitle = r'$\left[ \pm 45^\circ \right]_S$'
N = len(compositeLayup.columns)

# Import geometry of the analysed composite:
# Meaning:  1 - Mathematical layers P
#           2 - Order Psi of displacement-based approach
#           3 - Thickness ratio R/h
#           4 - Thickness dL of each physical layer
#           5 - Length ratio l/h
#           6 - Beginning angle of the composite laminate
#           7 - Opening angle of the composite laminate
#           8 - Evaluation angle of the composite laminate
compositeGeometry = importDatasetSA('compositeGeometry')

def normRadii(dataset):
    modDataset = dataset
    r_i = dataset[0,0]
    r_n = dataset[-1,0]
    for ii in range(0,np.shape(dataset)[1],2):
        for jj in range(np.shape(dataset)[0]):
            modDataset[jj,ii] = (dataset[jj,ii]-((r_i+r_n)/2))/(r_n-r_i)
    return modDataset
             

def importStressInnerAnalytical():
    stressComponent = ['Sigrr', 'Sigtt', 'Sigzz', 'Tautz', 'Taurz', 'Taurt']
    sumStresses = []
    sumStressesNormRadii = []
    numberOfEvalPoints = 20
    for ii in range(len(stressComponent)):
        modData = []
        for jj in range(N):
            file = [stressComponent[ii], 'AnaInner', str(jj+1)]
            dataset = importDatasetSA(''.join(file))
            ll = 0
            step = round(len(dataset.index)/numberOfEvalPoints)
            for kk in range(0,len(dataset.index),step):
                modData.append([dataset[0][kk],dataset[1][kk]])
                ll += 1
            modData.append([dataset[0].iloc[-1],dataset[1].iloc[-1]])
        if len(sumStresses) == 0:
            sumStresses = modData
        else:
            sumStresses = np.hstack((sumStresses,modData))
    sumStressesNormRadii =  normRadii(sumStresses)
    return sumStressesNormRadii


def importStressInnerSemiAnalytical():
    stressComponent = ['Sigrr', 'Sigtt', 'Sigzz','Tautz','Taurz','Taurt']
    sumStresses = []
    sumStressesNormRadii = []
    numberOfEvalPoints = 20
    M = compositeGeometry[0][0]
    for ii in range(len(stressComponent)):
        modData = []
        for jj in range(N*M):
            file = [stressComponent[ii], 'SAInner', str(jj+1)]
            dataset = importDatasetSA(''.join(file))
            ll = 0
            step = round(len(dataset.index)*M/numberOfEvalPoints)
            for kk in range(0,len(dataset.index),step):
                modData.append([dataset[0][kk],dataset[1][kk]])
                ll += 1
            modData.append([dataset[0].iloc[-1],dataset[1].iloc[-1]])
        if len(sumStresses) == 0:
            sumStresses = modData
        else:
            sumStresses = np.hstack((sumStresses,modData))
    sumStressesNormRadii =  normRadii(sumStresses)
    return sumStressesNormRadii


def importStressInnerFiniteElement():
    # [Sigrr, Sigtt, Sigzz, Taurt, Taurz, Tautz]
    stressComponent = [1,2,3,4,5,6]
    sumStresses = []
    sumStressesNormRadii = []
    for ii in range(N):
        dataset = []
        radii = []
        noStressComponents = 6
        circPos = ['TS','TG']
        axialPos = ['ZS','ZG']
        noMeanValue = len(circPos) + len(axialPos)
        noMeanValueRunVar = 0
        for jj in range(len(circPos)):
            for kk in range(len(axialPos)):
                file = ['FEInner', str(ii+1), '_', circPos[jj], '_', axialPos[kk]]
                dataset.append(importDatasetFE(''.join(file)))
                if len(radii) == 0:
                    radii.append(list(dataset[noMeanValueRunVar][9]))
                else:
                    pass
                noMeanValueRunVar += 1
        stressMean = np.zeros((len(dataset[0].index),noStressComponents))
        modData = np.zeros((len(dataset[0].index),2*len(stressComponent)))
        for mm in range(len(dataset[0].index)):
            for nn in range(noStressComponents):
                sum = 0.0
                for jk in range(noMeanValue):
                        sum += dataset[jk][nn][mm]
                stressMean[mm][nn] = sum/(noMeanValue)
        for mm in range(len(stressComponent)):
            modData[:,2*mm] = radii[0]
            modData[:,2*mm+1] = stressMean[:,stressComponent[mm]-1]
        sumStresses.append(modData)
    sumStressesNormRadii =  normRadii(np.concatenate(sumStresses,axis=0))
    return sumStressesNormRadii


linestyle = [':','-']
linewidth = [2.0, 1.0]
colour = ['k', 'k']
legend = ['FE', 'SA']
fontsize = 12

def plotInner(*argv):
    linestyle = [':', '-', '-']
    linewidth = [2.0, 1.0, 1.0]
    colour = ['k', 'k', 'r']
    legend = ['FE', 'SA', 'A']
    ylabel = [r'$\bar{\sigma}_{rr}$',r'$\bar{\sigma}_{\theta \theta}$',r'$\bar{\sigma}_{zz}$', 
              r'$\bar{\tau}_{\theta z}$', r'$\bar{\tau}_{r z}$', r'$\bar{\tau}_{r \theta}$']
    outputStressComp = ['Sigrr', 'Sigtt', 'Sigzz', 'Tautz', 'Taurz', 'Taurt']
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    for ii in range(len(ylabel[0:6])):
        for jj in range(len(argv)):
            plt.plot(argv[jj][:,2*ii], argv[jj][:,2*ii+1], ls = linestyle[jj], lw = linewidth[jj],
                     c = colour[jj], label = legend[jj])
        plt.legend(loc='upper right')
        plt.xlabel(r'$\bar{r}$',fontsize=fontsize)
        plt.xticks([-0.5, -0.25, 0, 0.25, 0.5],fontsize=fontsize)
        plt.ylabel(ylabel[ii],fontsize=fontsize)
        plt.locator_params(axis='y',nbins=7)
        plt.tick_params(axis='both', labelsize=fontsize)
        fileName = [outputStressComp[ii], 'InnerRadial.pdf']
        pdfFiles.append(''.join(fileName))
        plt.savefig(''.join(fileName))
        plt.show()
    return
       

def normLengthSA(dataset):
    modDataset = dataset
    z_0 = dataset[0,0]
    z_e = dataset[-1,0]
    for ii in range(0,np.shape(dataset)[1],2):
        for jj in range(np.shape(dataset)[0]):
            modDataset[jj,ii] = dataset[jj,ii]/(z_e-z_0)
    return modDataset


def normLengthFE(dataset):
    modDataset = dataset
    z_0 = dataset[-1,0]
    for ii in range(0,np.shape(dataset)[1],2):
        for jj in range(np.shape(dataset)[0]):
            modDataset[jj,ii] = (dataset[jj,ii]-(z_0/2))/z_0
    return modDataset


def importStressSemiAnalyticalFreeEdgeInterlaminar(interface):
    stressComponent = ['Sigrr', 'Taurz', 'Taurt']
    sumStresses = []
    sumStressesNormLength = []
    for ii in range(len(stressComponent)):
        file = [stressComponent[ii], 'SAInterfaceCirc_', str(interface)]
        dataset = importDatasetSA(''.join(file))
        if len(sumStresses) == 0:
            sumStresses = dataset
        else:
            sumStresses = np.hstack((sumStresses,dataset))
    sumStressesNormLength =  normLengthSA(sumStresses)
    return sumStressesNormLength


def importStressSemiAnalyticalFreeEdgeIntralaminar(layer):
    stressComponent = ['Sigtt', 'Sigzz', 'Tautz']
    sumStresses = []
    sumStressesNormLength = []
    for ii in range(len(stressComponent)):
        file = [stressComponent[ii], 'SALayerCirc', str(layer)]
        dataset = importDatasetSA(''.join(file))
        if len(sumStresses) == 0:
            sumStresses = dataset
        else:
            sumStresses = np.hstack((sumStresses,dataset))
    sumStressesNormLength =  normLengthSA(sumStresses)
    return sumStressesNormLength


def importStressFiniteElementFreeEdgeInterlaminar(interface):
    # [Sigrr, Sigtt, Sigzz, Tautz, Taurz, Taurt]
    stressComponent = [1,5,6]
    sumStresses = []
    sumStressesNormLength = []
    dataset = []
    axial = []
    noStressComponents = 6
    noMeanValue = 2
    layer = [compositeLayup[interface-1][0], compositeLayup[interface][0]]
    anglePos = ['ZS', 'ZG']
    for kk in range(len(anglePos)):
        for ll in range(len(layer)):
            if ll == 0:
                for mm in range(noMeanValue):
                    file = ['FEInterfaceCirc', str(interface), '_Layer_',
                            str(layer[ll]),  '_RS_', str(anglePos[kk]), '_', str(mm+1)]
                    dataset.append(importDatasetFE(''.join(file)))
                    axial.append(list(dataset[mm][10]))
            elif ll == 1:
                for mm in range(noMeanValue):
                    file = ['FEInterfaceCirc', str(interface), '_Layer_',
                            str(layer[ll]),  '_RG_', str(anglePos[kk]), '_', str(mm+1)]
                    dataset.append(importDatasetFE(''.join(file)))
            else:
                pass
    stressMean = np.zeros((len(dataset[0].index),noStressComponents))
    modData = np.zeros((len(dataset[0].index),2*len(stressComponent)))
    for nn in range(len(dataset[0].index)):
        for kk in range(noStressComponents):
            sum = 0.0
            for pp in range(len(layer)*len(anglePos)*noMeanValue):
                sum += dataset[pp][kk][nn]
            stressMean[nn][kk] = sum/(len(layer)*len(anglePos)*noMeanValue)
    for mm in range(len(stressComponent)):
        modData[:,2*mm] = axial[0]
        modData[:,2*mm+1] = stressMean[:,stressComponent[mm]-1]
    sumStresses.append(modData)
    sumStressesNormLength =  normLengthFE(np.concatenate(sumStresses,axis=0))
    return sumStressesNormLength


def importStressFiniteElementFreeEdgeIntralaminar(layer):
    # [Sigrr, Sigtt, Sigzz, Tautz, Taurz, Taurt]
    stressComponent = [2,3,4]
    sumStresses = []
    sumStressesNormLength = []
    dataset = []
    axial = []
    noStressComponents = 6
    noMeanValue = 2
    anglePos = ['ZS', 'ZG']
    for ll in range(len(anglePos)):
        if ll == 0:
            for mm in range(noMeanValue):
                file = ['FELayerCirc', str(layer), '_RS_', anglePos[ll]]
                dataset.append(importDatasetFE(''.join(file)))
                axial.append(list(dataset[mm][10]))
        elif ll == 1:
            for mm in range(noMeanValue):
                file = ['FELayerCirc', str(layer), '_RG_', anglePos[ll]]
                dataset.append(importDatasetFE(''.join(file)))
        else:
            pass
    stressMean = np.zeros((len(dataset[0].index),noStressComponents))
    modData = np.zeros((len(dataset[0].index),2*len(stressComponent)))
    for nn in range(len(dataset[0].index)):
        for kk in range(noStressComponents):
            sum = 0.0
            for pp in range(len(anglePos)*noMeanValue):
                sum += dataset[pp][kk][nn]
            stressMean[nn][kk] = sum/(len(anglePos)*noMeanValue)
    for mm in range(len(stressComponent)):
        modData[:,2*mm] = axial[0]
        modData[:,2*mm+1] = stressMean[:,stressComponent[mm]-1]
    sumStresses.append(modData)
    sumStressesNormLength =  normLengthFE(np.concatenate(sumStresses,axis=0))
    return sumStressesNormLength


def plotFreeEdgeInterlaminar(interface,*argv):
    titleName = ''.join([compositeLayupTitle, r' - Interface ', str(interface)])
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ylabel = [r'$\sigma_{rr}$', r'$\tau_{r z}$', r'$\tau_{r \theta}$']
    outputStressComp = ['Sigrr', 'Taurz', 'Taurt']
    for ii in range(len(ylabel)):
        for jj in range(len(argv)):
            plt.plot(argv[jj][:,2*ii], argv[jj][:,2*ii+1], ls = linestyle[jj], lw = linewidth[jj],
                     c = colour[jj], label = legend[jj])
        plt.legend(loc='upper center')
        plt.xlabel(r'$\bar{z}$')
        plt.xticks([-0.5, -0.25, 0, 0.25, 0.5])
        plt.ylabel(ylabel[ii])
        plt.title(titleName)
        fileName = [outputStressComp[ii], 'StressConcDecayInterface', str(interface), '.pdf']
        pdfFiles.append(''.join(fileName))
        plt.savefig(''.join(fileName))
        plt.show()
    return


def plotFreeEdgeIntralaminar(layer,*argv):
    titleName = ''.join([compositeLayupTitle, r' - Layer ', str(layer)])
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ylabel = [r'$\sigma_{\theta \theta}$', r'$\sigma_{z z}$', r'$\tau_{\theta z}$']
    outputStressComp = ['Sigtt', 'Sigzz', 'Tautz']
    for ii in range(len(ylabel)):
        for jj in range(len(argv)):
            plt.plot(argv[jj][:,2*ii], argv[jj][:,2*ii+1], ls = linestyle[jj], lw = linewidth[jj],
                     c = colour[jj], label = legend[jj])
        plt.legend(loc='upper center')
        plt.xlabel(r'$\bar{z}$')
        plt.xticks([-0.5, -0.25, 0, 0.25, 0.5])
        plt.ylabel(ylabel[ii])
        plt.title(titleName)
        fileName = [outputStressComp[ii], 'StressConcDecayLayer', str(layer), '.pdf']
        pdfFiles.append(''.join(fileName))
        plt.savefig(''.join(fileName))
        plt.show()
    return


def plotFreeEdgeInterlaminarMultipleInterfaces(interface):
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ylabel = [r'$\bar{\sigma}_{rr}$',r'$\bar{\tau}_{r z}$', r'$\bar{\tau}_{r \theta}$']
    outputStressComp = ['Sigrr', 'Taurz', 'Taurt']
    for ii in range(len(ylabel)):
        for jj in range(len(interface)):
            stress = [importStressFiniteElementFreeEdgeInterlaminar(jj+1),
                      importStressSemiAnalyticalFreeEdgeInterlaminar(jj+1)]
            for kk in range(len(stress)):
                if jj == 0:
                    plt.plot(stress[kk][:,2*ii], stress[kk][:,2*ii+1]*4/10, ls = linestyle[kk], lw = linewidth[kk],
                             c = colour[kk], label = legend[kk])
                else:
                    plt.plot(stress[kk][:,2*ii], stress[kk][:,2*ii+1]*4/10, ls = linestyle[kk], lw = linewidth[kk],
                                 c = colour[kk])
        plt.legend(loc='upper center')
        plt.xlabel(r'$\bar{z}$')
        plt.xticks([-0.5, -0.25, 0, 0.25, 0.5])
        plt.locator_params(axis='y', nbins=7)
        plt.ylabel(''.join([ylabel[ii]]))
        fileName = [outputStressComp[ii], 'StressConcDecayMultipleInterface.pdf']
        pdfFiles.append(''.join(fileName))
        plt.savefig(''.join(fileName))
        plt.show()
    return

def plotFreeEdgeIntralaminarMultipleLayers(layer):
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ylabel = [r'$\bar{\sigma}_{\theta \theta}$', r'$\bar{\sigma}_{z z}$', r'$\bar{\tau}_{\theta z}$']
    outputStressComp = ['Sigtt', 'Sigzz', 'Tautz']
    for ii in range(len(ylabel)):
        for jj in range(len(layer)):
            stress = [importStressFiniteElementFreeEdgeIntralaminar(jj+1),
                      importStressSemiAnalyticalFreeEdgeIntralaminar(jj+1)]
            for kk in range(len(stress)):
                if jj == 0:
                    plt.plot(stress[kk][:,2*ii], stress[kk][:,2*ii+1]*4/10, ls = linestyle[kk], lw = linewidth[kk],
                             c = colour[kk], label = legend[kk])
                else:
                    plt.plot(stress[kk][:,2*ii], stress[kk][:,2*ii+1]*4/10, ls = linestyle[kk], lw = linewidth[kk],
                                 c = colour[kk])
        plt.legend(loc='upper center')
        plt.xlabel(r'$\bar{z}$')
        plt.xticks([-0.5, -0.25, 0, 0.25, 0.5])
        plt.locator_params(axis='y', nbins=7)
        plt.ylabel(ylabel[ii])
        fileName = [outputStressComp[ii], 'StressConcDecayLayer', str(layer), '.pdf']
        pdfFiles.append(''.join(fileName))
        plt.savefig(''.join(fileName))
        plt.show()
    return


materialParameter = importDatasetSA('materialParameters')

# Create *.pdf-File with material configuration and composite geometry:
def createPdfMatConfigCompGeom(*argv):
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ylabel = [r'$\sigma_{rr}$']
    for ii in range(len(ylabel)):
        for jj in range(len(argv)):
            plt.plot(argv[jj][:,2*ii], argv[jj][:,2*ii+1], ls = linestyle[jj], lw = linewidth[jj],
                     c = colour[jj], label = legend[jj], alpha=0.1)
        plt.xlabel(r'$\bar{z}$')
        plt.xticks([-0.5, -0.25, 0, 0.25, 0.5])
        plt.ylabel(ylabel[ii])
        pltCenter = (max(argv[jj][:,2*ii+1])+min(argv[jj][:,2*ii+1]))/2
        configText = ''.join([r'Composite Layup: ', compositeLayupTitle, ',\n',
                             r' Math. layers: $P$ = ', str(compositeGeometry[0][0]), r', Order: $\Psi$ = ', str(compositeGeometry[1][0]), ',\n',
                             r' Thickness ratio: $R/h$ = ', str(compositeGeometry[2][0]), r', Length ratio: $L/h$ = ', str(compositeGeometry[4][0]), ',\n',
                             r' Layer thickness: $d_L$ = ', str(compositeGeometry[3][0]), r', Opening angle: $\theta_0$ = ', str(compositeGeometry[5][0]), '\n', '\n',
                             r' Load:', '\n', r' $\bar{M} = 10$ N, $\bar{F} \left( \alpha \right) = 0$ N/mm,  $\Delta T \left( r \right) = 0$ K' , '\n', '\n',
                             r' Material parameters:', '\n',
                             r' $E_1$ = ', str(materialParameter[1][0]), r' GPa,  $E_2$ = $E_3$ = ', str(materialParameter[0][0]), r' GPa,' , '\n',
                             r' $G_{23}$ = ', str(materialParameter[1][1]), r' GPa,  $G_{12}$ = $G_{13}$ = ', str(materialParameter[0][1]), r' GPa,' , '\n', 
                             r' $\nu_{12}$ = $\nu_{13}$ = $\nu_{23}$ =', str(materialParameter[0][2]), ', \n'
                             #r' $\alpha_{1}$ = ', str(format(materialParameter[1][3],'.2e')), r' 1/°C,   $\alpha_{2}$ = $\alpha_{3}$ = ', str(format(materialParameter[0][3],'.2e')), ' 1/°C'
                             ])
        plt.text(0.0, pltCenter, configText, fontsize = 11,
                 bbox=dict(facecolor='none',boxstyle='square',edgecolor='none'),
                 ha='center', va='center')
        pdfFiles.append(''.join(['MatConfigAndCompGeom.pdf']))
        plt.savefig(''.join(['MatConfigAndCompGeom.pdf']))
        plt.show()
    return


def mergeDeletePdfFiles():
    timeEval = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    mergeInstance = PdfFileMerger()
    for ii in pdfFiles:
        mergeInstance.append(ii)
    mergeInstance.write(''.join([str(timeEval),'_Summary.pdf']))
    mergeInstance.close()
    for jj in pdfFiles:
        file = [mainFolder, '\\', jj]
        if os.path.exists(''.join(file)):
            os.remove(''.join(file))
            

#createPdfMatConfigCompGeom(importStressSemiAnalyticalFreeEdgeInterlaminar(1),importStressFiniteElementFreeEdgeInterlaminar(1))


plotInner(importStressInnerFiniteElement(),importStressInnerSemiAnalytical(),importStressInnerAnalytical())


#for ii in range(N-1):
#    plotFreeEdgeInterlaminar(ii+1,importStressFiniteElementFreeEdgeInterlaminar(ii+1),importStressSemiAnalyticalFreeEdgeInterlaminar(ii+1))


#for ii in range(N):
#    plotFreeEdgeIntralaminar(ii+1,importStressFiniteElementFreeEdgeIntralaminar(ii+1),importStressSemiAnalyticalFreeEdgeIntralaminar(ii+1))


# plotFreeEdgeInterlaminarMultipleInterfaces([1,2,3])

# plotFreeEdgeIntralaminarMultipleLayers([1,2,3,4])


mergeDeletePdfFiles()
