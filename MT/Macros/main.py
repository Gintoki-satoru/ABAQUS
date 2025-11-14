''' File name: Script_CylindricallyCurvedComposite_SubModeling.py
	Date created: 2019-08-29
	Date last modified: 2020-05-13 - Andreas Kappel
'''

#---------------------------------------------------------------
# Python-Packages:
import __main__
import os
import time
import datetime as dt
import numpy as np
import operator

# Abaqus-Packages:
from abaqus import *
from abaqusConstants import *
import sketch
import part
import regionToolset
import assembly
import material 
import section
import step
import mesh
import interaction
import load
import job
import visualization
import xyPlot
import odbAccess
import displayGroupOdbToolset as dgot
import displayGroupMdbToolset as dgmt

path_modules = 'D:\\psingh\\MT\\ABAQUS\\MT\\Macros'
os.chdir(path_modules)

# Further packages:
import coordinateTransformation as ct

# Global variables
elementsLoc = [['r < r_eval','r > r_eval'],['t < t_eval','t > t_eval'],['z < z_eval','z > z_eval']]
elementNodesLoc = [['first three element nodes','last three element nodes'],['r < r_eval','r > r_eval'],['t < t_eval','t > t_eval'], ['z < z_eval','z > z_eval']]

#---------------------------------------------------------------
# Ausgabe der extrapolierten Spannungen an den jeweiligen Elementknoten des ersten und letzten Elements des Sets (in Breitenrichtung):
def getDispComp3D():
	return('x','y','z')

def getStressComp3D():
	return('Sigma_rr','Sigma_tt','Sigma_zz','Tau_rt','Tau_rz','Tau_tz')

def getStressesCoordinatesNode(step,FieldValueObject):
	stressComp = {stressCompStr:FieldValueObject.data[ii] for ii,stressCompStr in enumerate(getStressComp3D())}
	dispComp = {dispCompStr:compositeSubModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].\
				instance.getNodeFromLabel(FieldValueObject.nodeLabel).coordinates[ii] \
				for ii,dispCompStr in enumerate(getDispComp3D())}
	rCoord = ct.cart2pol_radius(dispComp['x'],dispComp['y'])
	phiCoord = (180*(ct.cart2pol_theta(dispComp['x'],dispComp['y']))/np.pi)
	
	return(stressComp[getStressComp3D()[0]], stressComp[getStressComp3D()[1]], stressComp[getStressComp3D()[2]],\
		   stressComp[getStressComp3D()[5]], stressComp[getStressComp3D()[4]], stressComp[getStressComp3D()[3]],\
		   dispComp[getDispComp3D()[0]], dispComp[getDispComp3D()[1]], dispComp[getDispComp3D()[2]],rCoord, phiCoord,\
		   FieldValueObject.nodeLabel, FieldValueObject.elementLabel)

# Ausgabe der Spannungen in den Integrationspunkten der jeweiligen Elemente des Sets:
def getStressesIntegrationPoint(FieldValueObject):
	stressComp = {stressCompStr:FieldValueObject.data[ii] for ii,stressCompStr in enumerate(getStressComp3D())}
	
	return (FieldValueObject.elementLabel, FieldValueObject.integrationPoint,\
			stressComp[getStressComp3D()[0]], stressComp[getStressComp3D()[1]], stressComp[getStressComp3D()[2]],\
			stressComp[getStressComp3D()[5]], stressComp[getStressComp3D()[4]], stressComp[getStressComp3D()[3]])

#---------------------------------------------------------------
# Funktionen zur Anwendung von Sortieralgorithmen:
def getCoordinateR(unsortedElementNodes):
	return ct.cart2pol_radius(unsortedElementNodes[1],unsortedElementNodes[2])

def getCoordinateT(unsortedElementNodes):
	return ct.cart2pol_theta(unsortedElementNodes[1],unsortedElementNodes[2])

def getCoordinateZ(unsortedElementNodes):
	return unsortedElementNodes[3]

def getNodeLabel(stressLastFrameElementNodalUnsorted):
	return stressLastFrameElementNodalUnsorted.nodeLabel

def getIntegrationPoint(stressGaussPointsUnsorted):
	return stressGaussPointsUnsorted[2]

#---------------------------------------------------------------
# Funktion zur Sortierung eines Element-Sets:
def sortElements(sortedElementNodes, FieldValueObject):
	FieldValueObjectSort = []
	for ii in range(len(sortedElementNodes)):
		NumberOfElementsSameNodeList = []
		for jj in range(len(FieldValueObject)):
			if FieldValueObject[jj].nodeLabel == sortedElementNodes[ii][0]:
				NumberOfElementsSameNodeList.append(jj)
		if len(NumberOfElementsSameNodeList) == 1:
			FieldValueObjectSort.append(FieldValueObject[NumberOfElementsSameNodeList[0]])
		elif len(NumberOfElementsSameNodeList) == 2:
			if FieldValueObject[NumberOfElementsSameNodeList[0]].elementLabel == FieldValueObjectSort[len(FieldValueObjectSort)-1].elementLabel:
				FieldValueObjectSort.append(FieldValueObject[NumberOfElementsSameNodeList[0]])
				FieldValueObjectSort.append(FieldValueObject[NumberOfElementsSameNodeList[1]])
			elif FieldValueObject[NumberOfElementsSameNodeList[1]].elementLabel == FieldValueObjectSort[len(FieldValueObjectSort)-1].elementLabel:
				FieldValueObjectSort.append(FieldValueObject[NumberOfElementsSameNodeList[1]])
				FieldValueObjectSort.append(FieldValueObject[NumberOfElementsSameNodeList[0]])
	return FieldValueObjectSort

#---------------------------------------------------------------
def getStressesDisplacements(set,step,sortDirection):
	# Einlesen des im Preprocessing definierten Knoten-/Element-Sets:
	compositeSetElements=compositeSubModelOdbObject.rootAssembly.elementSets[set]
	compositeSetElementNodes=compositeSubModelOdbObject.rootAssembly.nodeSets[set]
	
	# Bestimmung des letzten Frames der FE-Analyse:
	indexLastFrame = len(compositeSubModelOdbObject.steps[step].frames) - 1
	
	# Bestimmung der Knoten-Koordinaten/Knoten-Verschiebungen/Element-Spannungen im zylindrischen Koordinatensystem und abspeichern dieser in Arrays:
	coordsFirstFrame = compositeSubModelOdbObject.steps[step].frames[0].fieldOutputs['COORD'].getSubset(region=compositeSetElementNodes).values
	dispLastFrame = compositeSubModelOdbObject.steps[step].frames[indexLastFrame].fieldOutputs['U'].getSubset(region=compositeSetElementNodes).getTransformedField(datumCsys=postProc_Cyl_CS).values
	stressLastFrameElementNodalTemp = compositeSubModelOdbObject.steps[step].frames[-1].fieldOutputs['S'].getSubset(position=ELEMENT_NODAL, region=compositeSetElements).getTransformedField(datumCsys=postProc_Cyl_CS).values
	stressLastFrameGaussPoint = compositeSubModelOdbObject.steps[step].frames[-1].fieldOutputs['S'].getSubset(position=INTEGRATION_POINT, region=compositeSetElements).getTransformedField(datumCsys=postProc_Cyl_CS).values
	
	# Knotenkoordinaten des Knotensets in aufsteigender Reihenfolge sortieren - freier Rand (z = 0) bis entsprechend Ende (Knotenkoordinaten des letzten Elementknotens):
	unsortedElementNodes = []
	for ii in range(len(dispLastFrame)):
		unsortedElementNodes.append((dispLastFrame[ii].nodeLabel,coordsFirstFrame[ii].data[0],coordsFirstFrame[ii].data[1],coordsFirstFrame[ii].data[2]))
	
	if sortDirection == 'r':
		sortedElementNodes = sorted(unsortedElementNodes, key=getCoordinateR)
	elif sortDirection == 't':
		sortedElementNodes = sorted(unsortedElementNodes, key=getCoordinateT)
	elif sortDirection == 'z':
		sortedElementNodes = sorted(unsortedElementNodes, key=getCoordinateZ)
	elif sortDirection == 'tz':
		sortedElementNodes = sorted(unsortedElementNodes, key=operator.itemgetter(2,3))
	
	return(compositeSetElements,stressLastFrameElementNodalTemp,stressLastFrameGaussPoint,sortedElementNodes)

def elemLocCol(colDir):
	if colDir == 'r':
		elements = {elementsLocT:{elementsLocZ:[] for elementsLocZ in elementsLoc[2]} for elementsLocT in elementsLoc[1]}
		elementLabels = {elementsLocT:{elementsLocZ:[] for elementsLocZ in elementsLoc[2]} for elementsLocT in elementsLoc[1]}
	elif colDir == 't':
		elements = {elementsLocR:{elementsLocZ:[] for elementsLocZ in elementsLoc[2]} for elementsLocR in elementsLoc[0]}
		elementLabels = {elementsLocR:{elementsLocZ:[] for elementsLocZ in elementsLoc[2]} for elementsLocR in elementsLoc[0]}
	elif colDir == 'z':
		elements = {elementsLocR:{elementsLocT:[] for elementsLocT in elementsLoc[1]} for elementsLocR in elementsLoc[0]}
		elementLabels = {elementsLocR:{elementsLocT:[] for elementsLocT in elementsLoc[1]} for elementsLocR in elementsLoc[0]}
	return(elements,elementLabels)

def elemNodeLocCol(colDir):
	if colDir == 'r':
		stressesElementNodes = {elementNodesLocR:{elementNodesLocT:{elementNodesLocZ:[] for elementNodesLocZ in elementNodesLoc[3]} for elementNodesLocT in elementNodesLoc[2]} for elementNodesLocR in elementNodesLoc[0]}
	elif colDir == 't':
		stressesElementNodes = {elementNodesLocT:{elementNodesLocR:{elementNodesLocZ:[] for elementNodesLocZ in elementNodesLoc[3]} for elementNodesLocR in elementNodesLoc[1]} for elementNodesLocT in elementNodesLoc[0]}
	elif colDir == 'z':
		stressesElementNodes = {elementNodesLocZ:{elementNodesLocR:{elementNodesLocT:[] for elementNodesLocT in elementNodesLoc[2]} for elementNodesLocR in elementNodesLoc[1]} for elementNodesLocZ in elementNodesLoc[0]}
	return stressesElementNodes

def elemGaussPointsLocCol(colDir):
	if colDir == 'r':
		stressGaussPoints = {elementsLocT:{elementsLocZ:[] for elementsLocZ in elementsLoc[2]} for elementsLocT in elementsLoc[1]}
	elif colDir == 't':
		stressGaussPoints = {elementsLocR:{elementsLocZ:[] for elementsLocZ in elementsLoc[2]} for elementsLocR in elementsLoc[0]}
	elif colDir == 'z':
		stressGaussPoints = {elementsLocR:{elementsLocT:[] for elementsLocT in elementsLoc[1]} for elementsLocR in elementsLoc[0]}
	return stressGaussPoints

def dataGaussPointsLocCol(colDir):
	if colDir == 'r':
		dataOutputGaussPoints = {elementsLocT:{elementsLocZ:[] for elementsLocZ in elementsLoc[2]} for elementsLocT in elementsLoc[1]}
	elif colDir == 't':
		dataOutputGaussPoints = {elementsLocR:{elementsLocZ:[] for elementsLocZ in elementsLoc[2]} for elementsLocR in elementsLoc[0]}
	elif colDir == 'z':
		dataOutputGaussPoints = {elementsLocR:{elementsLocT:[] for elementsLocT in elementsLoc[1]} for elementsLocR in elementsLoc[0]}
	return dataOutputGaussPoints

def stressElementNodesF3NRadialInner(stressesElementNodes,step,elements,conVar):
	try:
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj]][conVar])
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj]][conVar])))
	except:
		raise KeyError
	return(stressesElementNodes)

def stressElementNodesF3NRadialFreeEdgeAxial(stressesElementNodes,step,elements,conVar):
	try:
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]-1):
				getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj+1]][conVar])
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]-1):
				stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii]][elementNodesLoc[3][jj+1]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj+1]][conVar])))
	except:
		try:
			for ii in range(np.array(elementsLoc).shape[1]):
				for jj in range(np.array(elementsLoc).shape[1]-1):
					getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj]][conVar])
			for ii in range(np.array(elementsLoc).shape[1]):
				for jj in range(np.array(elementsLoc).shape[1]-1):
					stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj]][conVar])))
		except:
			raise KeyError
	return(stressesElementNodes)

def stressElementNodesF3NRadialFreeEdgeCirc(stressesElementNodes,step,elements,conVar):
	try:
		for ii in range(np.array(elementsLoc).shape[1]-1):
			for jj in range(np.array(elementsLoc).shape[1]):
				getStressesCoordinatesNode(step,elements[elementsLoc[1][ii+1]][elementsLoc[2][jj]][conVar])
		for ii in range(np.array(elementsLoc).shape[1]-1):
			for jj in range(np.array(elementsLoc).shape[1]):
				stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii+1]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii+1]][elementsLoc[2][jj]][conVar])))
	except:
		try:
			for ii in range(np.array(elementsLoc).shape[1]-1):
				for jj in range(np.array(elementsLoc).shape[1]):
					getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj]][conVar])
			for ii in range(np.array(elementsLoc).shape[1]-1):
				for jj in range(np.array(elementsLoc).shape[1]):
					stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj]][conVar])))
		except:
			raise KeyError
	return(stressesElementNodes)

def stressElementNodesF3NRadialCorner(stressesElementNodes,step,elements,conVar):
	if cornerSubModelQuadrant == 1:
		try:
			for ii in range(np.array(elementsLoc).shape[1]-1):
				for jj in range(np.array(elementsLoc).shape[1]-1):
					stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii+1]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii+1]][elementsLoc[2][jj]][conVar])))
		except:
			raise KeyError
	elif cornerSubModelQuadrant == 2:
		try:
			for ii in range(np.array(elementsLoc).shape[1]-1):
				for jj in range(np.array(elementsLoc).shape[1]-1):
					stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj]][conVar])))
		except:
			raise KeyError
	elif cornerSubModelQuadrant == 3:
		try:
			for ii in range(np.array(elementsLoc).shape[1]-1):
				for jj in range(np.array(elementsLoc).shape[1]-1):
					stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii]][elementNodesLoc[3][jj+1]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj+1]][conVar])))
		except:
			raise KeyError
	else:
		try:
			for ii in range(np.array(elementsLoc).shape[1]-1):
				for jj in range(np.array(elementsLoc).shape[1]-1):
					stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii+1]][elementNodesLoc[3][jj+1]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii+1]][elementsLoc[2][jj+1]][conVar])))
		except:
			raise KeyError
	return(stressesElementNodes)

def stressElementNodesF3N(stressesElementNodes,step,elements,conVar,colDir):
	if colDir == 'r':
		try:
			stressesElementNodes = stressElementNodesF3NRadialInner(stressesElementNodes,step,elements,conVar)
		except:
			try:
				stressesElementNodes = stressElementNodesF3NRadialFreeEdgeAxial(stressesElementNodes,step,elements,conVar)
			except:
				try:
					stressesElementNodes = stressElementNodesF3NRadialFreeEdgeCirc(stressesElementNodes,step,elements,conVar)
				except:
					try:
						stressesElementNodes = stressElementNodesF3NRadialCorner(stressesElementNodes,step,elements,conVar)
					except:
						raise KeyError("Problem in method 'stressElementNodesF3N'")
	elif colDir == 't':
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[1][ii]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[0][ii]][elementsLoc[2][jj]][conVar])))
	elif colDir == 'z':
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[1][ii]][elementNodesLoc[2][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[0][ii]][elementsLoc[1][jj]][conVar])))
	return(stressesElementNodes)

def stressElementNodesL3N(stressesElementNodes,step,elements,conVar,colDir):
	if colDir == 'r':
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				stressesElementNodes[elementNodesLoc[0][1]][elementNodesLoc[2][ii]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][elementsLoc[2][jj]][conVar])))
	elif colDir == 't':
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				stressesElementNodes[elementNodesLoc[0][1]][elementNodesLoc[1][ii]][elementNodesLoc[3][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[0][ii]][elementsLoc[2][jj]][conVar])))
	elif colDir == 'z':
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				stressesElementNodes[elementNodesLoc[0][1]][elementNodesLoc[1][ii]][elementNodesLoc[2][jj]].append((getStressesCoordinatesNode(step,elements[elementsLoc[0][ii]][elementsLoc[1][jj]][conVar])))
	return(stressesElementNodes)

def elemDefLocCenGrav(rEval, tEval, zEval, stressLastFrameElementNodal,compositeSetElements):
	# Berechnung des Schwerpunktes aller Elemente des jeweiligen Elementknotens des Sets und aufgrund der Position, Abspeichern dieser in den jeweiligen Arrays:
	if rEval == 'None':
		colDir = 'r'
	elif tEval == 'None':
		colDir = 't'
	else:
		colDir = 'z'
	elements, elementLabels = elemLocCol(colDir)
	elementsExist = []
	if colDir == 'r':
		for ii in range(len(stressLastFrameElementNodal)):
			# Ueberpruefung, ob die Koordinaten des Elementschwerpunktes bereits berechnet wurden:
			if stressLastFrameElementNodal[ii].elementLabel in elementsExist:
				for jj in range(np.array(elementsLoc).shape[1]):
					for kk in range(np.array(elementsLoc).shape[1]):
						if stressLastFrameElementNodal[ii].elementLabel in elementLabels[elementsLoc[1][jj]][elementsLoc[2][kk]]:
							elements[elementsLoc[1][jj]][elementsLoc[2][kk]].append(stressLastFrameElementNodal[ii])
			else:
				elementsExist.append(stressLastFrameElementNodal[ii].elementLabel)
				# Identifikations des Elementes, dessen Schwerpunkt berechnet werden soll:
				for jj in range(len(compositeSetElements.elements[0])):
					if compositeSetElements.elements[0][jj].label == stressLastFrameElementNodal[ii].elementLabel:
						# Betrachte alle Elementknoten des jeweiligen Elements, bestimme mithilfe der kartesischen Koordinaten ueber eine Mittelung die Radien/Winkel der Gauss-Punkte
						# und speichere Elemente, deren Positionen der Gauss-Punkte gewuenschte Kriterien erfuellen, zur Weiterverarbeitung ab:
						nodeLabelConnect = np.array(list(compositeSetElements.elements[0][jj].connectivity))
						nodeLabelConnectInst = compositeSubModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].instance
						xCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[0] for kk in nodeLabelConnect])
						yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
						zCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[2] for kk in nodeLabelConnect])
						
						tCoord = ct.cart2pol_theta(xCoord,yCoord)
						sumNodesPosThetaMean,sumNodesPosZMean  = 180*np.mean(tCoord)/np.pi, np.mean(zCoord)
						
						# Abspeichern der Element-Labels der Elemente, welche die gewuenschten Kriterien erfuellen:
						if sumNodesPosThetaMean > tEval:
							if sumNodesPosZMean > zEval:
								elements['t > t_eval']['z > z_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['t > t_eval']['z > z_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
							else:
								elements['t > t_eval']['z < z_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['t > t_eval']['z < z_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
						else:
							if sumNodesPosZMean > zEval:
								elements['t < t_eval']['z > z_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['t < t_eval']['z > z_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
							else:
								elements['t < t_eval']['z < z_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['t < t_eval']['z < z_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
	elif colDir == 't':
		for ii in range(len(stressLastFrameElementNodal)):
			# Ueberpruefung, ob die Koordinaten des Elementschwerpunktes bereits berechnet wurden:
			if stressLastFrameElementNodal[ii].elementLabel in elementsExist:
				for jj in range(np.array(elementsLoc).shape[1]):
					for kk in range(np.array(elementsLoc).shape[1]):
						if stressLastFrameElementNodal[ii].elementLabel in elementLabels[elementsLoc[0][jj]][elementsLoc[2][kk]]:
							elements[elementsLoc[0][jj]][elementsLoc[2][kk]].append(stressLastFrameElementNodal[ii])
			else:
				elementsExist.append(stressLastFrameElementNodal[ii].elementLabel)
				# Identifikations des Elementes, dessen Schwerpunkt berechnet werden soll:
				for jj in range(len(compositeSetElements.elements[0])):
					if compositeSetElements.elements[0][jj].label == stressLastFrameElementNodal[ii].elementLabel:
						# Betrachte alle Elementknoten des jeweiligen Elements, bestimme mithilfe der kartesischen Koordinaten ueber eine Mittelung die Radien/Winkel der Gauss-Punkte
						# und speichere Elemente, deren Positionen der Gauss-Punkte gewuenschte Kriterien erfuellen, zur Weiterverarbeitung ab:
						nodeLabelConnect = np.array(list(compositeSetElements.elements[0][jj].connectivity))
						nodeLabelConnectInst = compositeSubModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].instance
						xCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[0] for kk in nodeLabelConnect])
						yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
						zCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[2] for kk in nodeLabelConnect])
						
						rCoord = ct.cart2pol_radius(xCoord,yCoord)
						sumNodesPosRMean,sumNodesPosZMean  = np.mean(rCoord), np.mean(zCoord)
						
						# Abspeichern der Element-Labels der Elemente, welche die gewuenschten Kriterien erfuellen:
						if sumNodesPosRMean > rEval:
							if sumNodesPosZMean > zEval:
								elements['r > r_eval']['z > z_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['r > r_eval']['z > z_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
							else:
								elements['r > r_eval']['z < z_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['r > r_eval']['z < z_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
						else:
							if sumNodesPosZMean > zEval:
								elements['r < r_eval']['z > z_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['r < r_eval']['z > z_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
							else:
								elements['r < r_eval']['z < z_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['r < r_eval']['z < z_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
	elif colDir == 'z':
		for ii in range(len(stressLastFrameElementNodal)):
			# Ueberpruefung, ob die Koordinaten des Elementschwerpunktes bereits berechnet wurden:
			if stressLastFrameElementNodal[ii].elementLabel in elementsExist:
				for jj in range(np.array(elementsLoc).shape[1]):
					for kk in range(np.array(elementsLoc).shape[1]):
						if stressLastFrameElementNodal[ii].elementLabel in elementLabels[elementsLoc[0][jj]][elementsLoc[1][kk]]:
							elements[elementsLoc[0][jj]][elementsLoc[1][kk]].append(stressLastFrameElementNodal[ii])
			else:
				elementsExist.append(stressLastFrameElementNodal[ii].elementLabel)
				# Identifikations des Elementes, dessen Schwerpunkt berechnet werden soll:
				for jj in range(len(compositeSetElements.elements[0])):
					if compositeSetElements.elements[0][jj].label == stressLastFrameElementNodal[ii].elementLabel:
						# Betrachte alle Elementknoten des jeweiligen Elements, bestimme mithilfe der kartesischen Koordinaten ueber eine Mittelung die Radien/Winkel der Gauss-Punkte
						# und speichere Elemente, deren Positionen der Gauss-Punkte gewuenschte Kriterien erfuellen, zur Weiterverarbeitung ab:
						nodeLabelConnect = np.array(list(compositeSetElements.elements[0][jj].connectivity))
						nodeLabelConnectInst = compositeSubModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].instance
						xCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[0] for kk in nodeLabelConnect])
						yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
						
						rCoord, tCoord = ct.cart2pol_radius(xCoord,yCoord), ct.cart2pol_theta(xCoord,yCoord)
						sumNodesPosRMean,sumNodesPosPhiMean  = np.mean(rCoord), 180*np.mean(tCoord)/np.pi
						
						# Abspeichern der Element-Labels der Elemente, welche die gewuenschten Kriterien erfuellen:
						if sumNodesPosPhiMean > tEval:
							if sumNodesPosRMean > rEval:
								elements['r > r_eval']['t > t_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['r > r_eval']['t > t_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
							else:
								elements['r < r_eval']['t > t_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['r < r_eval']['t > t_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
						else:
							if sumNodesPosRMean > rEval:
								elements['r > r_eval']['t < t_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['r > r_eval']['t < t_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
							else:
								elements['r < r_eval']['t < t_eval'].append(stressLastFrameElementNodal[ii])
								elementLabels['r < r_eval']['t < t_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
	return(elements, elementLabels)

def elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,boolInterface,colDir):
	# Bestimmung des Elementtyps:
	if colDir == 'r':
		try:
			elementTyp = elements[elementsLoc[1][0]][elementsLoc[2][0]][0].baseElementType
		except:
			try:
				elementTyp = elements[elementsLoc[1][0]][elementsLoc[2][1]][0].baseElementType
			except:
				try:
					elementTyp = elements[elementsLoc[1][1]][elementsLoc[2][0]][0].baseElementType
				except:
					elementTyp = elements[elementsLoc[1][1]][elementsLoc[2][1]][0].baseElementType
	elif colDir == 't':
		elementTyp = elements[elementsLoc[0][0]][elementsLoc[2][0]][0].baseElementType
	elif colDir == 'z':
		elementTyp = elements[elementsLoc[0][0]][elementsLoc[1][0]][0].baseElementType
	stressesElementNodes = elemNodeLocCol(colDir)
	laufVarElementNodesFieldOutput = 0
	for ii in range(len(sortedElementNodes)):
		if ii == 0:
				stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
				laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput + 1
		elif ii == len(sortedElementNodes)-1:
				stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
		else:
			if elementTyp == 'C3D8' or elementTyp == 'C3D8R':
				stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
				laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+2
			elif elementTyp == 'C3D20' or elementTyp == 'C3D20R':
				if ii % 2 == 0:
					stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
					laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+2
				else:
					stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
					laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+1
	if boolInterface:
		laufVarElementNodesFieldOutput = 0
		for ii in range(len(sortedElementNodes)):
			if ii == len(sortedElementNodes)-2:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
						laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+1
			elif ii == len(sortedElementNodes)-1:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
			else:
				if elementTyp == 'C3D8' or elementTyp == 'C3D8R':
					if ii == 0:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
						laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput + 2
					elif ii == len(sortedElementNodes)-1:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
					else:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
						laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+2
				elif elementTyp == 'C3D20' or elementTyp == 'C3D20R':
					if ii == 0:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
						laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput + 1
					elif ii == 1:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
						laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput + 2
					elif ii % 2 == 0:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
						laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+1
					else:
						stressesElementNodes = stressElementNodesL3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
						laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+2
	return stressesElementNodes

def elemLocInterface(stressLastFrameGaussPoint,elementLabels,colDir):
	# Bestimmung der Position der Elemente und abspeichern relevater Daten:
	stressGaussPoints = elemGaussPointsLocCol(colDir)
	dataOutputGaussPoints = dataGaussPointsLocCol(colDir)
	if colDir == 'r':
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				for mm in range(len(elementLabels[elementsLoc[1][ii]][elementsLoc[2][jj]])):
					stressGaussPointsTemp = elemGaussPointsLocCol(colDir)
					for nn in range(len(stressLastFrameGaussPoint)):
						if elementLabels[elementsLoc[1][ii]][elementsLoc[2][jj]][mm] == stressLastFrameGaussPoint[nn].elementLabel:
							stressGaussPointsTemp[elementsLoc[1][ii]][elementsLoc[2][jj]].append((stressLastFrameGaussPoint[nn], stressLastFrameGaussPoint[nn].elementLabel, stressLastFrameGaussPoint[nn].integrationPoint))
					stressGaussPointsTemp[elementsLoc[1][ii]][elementsLoc[2][jj]] = sorted(stressGaussPointsTemp[elementsLoc[1][ii]][elementsLoc[2][jj]],key=getIntegrationPoint)
					for nn in range(len(stressGaussPointsTemp[elementsLoc[1][ii]][elementsLoc[2][jj]])):
						stressGaussPoints[elementsLoc[1][ii]][elementsLoc[2][jj]].append(stressGaussPointsTemp[elementsLoc[1][ii]][elementsLoc[2][jj]][nn][0])
		# Abspeichern der Daten, die letztendlich ausgegeben werden sollen:
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				for kk in range(len(stressGaussPoints[elementsLoc[1][ii]][elementsLoc[2][jj]])):
					dataOutputGaussPoints[elementsLoc[1][ii]][elementsLoc[2][jj]].append((getStressesIntegrationPoint(stressGaussPoints[elementsLoc[1][ii]][elementsLoc[2][jj]][kk])))
	elif colDir == 't':
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				for mm in range(len(elementLabels[elementsLoc[0][ii]][elementsLoc[2][jj]])):
					stressGaussPointsTemp = elemGaussPointsLocCol(colDir)
					for nn in range(len(stressLastFrameGaussPoint)):
						if elementLabels[elementsLoc[0][ii]][elementsLoc[2][jj]][mm] == stressLastFrameGaussPoint[nn].elementLabel:
							stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[2][jj]].append((stressLastFrameGaussPoint[nn], stressLastFrameGaussPoint[nn].elementLabel, stressLastFrameGaussPoint[nn].integrationPoint))
					stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[2][jj]] = sorted(stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[2][jj]],key=getIntegrationPoint)
					for nn in range(len(stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[2][jj]])):
						stressGaussPoints[elementsLoc[0][ii]][elementsLoc[2][jj]].append(stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[2][jj]][nn][0])
		# Abspeichern der Daten, die letztendlich ausgegeben werden sollen:
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				for kk in range(len(stressGaussPoints[elementsLoc[0][ii]][elementsLoc[2][jj]])):
					dataOutputGaussPoints[elementsLoc[0][ii]][elementsLoc[2][jj]].append((getStressesIntegrationPoint(stressGaussPoints[elementsLoc[0][ii]][elementsLoc[2][jj]][kk])))
	elif colDir == 'z':
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				for mm in range(len(elementLabels[elementsLoc[0][ii]][elementsLoc[1][jj]])):
					stressGaussPointsTemp = elemGaussPointsLocCol(colDir)
					for nn in range(len(stressLastFrameGaussPoint)):
						if elementLabels[elementsLoc[0][ii]][elementsLoc[1][jj]][mm] == stressLastFrameGaussPoint[nn].elementLabel:
							stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[1][jj]].append((stressLastFrameGaussPoint[nn], stressLastFrameGaussPoint[nn].elementLabel, stressLastFrameGaussPoint[nn].integrationPoint))
					stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[1][jj]] = sorted(stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[1][jj]],key=getIntegrationPoint)
					for nn in range(len(stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[1][jj]])):
						stressGaussPoints[elementsLoc[0][ii]][elementsLoc[1][jj]].append(stressGaussPointsTemp[elementsLoc[0][ii]][elementsLoc[1][jj]][nn][0])
		# Abspeichern der Daten, die letztendlich ausgegeben werden sollen:
		for ii in range(np.array(elementsLoc).shape[1]):
			for jj in range(np.array(elementsLoc).shape[1]):
				for kk in range(len(stressGaussPoints[elementsLoc[0][ii]][elementsLoc[1][jj]])):
					dataOutputGaussPoints[elementsLoc[0][ii]][elementsLoc[1][jj]].append((getStressesIntegrationPoint(stressGaussPoints[elementsLoc[0][ii]][elementsLoc[1][jj]][kk])))
	return dataOutputGaussPoints

def dtypeOutputElementNodeFnc():
	return [('sigrr',float),('sigtt',float),('sigzz',float),('tautz',float),('taurz',float),('taurt',float),('xCoords',float),('yCoords',float),('zCoords',float),('r',float),('theta',float), ('labelNode',int),('labelElement',int)]

def dtypeOutputGaussPointFnc():
	return [('labelElement',int),('integrationPoint',int),('sigrr',float),('sigtt',float),('sigzz',float),('tautz',float),('taurz',float),('taurt',float)]

'''
Definierte Knoten-/Element-Sets werden mithilfe folgender Funktionen im PApply loadostprocessing verarbeitet, d.h. gewuenschte Daten, wie z. B. Spannungen o. Verschiebungen,  werden in einer *.txt-File abgespeichert.
Hinsichtlich der Spannungen werden diese sowohl an Elementknoten(durch Extrapolation) als auch in den Integrations-Punkten ausgegeben:
'''
def stressesExtrapolatedAndIntegrationPointInterfaceAxial(set, step, iInterfaceEval, rEval,  tEval):
	colDir = 'z'
	compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(),step,colDir)
	
	# sortedElementNodeLabels = np.array([ii[0] for ii in sortedElementNodes])
	stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
	
	elements, elementLabels = elemDefLocCenGrav(rEval, tEval, 'None', stressLastFrameElementNodal,compositeSetElements)
	
	# Axiale Sortierung der Elemente:
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			elements[elementsLoc[0][ii]][elementsLoc[1][jj]] = sortElements(sortedElementNodes,elements[elementsLoc[0][ii]][elementsLoc[1][jj]])
	
	stressesElementNodes = elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,True,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	# 1: GreaterAngleEval
	# 2: SmallerAngleEval
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:	
	dataOutputElemNodes = {elementNodesLocNo:{elementNodesLocR:{elementNodesLocT:[] for elementNodesLocT in elementNodesLoc[2]} for elementNodesLocR in elementNodesLoc[1]} for elementNodesLocNo in elementNodesLoc[0]}
	
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]][elementNodesLoc[2][0]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_TS_1'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]][elementNodesLoc[2][1]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_TG_1'
	
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]][elementNodesLoc[2][0]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_TS_1'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]][elementNodesLoc[2][1]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_TG_1'
	
	dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][0]][elementNodesLoc[2][0]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_TS_2'
	dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][0]][elementNodesLoc[2][1]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_TG_2'
	
	dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][1]][elementNodesLoc[2][0]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_TS_2'
	dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][1]][elementNodesLoc[2][1]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_TG_2'
	dtypeOutputElementNode = dtypeOutputElementNodeFnc()
	for ii in range(np.array(elementNodesLoc).shape[1]):
		for jj in range(np.array(elementNodesLoc).shape[1]):
			for kk in range(np.array(elementNodesLoc).shape[1]):
				np.savetxt(dataOutputElemNodes[elementNodesLoc[0][ii]][elementNodesLoc[1][jj]][elementNodesLoc[2][kk]]+'.txt', np.sort(np.array(stressesElementNodes[elementNodesLoc[0][ii]][elementNodesLoc[1][jj]][elementNodesLoc[2][kk]], dtype=dtypeOutputElementNode),order='zCoords'))
	
	#---------------------------------------------------------------
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:
	dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint, elementLabels,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
	dataOutputGaussPointsTxt[elementsLoc[0][1]][elementsLoc[1][1]]= set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_TG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][0]][elementsLoc[1][1]]= set + '_Layer_' + str(plyAngle[iInterfaceEval-1]) + '_TG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][1]][elementsLoc[1][0]]= set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_TS_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][0]][elementsLoc[1][0]]= set + '_Layer_' + str(plyAngle[iInterfaceEval-1]) + '_TS_IntPo'
	dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			np.savetxt(dataOutputGaussPointsTxt[elementsLoc[0][ii]][elementsLoc[1][jj]]+'.txt', np.array(dataOutputGaussPoints[elementsLoc[0][ii]][elementsLoc[1][jj]], dtype=dtypeOutputGaussPoint))

def stressesExtrapolatedAndIntegrationPointInterfaceCirc(set, step, iInterfaceEval, rEval,  zEval):
	colDir = 't'
	compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(),step,colDir)
	
	# sortedElementNodeLabels = np.array([ii[0] for ii in sortedElementNodes])
	stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
	
	elements, elementLabels = elemDefLocCenGrav(rEval, 'None', zEval, stressLastFrameElementNodal,compositeSetElements)
	
	# Sortierung der Elemente Umfangsrichtung:
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			elements[elementsLoc[0][ii]][elementsLoc[2][jj]] = sortElements(sortedElementNodes,elements[elementsLoc[0][ii]][elementsLoc[2][jj]])
	
	stressesElementNodes = elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,True,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	# 1: GreaterAngleEval
	# 2: SmallerAngleEval
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:	
	dataOutputElemNodes = {elementNodesLocNo:{elementNodesLocR:{elementNodesLocZ:[] for elementNodesLocZ in elementNodesLoc[3]} for elementNodesLocR in elementNodesLoc[1]} for elementNodesLocNo in elementNodesLoc[0]}
	
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]][elementNodesLoc[3][0]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_ZS_1'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]][elementNodesLoc[3][1]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_ZG_1'
	
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]][elementNodesLoc[3][0]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_ZS_1'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]][elementNodesLoc[3][1]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_ZG_1'
	
	dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][0]][elementNodesLoc[3][0]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_ZS_2'
	dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][0]][elementNodesLoc[3][1]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_ZG_2'
	
	dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][1]][elementNodesLoc[3][0]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_ZS_2'
	dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][1]][elementNodesLoc[3][1]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_ZG_2'
	dtypeOutputElementNode = dtypeOutputElementNodeFnc()
	for ii in range(np.array(elementNodesLoc).shape[1]):
		for jj in range(np.array(elementNodesLoc).shape[1]):
			for kk in range(np.array(elementNodesLoc).shape[1]):
				np.savetxt(dataOutputElemNodes[elementNodesLoc[0][ii]][elementNodesLoc[1][jj]][elementNodesLoc[3][kk]]+'.txt', np.sort(np.array(stressesElementNodes[elementNodesLoc[0][ii]][elementNodesLoc[1][jj]][elementNodesLoc[3][kk]], dtype=dtypeOutputElementNode),order='theta'))
	
	#---------------------------------------------------------------
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:
	dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint, elementLabels,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
	dataOutputGaussPointsTxt[elementsLoc[0][1]][elementsLoc[2][1]]= set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_ZG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][0]][elementsLoc[2][1]]= set + '_Layer_' + str(plyAngle[iInterfaceEval-1]) + '_ZG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][1]][elementsLoc[2][0]]= set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_ZS_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][0]][elementsLoc[2][0]]= set + '_Layer_' + str(plyAngle[iInterfaceEval-1]) + '_ZS_IntPo'
	dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			np.savetxt(dataOutputGaussPointsTxt[elementsLoc[0][ii]][elementsLoc[2][jj]]+'.txt', np.array(dataOutputGaussPoints[elementsLoc[0][ii]][elementsLoc[2][jj]], dtype=dtypeOutputGaussPoint))

def stressesExtrapolatedAndIntegrationPointLayerAxial(set, step, rEval, tEval):
	colDir = 'z'
	compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(),step,colDir)
	
	# sortedElementNodeLabels = np.array([ii[0] for ii in sortedElementNodes])
	stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
	
	elements, elementLabels = elemDefLocCenGrav(rEval, tEval, 'None', stressLastFrameElementNodal,compositeSetElements)
	
	# Axiale Sortierung der Elemente:
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			elements[elementsLoc[0][ii]][elementsLoc[1][jj]] = sortElements(sortedElementNodes,elements[elementsLoc[0][ii]][elementsLoc[1][jj]])
	
	stressesElementNodes = elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,False,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	# 1: GreaterAngleEval
	# 2: SmallerAngleEval
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:	
	dataOutputElemNodes = {elementNodesLoc[0][0]:{elementNodesLocR:{elementNodesLocT:[] for elementNodesLocT in elementNodesLoc[2]} for elementNodesLocR in elementNodesLoc[1]}}
	
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]][elementNodesLoc[2][1]] = set + '_RG_TG'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]][elementNodesLoc[2][1]] = set + '_RS_TG'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]][elementNodesLoc[2][0]] = set + '_RG_TS'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]][elementNodesLoc[2][0]] = set + '_RS_TS'
	
	dtypeOutputElementNode = dtypeOutputElementNodeFnc()
	for jj in range(np.array(elementNodesLoc).shape[1]):
		for kk in range(np.array(elementNodesLoc).shape[1]):
			np.savetxt(dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][jj]][elementNodesLoc[2][kk]]+'.txt', np.sort(np.array(stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[1][jj]][elementNodesLoc[2][kk]], dtype=dtypeOutputElementNode),order='zCoords'))
	
	
	#---------------------------------------------------------------
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:
	dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint, elementLabels,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
	dataOutputGaussPointsTxt[elementsLoc[0][1]][elementsLoc[1][1]]= set + '_RG_TG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][0]][elementsLoc[1][1]]= set + '_RS_TG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][1]][elementsLoc[1][0]]= set + '_RG_TS_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][0]][elementsLoc[1][0]]= set + '_RS_TS_IntPo'
	dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			np.savetxt(dataOutputGaussPointsTxt[elementsLoc[0][ii]][elementsLoc[1][jj]]+'.txt', np.array(dataOutputGaussPoints[elementsLoc[0][ii]][elementsLoc[1][jj]], dtype=dtypeOutputGaussPoint))

def stressesExtrapolatedAndIntegrationPointLayerCirc(set, step, rEval, zEval):
	colDir = 't'
	compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(),step,colDir)
	
	# sortedElementNodeLabels = np.array([ii[0] for ii in sortedElementNodes])
	stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
	
	elements, elementLabels = elemDefLocCenGrav(rEval, 'None', zEval, stressLastFrameElementNodal,compositeSetElements)
	
	# Axiale Sortierung der Elemente:
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			elements[elementsLoc[0][ii]][elementsLoc[1][jj]] = sortElements(sortedElementNodes,elements[elementsLoc[0][ii]][elementsLoc[2][jj]])
	
	stressesElementNodes = elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,False,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	# 1: GreaterAngleEval
	# 2: SmallerAngleEval
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:	
	dataOutputElemNodes = {elementNodesLoc[0][0]:{elementNodesLocR:{elementNodesLocZ:[] for elementNodesLocZ in elementNodesLoc[3]} for elementNodesLocR in elementNodesLoc[1]}}
	
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]][elementNodesLoc[3][1]] = set + '_RG_ZG'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]][elementNodesLoc[3][1]] = set + '_RS_ZG'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]][elementNodesLoc[3][0]] = set + '_RG_ZS'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]][elementNodesLoc[3][0]] = set + '_RS_ZS'
	
	dtypeOutputElementNode = dtypeOutputElementNodeFnc()
	for jj in range(np.array(elementNodesLoc).shape[1]):
		for kk in range(np.array(elementNodesLoc).shape[1]):
			np.savetxt(dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][jj]][elementNodesLoc[3][kk]]+'.txt', np.sort(np.array(stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[1][jj]][elementNodesLoc[3][kk]], dtype=dtypeOutputElementNode),order='theta'))
	
	
	#---------------------------------------------------------------
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:
	dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint, elementLabels,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
	dataOutputGaussPointsTxt[elementsLoc[0][1]][elementsLoc[2][1]]= set + '_RG_ZG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][0]][elementsLoc[2][1]]= set + '_RS_ZG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][1]][elementsLoc[2][0]]= set + '_RG_ZS_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[0][0]][elementsLoc[2][0]]= set + '_RS_ZS_IntPo'
	dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			np.savetxt(dataOutputGaussPointsTxt[elementsLoc[0][ii]][elementsLoc[2][jj]]+'.txt', np.array(dataOutputGaussPoints[elementsLoc[0][ii]][elementsLoc[2][jj]], dtype=dtypeOutputGaussPoint))

#---------------------------------------------------------------
'''
Definierte Knoten-/Element-Sets werden mithilfe folgender Funktionen im Postprocessing verarbeitet, d.h. gewuenschte Daten, wie z. B. Spannungen o. Verschiebungen,  werden in einer *.txt-File abgespeichert.
Hinsichtlich der Spannungen werden diese sowohl an Elementknoten(durch Extrapolation) als auch in den Integrations-Punkten ausgegeben:
'''
def stressesExtrapolatedAndIntegrationPointRadial(set, step, tEval, zEval):
	colDir = 'r'
	compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(),step,colDir)
	
	stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
	
	elements, elementLabels = elemDefLocCenGrav('None', tEval, zEval, stressLastFrameElementNodal,compositeSetElements)
	
	# Radiale Sortierung der Elemente:
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			elements[elementsLoc[1][ii]][elementsLoc[2][jj]] = sortElements(sortedElementNodes,elements[elementsLoc[1][ii]][elementsLoc[2][jj]])
	
	stressesElementNodes = elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,False,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	# 1: GreaterAngleEval
	# 2: SmallerAngleEval
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:	
	dataOutputElemNodes = {elementNodesLoc[0][0]:{elementNodesLocT:{elementNodesLocZ:[] for elementNodesLocZ in elementNodesLoc[3]} for elementNodesLocT in elementNodesLoc[2]}}
	
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][1]][elementNodesLoc[3][1]] = set + '_TG_ZG'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][0]][elementNodesLoc[3][1]] = set + '_TS_ZG'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][1]][elementNodesLoc[3][0]] = set + '_TG_ZS'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][0]][elementNodesLoc[3][0]] = set + '_TS_ZS'
	
	dtypeOutputElementNode = dtypeOutputElementNodeFnc()
	for jj in range(np.array(elementNodesLoc).shape[1]):
		for kk in range(np.array(elementNodesLoc).shape[1]):
			if bool(stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]][elementNodesLoc[3][kk]]):
				np.savetxt(dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]][elementNodesLoc[3][kk]]+'.txt', np.sort(np.array(stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]][elementNodesLoc[3][kk]], dtype=dtypeOutputElementNode),order='r'))
		
	#---------------------------------------------------------------
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:
	dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint,elementLabels,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
	dataOutputGaussPointsTxt[elementsLoc[1][1]][elementsLoc[2][1]]= set + '_TG_ZG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[1][0]][elementsLoc[2][1]]= set + '_TS_ZG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[1][1]][elementsLoc[2][0]]= set + '_TG_ZS_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[1][0]][elementsLoc[2][0]]= set + '_TS_ZS_IntPo'
	dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
	for ii in range(np.array(elementsLoc).shape[1]):
		for jj in range(np.array(elementsLoc).shape[1]):
			if bool(dataOutputGaussPoints[elementsLoc[1][ii]][elementsLoc[2][jj]]):
				np.savetxt(dataOutputGaussPointsTxt[elementsLoc[1][ii]][elementsLoc[2][jj]]+'.txt', np.array(dataOutputGaussPoints[elementsLoc[1][ii]][elementsLoc[2][jj]], dtype=dtypeOutputGaussPoint))

#---------------------------------------------------------------
'''
Definierte Knoten-/Element-Sets werden mithilfe folgender Funktionen im Postprocessing verarbeitet, d.h. gewuenschte Daten, wie z. B. Spannungen o. Verschiebungen,  werden in einer *.txt-File abgespeichert.
Hinsichtlich der Spannungen werden diese sowohl an Elementknoten(durch Extrapolation) als auch in den Integrations-Punkten ausgegeben:
'''
def stressesExtrapolatedAndIntegrationPointSurface(set, step):
	colDir = 'tz'
	compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(),step,colDir)
	
	stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
	
	stressesElementNodes,stressesElementNodesTemp = [],[]
	for ii in range(len(stressLastFrameElementNodal)):
		if not bool(stressesElementNodesTemp):
			stressesElementNodesTemp.append(list(getStressesCoordinatesNode(step,stressLastFrameElementNodal[ii])))
		elif ii == len(stressLastFrameElementNodal)-1:
			stressesElementNodesTemp.append(list(getStressesCoordinatesNode(step,stressLastFrameElementNodal[ii])))
			stressesElementNodesTemp = np.array(stressesElementNodesTemp)
			stressesElementNodesMeanStresses = [stressesElementNodesTemp[:,jj].mean() for jj in range(6)]
			for jj in range(6,13):
				if jj not in (11,12):
					stressesElementNodesMeanStresses.append(stressesElementNodesTemp[-1,jj])
				else:
					stressesElementNodesMeanStresses.append(int(stressesElementNodesTemp[-1,jj]))
			stressesElementNodes.append(stressesElementNodesMeanStresses)
		elif getStressesCoordinatesNode(step,stressLastFrameElementNodal[ii])[-2] == int(stressesElementNodesTemp[-1][-2]):
			stressesElementNodesTemp.append(list(getStressesCoordinatesNode(step,stressLastFrameElementNodal[ii])))
		else:
			stressesElementNodesTemp = np.array(stressesElementNodesTemp)
			stressesElementNodesMeanStresses = [stressesElementNodesTemp[:,jj].mean() for jj in range(6)]
			for jj in range(6,13):
				if jj not in (11,12):
					stressesElementNodesMeanStresses.append(stressesElementNodesTemp[-1,jj])
				else:
					stressesElementNodesMeanStresses.append(int(stressesElementNodesTemp[-1,jj]))
			stressesElementNodes.append(stressesElementNodesMeanStresses)
			stressesElementNodesTemp = [list(getStressesCoordinatesNode(step,stressLastFrameElementNodal[ii]))]
	
	outputStressesElementNode = [tuple(stressesElementNode) for stressesElementNode in stressesElementNodes]
	dataOutputElemNodes = set
	dtypeOutputElementNode = dtypeOutputElementNodeFnc()
	np.savetxt(dataOutputElemNodes+'.txt', np.array(outputStressesElementNode, dtype=dtypeOutputElementNode))
	

#---------------------------------------------------------------
# Names des Modells:
modelName = '0_90_S_R4_L4h_M1_AO90_AE45_CFK'

# Modellparameter:
# Schichtwinkel der einzelnen physikalischen Schichten:
plyAngle = [0,90]

# Definition des Pfades zur Sicherung der gesamten FE-Analyse: 
analysis_Path = 'D:\\psingh\\MT\\Fe'

# Bestimmung des Zeitpunktes der FE-Analyse:
analysis_currentDateTime = dt.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

# Erstellen eines neuen Ordners in dem angebeben Pfad und Definition dessen als Arbeitsverzeichnis:
analysis_newPath = analysis_Path + '\\' + modelName + '_' + analysis_currentDateTime
os.makedirs(analysis_newPath)
os.chdir(analysis_newPath)

# Anzahl physikalischer Schichten:
N = len(plyAngle)

# Auswertungsinterfaces (Interfaces zweier aufeinanderfolgender Schichten mit unterschiedlichen Faserorientierungswinkeln):
iInterfaceEval = list(range(1,N))

compositeMaterialName = 'cfk' # cfk, gfk, cfkDuro

if compositeMaterialName == 'cfk':
	E1, E2, E3 = 7460.0, 118148.0, 7460.0
	Nu12, Nu13, Nu23 = 0.021, 0.37, 0.34
	G12, G13, G23 = 4800.0, 2701.0, 4800.0
	alpha11,alpha22,alpha33=2.6e-5,-1.0e-6,2.6e-5
	dL = 0.125
elif compositeMaterialName == 'gfk':
	E1, E2, E3 = 9552.6, 39296.0, 9552.6
	Nu21, Nu13, Nu23 = 0.29, 0.38, 0.29
	Nu12=E1/E2*Nu21
	G12, G13, G23 = 3080.5, 3449.0, 3080.5
	alpha11,alpha22,alpha33=2.6e-5,8.6e-6,2.6e-5
	dL = 0.190
elif compositeMaterialName == 'cfkDuro':
	E1, E2, E3 = 6895.0, 172375.0, 6895.0
	Nu21, Nu13, Nu23 = 0.25, 0.25, 0.25
	Nu12=E1/E2*Nu21
	G12,G13,G23=3448.0,1379.0,3448.0
	alpha11,alpha22,alpha33=2.6e-5,-1.0e-6,2.6e-5
	dL = 0.125
else:
	pass

# R/h:
rt = 4.0

# L/h:
lt = 4.0

# Laminatdicke:
h = N*dL

# Innenradius:
r0 = rt*h-h/2

# Laminattiefe:
length = lt*h

# Schichtdicke - Aluminium 2024-T3:
dL2024 = 0.3

# Aluminium 2024-T3:
boolHybrid = False
metalMaterialName='Aluminum 2024-T3'
E2024 = 71700
Nu2024 = 0.33
alpha2024 = 23.1e-6
beta2024 = 0.0

# Betrag des angreifenden Momentenflusses:
bendingMoment = 0.0
if bendingMoment != 0.0:
	boolBendingMoment = True
else:
	boolBendingMoment = False

# Betrag der angreifenden Radiallast:
radialForce = 0.0
if radialForce != 0.0:
	boolRadialForce = True
else:
	boolRadialForce = False

# Betrag des angreifenden Normalkraftflusses:
circForce = 0.0
if circForce != 0.0:
	boolCircForce = True
else:
	boolCircForce = False

# Betrag des angreifenden Innen-/Aussendruckes:
OuterPressure, InnerPressure = 0.0, 0.0
if OuterPressure != 0.0 or InnerPressure != 0.0:
	boolPressure = True
else:
	boolPressure = False

# Betrag des angreifenden zylindrischen Innen-/Aussendruckes:
cylindricalOuterPressure, cylindricalInnerPressure = 0.0, 0.0
if cylindricalOuterPressure != 0.0 or cylindricalInnerPressure != 0.0:
	boolCylindricalPressure = True
else:
	boolCylindricalPressure = False

# Temperaturdifferenz:
tempDif = 20.0
if tempDif != 0.0:
	boolTempDif = True
else:
	boolTempDif = False

# Feuchtekonzentration:
moistDif = 0.0
if moistDif != 0.0:
	boolMoistDif = True
else:
	boolMoistDif = False

if (any([boolBendingMoment, boolRadialForce, boolCircForce]) and any([boolPressure, boolCylindricalPressure])) or (any([boolBendingMoment, boolRadialForce, boolCircForce]) and any([boolTempDif, boolMoistDif])) or (any([boolPressure, boolCylindricalPressure]) and any([boolTempDif, boolMoistDif])):
	raise ValueError('Please check the modelled structural situation!')
else:
	pass

if any([boolTempDif, boolMoistDif]):
	boolAxialSubModel = False
	boolCircSubModel = not boolAxialSubModel
	boolCornerSubModel, cornerSubModelQuadrant = False, 1
	boolSurfaceSubModel, surfaceSubModelInterface = True, 1
	# Startwinkel des Modells, Winkelbereich des Sub-Modells, Oeffnungswinkel, Auswertungswinkel 
	angleStart,angleSubModel,angleOpening = 0.0, np.pi/8,np.pi/3
	# Auswertungswinkel:
	angleEval = 1*angleOpening/4
	# Breitenpartitionierungen:
	axialPartitionPercent = [0.25,0.5,0.75]
	# Submodel-Winkel:
	if (angleStart+angleEval-angleSubModel/2  < angleStart) or (angleStart+angleEval+angleSubModel/2  > angleStart+angleOpening):
		raise ValueError('Check the circumferential partition!')
	else:
		angleSubModelStart = angleStart
		angleSubModelEnd = angleStart+angleOpening
else:
	boolAxialSubModel = True
	boolCircSubModel, boolCornerSubModel, boolSurfaceSubModel  = False, False, False
	# Startwinkel des Modells, Winkelbereich des Sub-Modells, Oeffnungswinkel, Auswertungswinkel 
	angleStart,angleSubModel,angleOpening = 0.0, np.pi/8,np.pi/3
	# Auswertungswinkel:
	angleEval = 2*angleOpening/4
	# Submodel-Winkel:
	angleSubModelStart = angleStart+angleEval-angleSubModel/2
	angleSubModelEnd = angleStart+angleEval+angleSubModel/2
	# Breitenpartitionierungen:
	axialPartitionPercent = [0.1,0.5,0.9]

axialEval = axialPartitionPercent[int(len(axialPartitionPercent)/2)]*length

# Start der FE-(Konvergenz-)Analyse:
#---------------------------------------------------------------

# Definiere eine neue Model-Database:
Mdb()
session.viewports['Viewport: 1'].setValues(displayedObject=None)

# Definiere das Abaqus-Model:
mdb.models.changeKey(fromName='Model-1', toName=modelName)
compositeModel = mdb.models[modelName]

def radialGeometryParameters(N,rk,iInterfaceEval):
	# Radien der Mittelflaechen der Laminat-Einzelschichten:
	rm = [(rk[ii]+rk[ii+1])/2 for ii in range(N)]
	# Radien aller Interfaces:
	rInterfaceEval =[rk[iInterEval] for iInterEval in iInterfaceEval]
	# Radien aller Partitionierungsinterfaces:
	rInterfaceAll =[]
	for ii in range(N+1):
		if ii < N:
			rInterfaceAll.append(rk[ii])
			rInterfaceAll.append(rm[ii])
		else:
			rInterfaceAll.append(rk[ii])
	# Mittlere Radien aller partitionierten physikalischen Schichten:
	rFacesAll = [(rInterfaceAll[ii]+rInterfaceAll[ii+1])/2 for ii in range(2*N)]
	# Vektor zur Definition der Interfaces, die eine Element-Konzentration erfahren:
	rInterfacesBias =[]
	ii = 1
	kk = 0
	for jj in range(N+1):
		if jj == 0:
			rInterfacesBias.append(0)
		elif jj == (N):
			rInterfacesBias.append(0)
		else:
			if kk < len(iInterfaceEval):
				if rk[ii] == rk[iInterfaceEval[kk]]:
					rInterfacesBias.append(1)
					rInterfacesBias.append(-1)
					ii = ii + 1
					kk = kk + 1
				else:
					rInterfacesBias.append(0)
					rInterfacesBias.append(0)
					ii = ii + 1
			else:
				rInterfacesBias.append(0)
				rInterfacesBias.append(0)
				ii = ii + 1
	return(rm,rInterfaceEval,rInterfaceAll,rFacesAll,rInterfacesBias)

def circGeometryParameters():
	if boolCircSubModel:
		# Umfangswinkel-Kanten:
		thetaEdges = [angleSubModelStart, angleSubModelStart+angleSubModel, angleSubModelEnd-angleSubModel, angleSubModelEnd]
		# Umfangswinkel-Partitionierungen:
		thetaPartition = [angleSubModelStart+angleSubModel, angleSubModelEnd-angleSubModel]
		# Umfangswinkel-Kanten - Submodel:
		thetaEdgesSubModel = [thetaEdge for thetaEdge in thetaEdges]
	else:
		# Umfangswinkel-Kanten:
		thetaEdges = [angleStart, angleSubModelStart, angleStart+angleEval, angleSubModelEnd, angleStart+angleOpening]
		# Umfangswinkel-Partitionierungen:
		thetaPartition = [angleSubModelStart, angleStart+angleEval, angleSubModelEnd]
		# Umfangswinkel-Kanten - Submodel:
		thetaEdgesSubModel = [thetaEdge for thetaEdge in thetaPartition]
	# Umfangswinkel-Mittenkanten:
	thetaFaces = [(thetaEdges[ii]+thetaEdges[ii+1])/2 for ii in range(len(thetaEdges)-1)]
	thetaFacesSubModel = [(thetaEdgesSubModel[ii]+thetaEdgesSubModel[ii+1])/2 for ii in range(len(thetaEdgesSubModel)-1)]
	return(thetaPartition,thetaEdges,thetaFaces,thetaEdgesSubModel,thetaFacesSubModel)

def axialGeometryParameters():
	# Vektor zur Definition der Breitenpartitionierungen:
	axialPartition = [axialPartitionPercent[ii]*length for ii in range(len(axialPartitionPercent))]
	# Vektor zur Definition aller Edges in Breitenrichtung:
	axialEdges = [0.0]
	for jj in range(len(axialPartition)):
		axialEdges.append(axialPartition[jj])
	axialEdges.append(length)
	# Vektor zur Definition aller Mittenkanten in Breitenrichtung:
	axialFaces = [(axialEdges[ii]+axialEdges[ii+1])/2 for ii in range(len(axialEdges)-1)]
	if boolCircSubModel:
		# Vektor zur Definition aller Edges des SubModels in Breitenrichtung:
		axialEdgesSubModel = [axialPartition[int(len(axialPartition)/2)-1], axialPartition[int(len(axialPartition)/2)], axialPartition[int(len(axialPartition)/2)+1]]
	else:
		# Vektor zur Definition aller Edges des SubModels in Breitenrichtung:
		axialEdgesSubModel = [axialEdge for axialEdge in axialEdges]
	# Vektor zur Definition aller Faces des SubModels in Breitenrichtung:
	axialFacesSubModel = [(axialEdgesSubModel[ii]+axialEdgesSubModel[ii+1])/2 for ii in range(len(axialEdgesSubModel)-1)]
	return(axialPartition,axialEdges,axialFaces,axialEdgesSubModel,axialFacesSubModel)

def surfaceSubModelGeometryParameters(rI):
	rkSurfaceSubModelEdges = [rm[rI-1],rk[rI],rm[rI]]
	rkSurfaceSubModelFaces = [(rkSurfaceSubModelEdges[ii]+rkSurfaceSubModelEdges[ii+1])/2 for ii in range(len(rkSurfaceSubModelEdges)-1)]
	rInterfacesSubModelSurfaceBias = [1,-1]
	return(rkSurfaceSubModelEdges,rkSurfaceSubModelFaces,rInterfacesSubModelSurfaceBias)

rk =[r0]
if boolHybrid:
	for ii in range (N):
		if plyAngle[ii] == 2024:
			rk.append (rk[ii]+dL2024)
		else:
			rk.append (rk[ii]+dL)
else:
	for ii in range (N):
		rk.append (rk[ii]+dL)

rN = rk[-1]

rm,rInterfaceEval,rInterfaceAll,rFacesAll,rInterfacesBias = radialGeometryParameters(N,rk,iInterfaceEval)
thetaPartition,thetaEdges,thetaFaces,thetaEdgesSubModel,thetaFacesSubModel = circGeometryParameters()
axialPartition,axialEdges,axialFaces,axialEdgesSubModel,axialFacesSubModel = axialGeometryParameters()
if boolSurfaceSubModel:
	rkSurfaceSubModelEdges,rkSurfaceSubModelFaces,rInterfacesSubModelSurfaceBias = surfaceSubModelGeometryParameters(surfaceSubModelInterface)


#---------------------------------------------------------------
def sketchCurvedComposite():
	# Skizze/Koerper des zylindrisch gekruemmten Laminats:
	curvedCompositeSketch = compositeModel.ConstrainedSketch(name='Cylindrical_Composite_Sketch', sheetSize=2*rk[-1])
	curvedCompositeSketch.ArcByCenterEnds(center=(0.0, 0.0), point1=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0)), point2=(ct.pol2cart_x(rk[0],  thetaEdges[-1]), ct.pol2cart_y(rk[0], thetaEdges[-1])), direction=COUNTERCLOCKWISE)
	curvedCompositeSketch.ArcByCenterEnds(center=(0.0, 0.0), point1=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0)), point2=(ct.pol2cart_x(rk[-1],  thetaEdges[-1]), ct.pol2cart_y(rk[-1],  thetaEdges[-1])), direction=COUNTERCLOCKWISE)
	
	curvedCompositeSketch.Line(point1=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0)), point2=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0)))
	curvedCompositeSketch.Line(point1=(ct.pol2cart_x(rk[0],  thetaEdges[-1]), ct.pol2cart_y(rk[0],  thetaEdges[-1])), point2=(ct.pol2cart_x(rk[-1],  thetaEdges[-1]), ct.pol2cart_y(rk[-1],  thetaEdges[-1])))
	
	curvedCompositePart = compositeModel.Part(name='Cylindrically_Curved_Composite', dimensionality=THREE_D, type=DEFORMABLE_BODY)
	curvedCompositePart.BaseSolidExtrude(sketch=curvedCompositeSketch, depth=length)
	return(curvedCompositePart)		

curvedCompositePart = sketchCurvedComposite()

#---------------------------------------------------------------
# Assembly:
compositeAssembly = compositeModel.rootAssembly
def instanceCurvedComposite():
	curvedCompositeInstance = compositeAssembly.Instance(name='Curved_Composite_Instance', part=curvedCompositePart, dependent=ON)
	return curvedCompositeInstance

curvedCompositeInstance = instanceCurvedComposite()

#---------------------------------------------------------------
# Zylindrisches Koordinatensystem:
def coordinateSystemCurvedComposite():
	# Erzeugen des zylindrischen Koordinatensystems:
	curvedCompositePart.DatumCsysByThreePoints(coordSysType = CYLINDRICAL, name='Cylindrical_CoordinateSystem',
										point1 = (rk[-1], 0, length/2), 
										point2 = (0, rk[-1], length/2),
										origin= (0, 0, length/2))
	# Zylindrisches Koordinatensystem fuer Randbedingungen:
	curvedCompositeCSCylId = curvedCompositePart.features['Cylindrical_CoordinateSystem'].id
	curvedCompositeCylCoordSys = compositeModel.rootAssembly.instances['Curved_Composite_Instance'].datums[curvedCompositeCSCylId]
	return curvedCompositeCylCoordSys

curvedCompositeCylCoordSys = coordinateSystemCurvedComposite()
#---------------------------------------------------------------
# Partitionierung in radialer Richtung:
# Einteilung des Laminats in physikalische Schichten:
def partitionCurvedCompositeRadial():
	for ii in range (N-1):
		curvedCompositeMakeSketchTransformPlane = curvedCompositePart.faces.findAt((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii],  thetaEdges[-1]/2), length),)
		curvedCompositeMakeSketchTransformUpEdge = curvedCompositePart.edges.findAt((ct.pol2cart_x(rk[0],  thetaEdges[-1]/2), ct.pol2cart_y(rk[0],  thetaEdges[-1]/2), length),)
		curvedCompositeMakeSketchTransform = curvedCompositePart.MakeSketchTransform(sketchPlane=curvedCompositeMakeSketchTransformPlane, sketchPlaneSide=SIDE1, sketchUpEdge=curvedCompositeMakeSketchTransformUpEdge, origin=(0, 0, length))
		curvedCompositeTransformedSketch = compositeModel.ConstrainedSketch(name='Transformed_Sketch', sheetSize=2*rk[-1], transform=curvedCompositeMakeSketchTransform)
		
		curvedCompositePart.projectReferencesOntoSketch(sketch=curvedCompositeTransformedSketch, filter=COPLANAR_EDGES)
		curvedCompositeTransformedSketch.ArcByCenterEnds(center=(0.0,0.0), point1=(ct.pol2cart_x(rk[ii+1], 0), ct.pol2cart_y(rk[ii+1], 0)),point2=(ct.pol2cart_x(rk[ii+1],  thetaEdges[-1]), ct.pol2cart_y(rk[ii+1],  thetaEdges[-1])),direction=COUNTERCLOCKWISE)
		
		curvedCompositePickedFaces = curvedCompositePart.faces.findAt(((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii], thetaEdges[-1]/2), length),),)
		curvedCompositePart.PartitionFaceBySketch(faces=curvedCompositePickedFaces, sketch=curvedCompositeTransformedSketch, sketchUpEdge=curvedCompositeMakeSketchTransformUpEdge)
		
		curvedCompositePickedCells = curvedCompositePart.cells.findAt(((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii],  thetaEdges[-1]/2), length/2),),)
		curvedCompositePickedEdge = (curvedCompositePart.edges.findAt((ct.pol2cart_x(rk[ii+1],  thetaEdges[-1]/2), ct.pol2cart_y(rk[ii+1],  thetaEdges[-1]/2), length),),)
		curvedCompositePart.PartitionCellByExtrudeEdge(line=curvedCompositePart.edges.findAt((ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0), length/4)), cells=curvedCompositePickedCells, edges=curvedCompositePickedEdge, sense=REVERSE)
	
	# Zusaetzliche Einteilung der physikalischen Schichten hinsichtlich ihrer mittleren Radien:
	for ii in range (N):
		curvedCompositeMakeSketchTransformPlane = curvedCompositePart.faces.findAt((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii],  thetaEdges[-1]/2), length),)
		curvedCompositeMakeSketchTransformUpEdge = curvedCompositePart.edges.findAt((ct.pol2cart_x(rk[0],  thetaEdges[-1]/2), ct.pol2cart_y(rk[0],  thetaEdges[-1]/2), length),)
		curvedCompositeMakeSketchTransform = curvedCompositePart.MakeSketchTransform(sketchPlane=curvedCompositeMakeSketchTransformPlane, sketchPlaneSide=SIDE1, sketchUpEdge=curvedCompositeMakeSketchTransformUpEdge, origin=(0, 0, length))
		curvedCompositeTransformedSketch = compositeModel.ConstrainedSketch(name='Transformed_Sketch', sheetSize=2*rk[-1], transform=curvedCompositeMakeSketchTransform)
		
		curvedCompositePart.projectReferencesOntoSketch(sketch=curvedCompositeTransformedSketch, filter=COPLANAR_EDGES)
		curvedCompositeTransformedSketch.ArcByCenterEnds(center=(0.0,0.0), point1=(ct.pol2cart_x(rm[ii], 0), ct.pol2cart_y(rm[ii], 0)),point2=(ct.pol2cart_x(rm[ii],  thetaEdges[-1]), ct.pol2cart_y(rm[ii],  thetaEdges[-1])),direction=COUNTERCLOCKWISE)
		
		curvedCompositePickedFaces = curvedCompositePart.faces.findAt(((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii], thetaEdges[-1]/2), length),),)
		curvedCompositePart.PartitionFaceBySketch(faces=curvedCompositePickedFaces, sketch=curvedCompositeTransformedSketch, sketchUpEdge=curvedCompositeMakeSketchTransformUpEdge)
		
		curvedCompositePickedCells = curvedCompositePart.cells.findAt(((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii],  thetaEdges[-1]/2), length/2),),)
		curvedCompositePickedEdge = (curvedCompositePart.edges.findAt((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii],  thetaEdges[-1]/2), length),),)
		curvedCompositePart.PartitionCellByExtrudeEdge(line=curvedCompositePart.edges.findAt((ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0), length/4)), cells=curvedCompositePickedCells, edges=curvedCompositePickedEdge, sense=REVERSE)

def partitionCurvedCompositeCirc():
	compositePartitionPlaneLaufVar = 1
	# Partitionierung in tagentialer Richtung:
	for ii in range(len(thetaPartition)):
		compositePartitionDatumPointLaufVar = len(curvedCompositePart.datums)
		
		curvedCompositePart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[0], thetaPartition[ii]), ct.pol2cart_y(rk[0], thetaPartition[ii]), length))
		curvedCompositePart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1], thetaPartition[ii]), ct.pol2cart_y(rk[-1], thetaPartition[ii]), length))
		curvedCompositePart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1], thetaPartition[ii]), ct.pol2cart_y(rk[-1], thetaPartition[ii]), 0))
		
		compositePartDatumsKeys = curvedCompositePart.datums.keys()
		compositePartDatumsKeys.sort()
		
		compositePartitionCell = curvedCompositePart.cells
		curvedCompositePart.DatumPlaneByThreePoints(point1=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar]], 
																		point2=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar+1]],
																		point3=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar+2]])
		compositePartitionPlaneId = curvedCompositePart.features['Datum plane-'+ str(compositePartitionPlaneLaufVar)].id
		curvedCompositePart.PartitionCellByDatumPlane(cells=compositePartitionCell, datumPlane=curvedCompositePart.datums[compositePartitionPlaneId])
		
		compositePartitionPlaneLaufVar += 1

def partitionCurvedCompositeAxial():
	compositePartitionPlaneLaufVar = len(thetaPartition) + 1
	# Partitionierung in Breitenrichtung:
	for kk in range(len(axialPartition)):
		compositePartitionDatumPointLaufVar = len(curvedCompositePart.datums)
		
		curvedCompositePart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0), axialPartition[kk]))
		curvedCompositePart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0), axialPartition[kk]))
		curvedCompositePart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1],  thetaEdges[-1]), ct.pol2cart_y(rk[-1],  thetaEdges[-1]), axialPartition[kk]))
		
		compositePartDatumsKeys = curvedCompositePart.datums.keys()
		compositePartDatumsKeys.sort()
		
		compositePartitionCell = curvedCompositePart.cells
		curvedCompositePart.DatumPlaneByThreePoints(point1=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar]], 
																		point2=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar+1]],
																		point3=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar+2]])
		compositePartitionPlaneId = curvedCompositePart.features['Datum plane-'+ str(compositePartitionPlaneLaufVar)].id
		curvedCompositePart.PartitionCellByDatumPlane(cells=compositePartitionCell, datumPlane=curvedCompositePart.datums[compositePartitionPlaneId])
		
		compositePartitionPlaneLaufVar += 1

partitionCurvedCompositeRadial()
partitionCurvedCompositeCirc()
partitionCurvedCompositeAxial()
#---------------------------------------------------------------
def materialMetal():
	# Materialdefinition ueber Ingenieurkonstanten:
	aluminumMaterial = compositeModel.Material(metalMaterialName) 
	aluminumMaterial.Elastic(table=((E2024, Nu2024), ))
	compositeModel.HomogeneousSolidSection(name='Metal_Section', material=metalMaterialName, thickness=None)

def materialComposite():
	compositeMaterial = compositeModel.Material(compositeMaterialName)
	compositeMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=((E1, E2, E3, Nu12, Nu13, Nu23, G23, G13, G12), ))
	compositeMaterial.Expansion(type=ORTHOTROPIC, table=((alpha11, alpha22, alpha33), ))
	compositeModel.HomogeneousSolidSection(name='Composite_Section', material=compositeMaterialName, thickness=None)

materialComposite()
if boolHybrid:
	materialMetal()

def sectionCurvedComposite():
	# Section-Zuweisung bzw. Definition der Faserorientierungswinkel der einzelnen physikalischen Schichten:
	jj = 0
	for ii in range(N):
		curvedCompositeAngleInc = 0
		for mm in range(len(thetaFaces)):
			for kk in range(len(axialEdges)-1):
				curvedCompositeCells = curvedCompositePart.cells.findAt(((ct.pol2cart_x(rFacesAll[jj], thetaFaces[mm]), ct.pol2cart_y(rFacesAll[jj], thetaFaces[mm]), (axialEdges[kk]+axialEdges[kk+1])/2),),
															((ct.pol2cart_x(rFacesAll[jj+1], thetaFaces[mm]), ct.pol2cart_y(rFacesAll[jj+1], thetaFaces[mm]), (axialEdges[kk]+axialEdges[kk+1])/2),))
				
				curvedCompositeRegion = regionToolset.Region(cells=curvedCompositeCells)
				if plyAngle[ii] == 2024:
					curvedCompositePart.SectionAssignment(region=curvedCompositeRegion, sectionName='Metal_Section', offsetType=MIDDLE_SURFACE)
				else:
					curvedCompositePart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section', offsetType=MIDDLE_SURFACE)
					curvedCompositePart.MaterialOrientation(region=curvedCompositeRegion, orientationType=SYSTEM, localCsys=mdb.models[modelName].parts['Cylindrically_Curved_Composite'].datums[2], axis=AXIS_1, additionalRotationType=ROTATION_ANGLE, angle=plyAngle[ii], stackDirection=STACK_1)
			curvedCompositeAngleInc = curvedCompositeAngleInc + 1
		jj = jj + 2

sectionCurvedComposite()
#---------------------------------------------------------------
# Vernetzung in radialer Richtung:
mR,mT,mZ = 5,50,10
mRRatio,mTRatio,mZRatio = 1,1,1
def meshCurvedCompositeGlobal():
	for ii in range(len(rFacesAll)):
		for jj in range(len(thetaEdges)):
			for kk in range(len(axialEdges)):
				curvedCompositeMeshEdge = curvedCompositePart.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj]), ct.pol2cart_y(rFacesAll[ii], thetaEdges[jj]), axialEdges[kk]),))
				curvedCompositePart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRRatio, number=mR, constraint=FINER)
	# Vernetzung in tangentialer Richtung:
	for ii in range (len(rInterfaceAll)):
		for jj in range(len(thetaFaces)):
			for kk in range(len(axialEdges)):
				curvedCompositeMeshEdge = curvedCompositePart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaFaces[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaFaces[jj]), axialEdges[kk]),))
				if int(round(mT*(thetaEdges[jj+1]-thetaEdges[jj])/thetaEdges[-1])) < 1:
					curvedCompositePart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mTRatio, number=1, constraint=FINER)
				else:
					curvedCompositePart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mTRatio, number=int(round(mT*(thetaEdges[jj+1]-thetaEdges[jj])/thetaEdges[-1])), constraint=FINER)
	# Vernetzung in Breitenrichtung:
	for ii in range (len(rInterfaceAll)):
		for jj in range(len(thetaEdges)):
			for kk in range(len(axialFaces)):
				curvedCompositeMeshEdge = curvedCompositePart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaEdges[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdges[jj]), axialFaces[kk]),))
				if int(round(mZ*(axialEdges[kk+1]-axialEdges[kk])/length)) < 1:
					curvedCompositePart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mZRatio, number=1, constraint=FINER)
				else:
					curvedCompositePart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mZRatio, number=int(round(mZ*(axialEdges[kk+1]-axialEdges[kk])/length)), constraint=FINER)
	# Element types: C3D8, C3D8R, C3D20, C3D20R
	compositeElementType = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
	curvedCompositePart.setMeshControls(regions=curvedCompositePart.cells, elemShape=HEX, technique=STRUCTURED)
	curvedCompositePart.setElementType(regions=(curvedCompositePart.cells[:], ), elemTypes=(compositeElementType, ))
	curvedCompositePart.seedPart(size=0.05, deviationFactor=0.1, constraint=FINER)
	curvedCompositePart.generateMesh()

meshCurvedCompositeGlobal()
#---------------------------------------------------------------
# Step:
def stepMechanicalLoads():
	stepPrevious,step,stepDescription = 'Initial','Mechanical_Loading','Introduce mechanical loadings'
	compositeModel.StaticStep(name=step, previous=stepPrevious, description=stepDescription, nlgeom=OFF, initialInc=1.0)
	return (stepPrevious,step)

def stepHygrothermalLoads():
	stepPrevious,step,stepDescription = 'Initial','Hygrothermal_Loading','Introduce hygrothermal loadings'
	compositeModel.StaticStep(name=step, previous=stepPrevious, description=stepDescription, nlgeom=OFF, initialInc=1.0)
	return (stepPrevious,step)

if boolBendingMoment or boolRadialForce or boolCircForce or boolPressure or boolCylindricalPressure:
	stepPrevious,step = stepMechanicalLoads()
elif boolTempDif or boolMoistDif:
	stepPrevious,step = stepHygrothermalLoads()

#---------------------------------------------------------------
# Randbedingungen:
def boundaryConditionBendingMoment():
	pass

def boundaryConditionEdgeLoad():
	# Sperren der radialen Bewegung:
	curvedCompositeStringFixRadialEdge = ''
	for kk in range(len(axialFaces)):
		curvedCompositeStringFixRadialEdge = curvedCompositeStringFixRadialEdge + str(((ct.pol2cart_x(r0,  thetaEdges[-1]), ct.pol2cart_y(r0, thetaEdges[-1]), axialFaces[kk]),),) + ','
		curvedCompositeFixRadialEdgesExec = 'curvedCompositeFixRadialEdge = curvedCompositeInstance.edges.findAt(' + curvedCompositeStringFixRadialEdge + ')'
		exec(curvedCompositeFixRadialEdgesExec)
	curvedCompositeFixRadialSet = compositeAssembly.Set(edges=curvedCompositeFixRadialEdge, name='Set_Fix_U_Interface_'+ str(0))
	compositeModel.DisplacementBC(name='Fix_U_Interface_'+ str(0), createStepName=stepPrevious, region=curvedCompositeFixRadialSet, localCsys=curvedCompositeCylCoordSys, u1=0.0)
	# Sperren der tangentialen Bewegung
	for ii in range (N):
		curvedCompositeStringFixCircFaces = ''
		for kk in range(len(axialFaces)):
			curvedCompositeStringFixCircFaces = curvedCompositeStringFixCircFaces + str(((ct.pol2cart_x((rm[ii]+rk[ii])/2, thetaEdges[-1]), ct.pol2cart_y((rm[ii]+rk[ii])/2, thetaEdges[-1]), axialFaces[kk]),),) + ','
			curvedCompositeStringFixCircFaces = curvedCompositeStringFixCircFaces + str(((ct.pol2cart_x((rm[ii]+rk[ii+1])/2, thetaEdges[-1]), ct.pol2cart_y((rm[ii]+rk[ii+1])/2, thetaEdges[-1]), axialFaces[kk]),),) + ','
			curvedCompositeFixCircFacesExec = 'curvedCompositeFixCircFaces = curvedCompositeInstance.faces.findAt(' + curvedCompositeStringFixCircFaces + ')'
			exec(curvedCompositeFixCircFacesExec)
		
		curvedCompositeFixCircSet = compositeAssembly.Set(faces=curvedCompositeFixCircFaces, name='Set_Fix_V_Layer_'+ str(ii+1))
		compositeModel.DisplacementBC(name='Fix_V_Layer_'+ str(ii+1), createStepName=stepPrevious, region=curvedCompositeFixCircSet, localCsys=curvedCompositeCylCoordSys, u2 = 0.0)
	# Sperren der Bewegung in Axialrichtung:
	for ii in range (N):
		curvedCompositeFixAxialPoint = curvedCompositeInstance.edges.findAt(((ct.pol2cart_x((rm[ii]+rk[ii])/2,  thetaEdges[-1]), ct.pol2cart_y((rm[ii]+rk[ii])/2,  thetaEdges[-1]), length/2),),((ct.pol2cart_x((rm[ii]+rk[ii+1])/2,  thetaEdges[-1]), ct.pol2cart_y((rm[ii]+rk[ii+1])/2,  thetaEdges[-1]), length/2),))
		curvedCompositeFixAxialSet = compositeAssembly.Set(edges=curvedCompositeFixAxialPoint, name='Set_Fix_W_Layer_'+ str(ii+1))
		compositeModel.DisplacementBC(name='Fix_W_Layer_'+ str(ii+1), createStepName=stepPrevious, region=curvedCompositeFixAxialSet, localCsys=curvedCompositeCylCoordSys, u3=0.0)

def boundaryConditionSurfaceLoad():
	for jj,ll in enumerate([thetaEdges[0],thetaEdges[-1]]):
		curvedCompositeStringFixCircFaces = ''
		for ii in range(len(rFacesAll)):
			for kk in range(len(axialFaces)):
				curvedCompositeStringFixCircFaces = curvedCompositeStringFixCircFaces + str(((ct.pol2cart_x(rFacesAll[ii], ll), ct.pol2cart_y(rFacesAll[ii], ll), axialFaces[kk]),),) + ','
				curvedCompositeFixCircFacesExec = 'curvedCompositeFixCircFaces = curvedCompositeInstance.faces.findAt(' + curvedCompositeStringFixCircFaces + ')'
				exec(curvedCompositeFixCircFacesExec)
		curvedCompositeFixCircSet = compositeAssembly.Set(faces=curvedCompositeFixCircFaces, name='Set_Fix_U_SS_'+ str(jj))
		compositeModel.DisplacementBC(name='Fix_U_SS_'+ str(jj), createStepName=stepPrevious, region=curvedCompositeFixCircSet, localCsys=curvedCompositeCylCoordSys, u1=0.0)

def boundaryConditionHygrothermal():
	pass

if boolBendingMoment and not any([boolRadialForce,boolCircForce]):
	boundaryConditionBendingMoment()
elif any([boolBendingMoment,boolRadialForce,boolCircForce]):
	boundaryConditionEdgeLoad()
elif any([boolPressure,boolCylindricalPressure]):
	boundaryConditionSurfaceLoad()
elif any([boolTempDif,boolMoistDif]):
	boundaryConditionHygrothermal()
else:
	pass

#---------------------------------------------------------------
# Lastfall:
def loadsCreateKinCoup(jj):
	# Erstellen eines Referenzpunktes und Fixierung dessen:
	compositeAssembly.ReferencePoint(point=(ct.pol2cart_x((rk[0]+rk[-1])/2, jj*angleOpening), ct.pol2cart_y((rk[0]+rk[-1])/2, jj*angleOpening), length/2))
	curvedCompositeLoadRefPoint = (compositeAssembly.referencePoints.findAt((ct.pol2cart_x((rk[0]+rk[-1])/2, jj*angleOpening), ct.pol2cart_y((rk[0]+rk[-1])/2, jj*angleOpening), length/2),),)
	curvedCompositeLoadRefPointRegion = regionToolset.Region(referencePoints=curvedCompositeLoadRefPoint)
	curvedCompositeFixRefPointSet = compositeAssembly.Set(region=curvedCompositeLoadRefPointRegion, name='Set_Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)))
	compositeModel.DisplacementBC(name='Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeFixRefPointSet, localCsys=curvedCompositeCylCoordSys, u3=0.0, ur1=0.0, ur2=0.0)
	# Kinematic Coupling:
	curvedCompositeStringKinCoupFaces =''
	for ii in range (N):
		for kk in range(len(axialFaces)):
			curvedCompositeStringKinCoupFaces = curvedCompositeStringKinCoupFaces + str(((ct.pol2cart_x((rm[ii]+rk[ii])/2, jj*angleOpening), ct.pol2cart_y((rm[ii]+rk[ii])/2, jj*angleOpening), axialFaces[kk]),),) + ','
			curvedCompositeStringKinCoupFaces = curvedCompositeStringKinCoupFaces + str(((ct.pol2cart_x((rm[ii]+rk[ii+1])/2, jj*angleOpening), ct.pol2cart_y((rm[ii]+rk[ii+1])/2, jj*angleOpening), axialFaces[kk]),),) + ','
	curvedCompositeKinCoupFacesExec = 'curvedCompositeKinCoupFaces = curvedCompositeInstance.faces.findAt(' + curvedCompositeStringKinCoupFaces + ')'
	exec(curvedCompositeKinCoupFacesExec)
	curvedCompositeKinCoupRegion = regionToolset.Region(side1Faces=curvedCompositeKinCoupFaces)
	if boolBendingMoment:
		compositeModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion,influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=ON)
	else:
		compositeModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion,influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=ON, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=ON)
	compositeModel.constraints['Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi))].setValues(localCsys=curvedCompositeCylCoordSys)
	return curvedCompositeLoadRefPointRegion

def loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion):
	# Aeussere Last - Moment:
	compositeModel.Moment(name='Load_BendingMoment_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeLoadRefPointRegion, cm3=((-1)**(jj))*length*bendingMoment, localCsys=curvedCompositeCylCoordSys, distributionType=UNIFORM)

def loadsRadialForce():
	# Aeussere Last - Querkraft:
	compositeModel.ConcentratedForce(name='Load_RadialForce', createStepName=step, region=curvedCompositeLoadRefPointRegion, cf1=length*radialForce, localCsys=curvedCompositeCylCoordSys, distributionType=UNIFORM)

def loadsCircForce():
	# Aeussere Last - Normalkraft:
	compositeModel.ConcentratedForce(name='Load_CircForce', createStepName=step, region=curvedCompositeLoadRefPointRegion, cf2=-length*circForce, localCsys=curvedCompositeCylCoordSys, distributionType=UNIFORM)

def loadsSurfaceLoad():
	if OuterPressure != 0.0 or cylindricalOuterPressure != 0.0:
		# Aeussere Last - Aussendruck:
		curvedCompositeStringFaces = ''
		for jj in range(len(thetaFaces)):
			for kk in range(len(axialFaces)):
				curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[-1], thetaFaces[jj]), ct.pol2cart_y(rk[-1], thetaFaces[jj]), axialFaces[kk]),),) + ','
		
		curvedCompositeStringFacesExec = 'curvedCompositeOuterPressureFaces = curvedCompositeInstance.faces.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		
		if boolCylindricalPressure:
			compositeModel.ExpressionField(name='Load_OuterPressure_CylindricalBending', localCsys=curvedCompositeCylCoordSys, description='',expression='sin(pi*Th/angleOpening)')
			compositeModel.Pressure(name='Load_OuterPressure_CylindricalBending', createStepName=step, region=regionToolset.Region(side1Faces=curvedCompositeOuterPressureFaces), magnitude=cylindricalOuterPressure, distributionType=FIELD, field='Load_OuterPressure_CylindricalBending', amplitude=UNSET)
		else:
			compositeModel.Pressure(name='Load_OuterPressure', createStepName=step, region=regionToolset.Region(side1Faces=curvedCompositeOuterPressureFaces), magnitude=OuterPressure, distributionType=UNIFORM)
	if InnerPressure != 0.0 or cylindricalInnerPressure != 0.0:
		# Aeussere Last - Innendruck:
		curvedCompositeStringFaces = ''
		for jj in range(len(thetaFaces)):
			for kk in range(len(axialFaces)):
				curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[0], thetaFaces[jj]), ct.pol2cart_y(rk[0], thetaFaces[jj]), axialFaces[kk]),),) + ','
		
		curvedCompositeStringFacesExec = 'curvedCompositeInnerPressureFaces = curvedCompositeInstance.faces.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		
		if boolCylindricalPressure:
			compositeModel.ExpressionField(name='Load_InnerPressure_CylindricalBending', localCsys=curvedCompositeCylCoordSys, description='',expression='sin(pi*Th/angleOpening)')
			compositeModel.Pressure(name='Load_InnerPressure_CylindricalBending', createStepName=step, region=regionToolset.Region(side1Faces=curvedCompositeInnerPressureFaces), magnitude=cylindricalInnerPressure, distributionType=FIELD, field='Load_InnerPressure_CylindricalBending', amplitude=UNSET)
		else:
			compositeModel.Pressure(name='Load_InnerPressure', createStepName=step, region=regionToolset.Region(side1Faces=curvedCompositeInnerPressureFaces), magnitude=InnerPressure, distributionType=UNIFORM)

def loadsHygrothermal():
	# Definition der Ausgangstemperatur:
	compositeModel.Temperature(name='Temp_Initial', createStepName=stepPrevious, region=regionToolset.Region(cells=curvedCompositeInstance.cells), distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0,))
	compositeModel.Temperature(name='Temperature_Difference', createStepName=step, region=regionToolset.Region(cells=curvedCompositeInstance.cells), distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(tempDif,))

if boolBendingMoment and not any([boolRadialForce,boolCircForce]):
	noKinCoup = 2
	for jj in range(noKinCoup):
		curvedCompositeLoadRefPointRegion = loadsCreateKinCoup(jj)
		loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion)
elif any([boolBendingMoment,boolRadialForce,boolCircForce]):
	noKinCoup = 1
	for jj in range(noKinCoup):
		curvedCompositeLoadRefPointRegion = loadsCreateKinCoup(jj)
		if boolBendingMoment:
			loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion)
	if boolRadialForce:
		loadsRadialForce()
	if boolCircForce:
		loadsCircForce()
elif any([boolPressure,boolCylindricalPressure]):
	noKinCoup = 2
	for jj in range(noKinCoup):
		curvedCompositeLoadRefPointRegion = loadsCreateKinCoup(jj)
	loadsSurfaceLoad()
elif any([boolTempDif,boolMoistDif]):
	loadsHygrothermal()
else:
	pass

# ---------------------------------------------------------------
# Field Output:
compositeModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'E', 'COORD'))

# ---------------------------------------------------------------
# History Output

#-------------------------------------------------------------------
# Job
# Erstellen des Jobs:
jobName = modelName + '_Job'
mdb.Job(name=jobName, model=modelName,
		description='Run FE-analysis', parallelizationMethodExplicit=DOMAIN, numDomains=12, numCpus=12, memory=85,echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)

# Job starten:
mdb.jobs[jobName].submit(consistencyChecking=OFF)
mdb.jobs[jobName].waitForCompletion()

#-------------------------------------------------------------------
#Sub-Modellierung
#-------------------------------------------------------------------
# Kopie des Main-Modells:
compositeSubModel = mdb.Model(name=modelName+'_SubModel', objectToCopy=compositeModel)
curvedCompositeSubModelPart = mdb.models[modelName+'_SubModel'].parts['Cylindrically_Curved_Composite']
session.viewports['Viewport: 1'].setValues(displayedObject=curvedCompositeSubModelPart)

compositeSubModelAssembly = compositeSubModel.rootAssembly
compositeSubModelInstance = compositeSubModelAssembly.Instance(name='Curved_Composite_Instance', part=curvedCompositeSubModelPart, dependent=ON)

#---------------------------------------------------------------
# Remove all faces not needed:
def subModelAxialRemoveGeometry():
	for ii in range (N):
		compositeEvalAngleInc = np.pi/180
		for kk in range(len(axialEdges)):
			compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
			((ct.pol2cart_x((rk[ii]+rm[ii])/2, compositeEvalAngleInc), ct.pol2cart_y((rk[ii]+rm[ii])/2, compositeEvalAngleInc), axialEdges[kk] ),),
			((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, compositeEvalAngleInc), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, compositeEvalAngleInc), axialEdges[kk] ),),
			((ct.pol2cart_x((rk[ii]+rm[ii])/2,  thetaEdges[-1]-compositeEvalAngleInc), ct.pol2cart_y((rk[ii]+rm[ii])/2,  thetaEdges[-1]-compositeEvalAngleInc), axialEdges[kk] ),),
			((ct.pol2cart_x((rk[ii+1]+rm[ii])/2,  thetaEdges[-1]-compositeEvalAngleInc), ct.pol2cart_y((rk[ii+1]+rm[ii])/2,  thetaEdges[-1]-compositeEvalAngleInc), axialEdges[kk] ),))
			curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		
		compositeEvalAngleInc = 0.0
		for kk in range(len(axialEdges)-1):
			compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
			((ct.pol2cart_x((rk[ii]+rm[ii])/2, compositeEvalAngleInc), ct.pol2cart_y((rk[ii]+rm[ii])/2, compositeEvalAngleInc), (axialEdges[kk]+axialEdges[kk+1])/2 ),),
			((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, compositeEvalAngleInc), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, compositeEvalAngleInc), (axialEdges[kk]+axialEdges[kk+1])/2 ),),
			((ct.pol2cart_x((rk[ii]+rm[ii])/2,  thetaEdges[-1]), ct.pol2cart_y((rk[ii]+rm[ii])/2,  thetaEdges[-1]), (axialEdges[kk]+axialEdges[kk+1])/2 ),),
			((ct.pol2cart_x((rk[ii+1]+rm[ii])/2,  thetaEdges[-1]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2,  thetaEdges[-1]), (axialEdges[kk]+axialEdges[kk+1])/2 ),))
			curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		
		compositeEvalAngleInc = np.pi/180
		for kk in range(len(axialEdges)-1):
			compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
			((ct.pol2cart_x(rk[ii], compositeEvalAngleInc), ct.pol2cart_y(rk[ii], compositeEvalAngleInc), (axialEdges[kk]+axialEdges[kk+1])/2 ),),
			((ct.pol2cart_x(rm[ii], compositeEvalAngleInc), ct.pol2cart_y(rm[ii], compositeEvalAngleInc), (axialEdges[kk]+axialEdges[kk+1])/2 ),),
			((ct.pol2cart_x(rk[ii],  thetaEdges[-1]-compositeEvalAngleInc), ct.pol2cart_y(rk[ii],  thetaEdges[-1]-compositeEvalAngleInc), (axialEdges[kk]+axialEdges[kk+1])/2 ),),
			((ct.pol2cart_x(rm[ii],  thetaEdges[-1]-compositeEvalAngleInc), ct.pol2cart_y(rm[ii],  thetaEdges[-1]-compositeEvalAngleInc), (axialEdges[kk]+axialEdges[kk+1])/2 ),))
			curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		
		if ii == (N-1):
			for kk in range(len(axialEdges)-1):
				compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
				((ct.pol2cart_x(rk[ii+1], compositeEvalAngleInc), ct.pol2cart_y(rk[ii+1], compositeEvalAngleInc), (axialEdges[kk]+axialEdges[kk+1])/2 ),),
				((ct.pol2cart_x(rk[ii+1],  thetaEdges[-1]-compositeEvalAngleInc), ct.pol2cart_y(rk[ii+1],  thetaEdges[-1]-compositeEvalAngleInc), (axialEdges[kk]+axialEdges[kk+1])/2 ),))
				curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
	# Connect the results of the submodel and main model:
	compositeSubModel.setValues(globalJob=jobName)

def subModelCircRemoveGeometry():
	for ii in range(N):
		for jj in range(len(thetaEdges)):
			compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
			((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaEdges[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2, thetaEdges[jj]), axialFaces[0] ),),
			((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaEdges[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, thetaEdges[jj]), axialFaces[0] ),),
			((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaEdges[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2,  thetaEdges[jj]), axialFaces[-1] ),),
			((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaEdges[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2,  thetaEdges[jj]), axialFaces[-1] ),))
			curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		for jj in range(len(thetaFaces)):
			compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
			((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2, thetaFaces[jj]), axialEdges[0] ),),
			((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), axialEdges[0] ),),
			((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2,  thetaFaces[jj]), axialEdges[-1] ),),
			((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2,  thetaFaces[jj]), axialEdges[-1] ),))
			curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		for jj in range(len(thetaFaces)):
			compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
			((ct.pol2cart_x(rk[ii], thetaFaces[jj]), ct.pol2cart_y(rk[ii], thetaFaces[jj]), axialFaces[0] ),),
			((ct.pol2cart_x(rm[ii], thetaFaces[jj]), ct.pol2cart_y(rm[ii], thetaFaces[jj]), axialFaces[0] ),),
			((ct.pol2cart_x(rk[ii], thetaFaces[jj]), ct.pol2cart_y(rk[ii], thetaFaces[jj]), axialFaces[-1] ),),
			((ct.pol2cart_x(rm[ii], thetaFaces[jj]), ct.pol2cart_y(rm[ii], thetaFaces[jj]), axialFaces[-1] ),))
			curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		if ii == (N-1):
			for jj in range(len(thetaFaces)):
				compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
				((ct.pol2cart_x(rk[ii+1], thetaFaces[jj]), ct.pol2cart_y(rk[ii+1], thetaFaces[jj]), axialFaces[0] ),),
				((ct.pol2cart_x(rk[ii+1], thetaFaces[jj]), ct.pol2cart_y(rk[ii+1], thetaFaces[jj]), axialFaces[-1] ),))
				curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
	# Connect the results of the submodel and main model:
	compositeSubModel.setValues(globalJob=jobName)

def subModelCheckCorner(quadrant):
	if quadrant == 1:
		notRemoveAxialFace,notRemoveThetaFace = axialFaces[-1],thetaFaces[0]
		notRemoveAxialEdges,notRemoveThetaEdges = axialEdges[-2:],thetaEdges[:2]
	elif quadrant == 2:
		notRemoveAxialFace,notRemoveThetaFace = axialFaces[-1],thetaFaces[-1]
		notRemoveAxialEdges,notRemoveThetaEdges = axialEdges[-2:],thetaEdges[-2:]
	elif quadrant == 3:
		notRemoveAxialFace,notRemoveThetaFace = axialFaces[0],thetaFaces[-1]
		notRemoveAxialEdges,notRemoveThetaEdges = axialEdges[:2],thetaEdges[-2:]
	elif quadrant == 4:
		notRemoveAxialFace,notRemoveThetaFace = axialFaces[0],thetaFaces[0]
		notRemoveAxialEdges,notRemoveThetaEdges = axialEdges[:2],thetaEdges[:2]
	else:
		raise ValueError('Please check the input value of the corner you want to investigate in detail!')
	return(notRemoveAxialFace,notRemoveThetaFace,notRemoveAxialEdges,notRemoveThetaEdges)

def subModelCornerRemoveGeometry(quadrant=1):
	notRemoveAxialFace,notRemoveThetaFace,notRemoveAxialEdges,notRemoveThetaEdges = subModelCheckCorner(quadrant)
	for ii in range(N):
		for kk in range(len(axialFaces)):
			for jj in range(len(thetaEdges)):
				if (axialFaces[kk] == notRemoveAxialFace and thetaEdges[jj] in notRemoveThetaEdges):
					pass
				else:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaEdges[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2, thetaEdges[jj]), axialFaces[kk] ),),
					((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaEdges[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, thetaEdges[jj]), axialFaces[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		for kk in range(len(axialEdges)):
			for jj in range(len(thetaFaces)):
				if (thetaFaces[jj] == notRemoveThetaFace and axialEdges[kk] in notRemoveAxialEdges):
					pass
				else:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2, thetaFaces[jj]), axialEdges[kk] ),),
					((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), axialEdges[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		for kk in range(len(axialFaces)):
			for jj in range(len(thetaFaces)):
				if (axialFaces[kk] == notRemoveAxialFace and thetaFaces[jj] == notRemoveThetaFace):
					pass
				else:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x(rk[ii], thetaFaces[jj]), ct.pol2cart_y(rk[ii], thetaFaces[jj]), axialFaces[kk] ),),
					((ct.pol2cart_x(rm[ii], thetaFaces[jj]), ct.pol2cart_y(rm[ii], thetaFaces[jj]), axialFaces[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
		if ii == (N-1):
			for kk in range(len(axialFaces)):
				for jj in range(len(thetaFaces)):
					if (axialFaces[kk] == notRemoveAxialFace and thetaFaces[jj] == notRemoveThetaFace):
						pass
					else:
						compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
						((ct.pol2cart_x(rk[ii+1], thetaFaces[jj]), ct.pol2cart_y(rk[ii+1], thetaFaces[jj]), axialFaces[kk] ),))
						curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
	# Connect the results of the submodel and main model:
	compositeSubModel.setValues(globalJob=jobName)

def subModelSurfaceRemoveGeometry():
	for ii in range(N):
		for jj in range(len(thetaEdges)):
			for kk in range(len(axialFaces)):
				if (rk[ii]+rm[ii])/2 not in rkSurfaceSubModelFaces:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaEdges[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2, thetaEdges[jj]), axialFaces[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
				else:
					pass
				if (rk[ii+1]+rm[ii])/2 not in rkSurfaceSubModelFaces:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaEdges[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, thetaEdges[jj]), axialFaces[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
				else:
					pass
		for jj in range(len(thetaFaces)):
			for kk in range(len(axialEdges)):
				if (rk[ii]+rm[ii])/2 not in rkSurfaceSubModelFaces:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2, thetaFaces[jj]), axialEdges[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
				else:
					pass
				if (rk[ii+1]+rm[ii])/2 not in rkSurfaceSubModelFaces:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), axialEdges[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
				else:
					pass
		for jj in range(len(thetaFaces)):
			for kk in range(len(axialFaces)):
				if rk[ii] not in rkSurfaceSubModelEdges:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x(rk[ii], thetaFaces[jj]), ct.pol2cart_y(rk[ii], thetaFaces[jj]), axialFaces[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
				else:
					pass
				if rm[ii] not in rkSurfaceSubModelEdges:
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x(rm[ii], thetaFaces[jj]), ct.pol2cart_y(rm[ii], thetaFaces[jj]), axialFaces[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
				else:
					pass
		if ii == (N-1):
			for jj in range(len(thetaFaces)):
				for kk in range(len(axialFaces)):
					compositeFacesDelete = curvedCompositeSubModelPart.faces.findAt(
					((ct.pol2cart_x(rk[ii+1], thetaFaces[jj]), ct.pol2cart_y(rk[ii+1], thetaFaces[jj]), axialFaces[kk] ),))
					curvedCompositeSubModelPart.RemoveFaces(faceList = compositeFacesDelete, deleteCells=False)
	# Connect the results of the submodel and main model:
	compositeSubModel.setValues(globalJob=jobName)

if boolAxialSubModel:
	subModelAxialRemoveGeometry()
elif boolCircSubModel and not boolCornerSubModel and not boolSurfaceSubModel:
	subModelCircRemoveGeometry()
elif boolCornerSubModel and not boolSurfaceSubModel:
	subModelCornerRemoveGeometry(cornerSubModelQuadrant)
elif boolSurfaceSubModel:
	subModelSurfaceRemoveGeometry()
else:
	pass

#---------------------------------------------------------------
def subModelAxialRemoveBCEdgeLoad():
	compositeSubModel.boundaryConditions.delete(('Fix_U_Interface_'+ str(0),))
	for ii in range (N):
		compositeSubModel.boundaryConditions.delete(('Fix_V_Layer_'+str(ii+1),))
		compositeSubModel.boundaryConditions.delete(('Fix_W_Layer_'+str(ii+1),))

def subModelAxialRemoveBCSurfaceLoad():
	for jj in range(2):
		compositeSubModel.boundaryConditions.delete(('Fix_U_SS_'+ str(jj),))

def subModelAxialRemoveLoadsBendingMoment(noKinCoup):
	for jj in range(noKinCoup):
		del compositeSubModel.loads['Load_BendingMoment_'+ str(jj*int(180*angleOpening/np.pi))]

def subModelAxialRemoveLoadsEdgeLoad(noKinCoup):
	if boolBendingMoment:
		subModelAxialRemoveLoadsBendingMoment(noKinCoup)
	if boolRadialForce:
		del compositeSubModel.loads['Load_RadialForce']
	if boolCircForce:
		del compositeSubModel.loads['Load_CircForce']

def subModelAxialRemoveLoadsSurfaceLoad():
	if boolPressure:
		if OuterPressure != 0.0:
			del compositeSubModel.loads['Load_OuterPressure']
		if InnerPressure != 0.0:
			del compositeSubModel.loads['Load_InnerPressure']
	elif boolCylindricalPressure:
		if cylindricalOuterPressure != 0.0:
			del compositeSubModel.loads['Load_OuterPressure_CylindricalBending']
		if cylindricalInnerPressure != 0.0:
			del compositeSubModel.loads['Load_InnerPressure_CylindricalBending']

def subModelAxialRemoveLoadsHygrothermal():
	del compositeSubModel.predefinedFields['Temp_Initial']
	del compositeSubModel.predefinedFields['Temperature_Difference']

def subModelAxialRemoveKinCoup(noKinCoup):
	# Remove the kinematic coupling:
	for jj in range(noKinCoup):
		compositeSubModel.boundaryConditions.delete(('Fix_RefPoint_Pos_'+  str(jj*int(180*angleOpening/np.pi)),))
		del compositeSubModel.constraints['Kinematic_Coupling_'+  str(jj*int(180*angleOpening/np.pi))]
		compositeSubModelAssembly.deleteSets(setNames=('Set_Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)),))
		del compositeSubModelAssembly.features['RP-'+str(jj+1)]

def subModelAxialRemoveSetsEdgeLoad():
	compositeSubModelAssembly.deleteSets(('Set_Fix_U_Interface_'+ str(0),))
	for ii in range(N):
		compositeSubModelAssembly.deleteSets(('Set_Fix_V_Layer_'+str(ii+1),))
		compositeSubModelAssembly.deleteSets(('Set_Fix_W_Layer_'+str(ii+1),))

def subModelAxialRemoveSetsSurfaceLoad():
	for jj in range(2):
		compositeSubModelAssembly.deleteSets(('Set_Fix_U_SS_Interface_'+ str(jj),))

if boolBendingMoment and not any([boolRadialForce,boolCircForce]):
	subModelAxialRemoveLoadsBendingMoment(noKinCoup)
	subModelAxialRemoveKinCoup(noKinCoup)
elif any([boolBendingMoment,boolRadialForce,boolCircForce]):
	subModelAxialRemoveBCEdgeLoad()
	subModelAxialRemoveLoadsEdgeLoad(noKinCoup)
	subModelAxialRemoveKinCoup(noKinCoup)
	subModelAxialRemoveSetsEdgeLoad()
elif any([boolPressure,boolCylindricalPressure]):
	subModelAxialRemoveKinCoup(noKinCoup)
	subModelAxialRemoveBCSurfaceLoad()
	subModelAxialRemoveLoadsSurfaceLoad()
	subModelAxialRemoveSetsSurfaceLoad()
elif any([boolTempDif,boolMoistDif]):
	subModelAxialRemoveLoadsHygrothermal()
else:
	pass

#---------------------------------------------------------------
def meshCurvedCompositeAxialSubModel():
	# Remesh submodel:
	# Vernetzung in radialer Richtung:
	mRSubModel = 6
	mRInterfaceRatioSubModel = 2.0
	mRInterfaceRatioSubModelIsotropic = 10.0
	# mRConstSubModel = int(round(3*mRSubModel/5))
	mRConstSubModel = mRSubModel
	mRConstRatioSubModel = 1
	for ii in range (len(rFacesAll)):
		for kk in range(len(axialEdgesSubModel)):
			for jj in range(len(thetaEdgesSubModel)):
				compositeMeshEdgeVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rFacesAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rFacesAll[ii], thetaEdgesSubModel[jj]), axialEdges[kk]))
				curvedCompositeMeshEdge = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rFacesAll[ii], thetaEdgesSubModel[jj]), axialEdges[kk] ),))
				
				rI = ct.cart2pol_radius(compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[0]].pointOn[0][1])
				rII = ct.cart2pol_radius(compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[1]].pointOn[0][1])
				
				if plyAngle[int(ii/2)] == 2024:
					if rI > rII:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
					else:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
				else:
					if rI > rII:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
					else:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
	# Vernetzung in tangentialer Richtung:
	mTSubModel = int(2.5*mRSubModel)
	mTRatioSubModel = 15
	for ii in range (len(rInterfaceAll)):
		for kk in range(len(axialEdges)):
			compositeMeshEdgeIVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[0]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[0]), axialEdges[kk] ))
			compositeMeshEdgeI = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[0]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[0]), axialEdges[kk] ),))
			
			thetaI = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][1])
			thetaII = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][1])
			if thetaI < thetaII:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
			else:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
			
			compositeMeshEdgeIIVertice = curvedCompositeSubModelPart.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[-1]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[-1]), axialEdges[kk] ))
			compositeMeshEdgeII = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[-1]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[-1]), axialEdges[kk] ),))
			
			thetaI = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[0]].pointOn[0][1])
			thetaII = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[1]].pointOn[0][1])
			if thetaI < thetaII:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeII, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
			else:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeII, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
	# Vernetzung in Breitenrichtung:
	mZSubModel = int(3*mRSubModel)
	mZRatioSubModel = 10
	for ii in range(len(rInterfaceAll)):
		for jj in range(len(thetaEdgesSubModel)):
			compositeMeshEdgeIVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[0]))
			compositeMeshEdgeI = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[0]),))
			
			if compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][2] < compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][2]:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
			else:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
			
			compositeMeshEdgeIIVertice = curvedCompositeSubModelPart.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[-1] ))
			compositeMeshEdgeII = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[-1] ),))
			
			if compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[0]].pointOn[0][2] > compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[1]].pointOn[0][2]:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeII, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
			else:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeII, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
			
			for kk in range(1,len(axialFacesSubModel)-1):
				curvedCompositeMeshEdge = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[kk] ),))
				if int(round(mZSubModel*(axialEdgesSubModel[kk+1]-axialEdgesSubModel[kk])/length)) < 1:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mZRatio, number=1, constraint=FINER)
				else:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=1, number=int(round(mZSubModel*(axialEdgesSubModel[kk+1]-axialEdgesSubModel[kk])/length)), constraint=FINER)

def meshCurvedCompositeCircSubModel():
	# Remesh submodel:
	# Vernetzung in radialer Richtung:
	mRSubModel = 6
	mRInterfaceRatioSubModel = 2.0
	mRInterfaceRatioSubModelIsotropic = 10.0
	# mRConstSubModel = int(round(3*mRSubModel/5))
	mRConstSubModel = mRSubModel
	mRConstRatioSubModel = 1
	
	for ii in range (len(rFacesAll)):
		for kk in range(len(axialEdgesSubModel)):
			for jj in range(len(thetaEdgesSubModel)):
				compositeMeshEdgeVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rFacesAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rFacesAll[ii], thetaEdgesSubModel[jj]), axialEdgesSubModel[kk]))
				curvedCompositeMeshEdge = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rFacesAll[ii], thetaEdgesSubModel[jj]), axialEdgesSubModel[kk] ),))
				
				rI = ct.cart2pol_radius(compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[0]].pointOn[0][1])
				rII = ct.cart2pol_radius(compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[1]].pointOn[0][1])
				
				if plyAngle[int(ii/2)] == 2024:
					if rI > rII:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
					else:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
				else:
					if rI > rII:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
					else:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
	
	# Vernetzung in tangentialer Richtung:
	mTSubModel = int(5*mRSubModel)
	mTRatioSubModel = 15
	for ii in range (len(rInterfaceAll)):
		for kk in range(len(axialEdgesSubModel)):
			compositeMeshEdgeIVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[0]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[0]), axialEdgesSubModel[kk] ))
			compositeMeshEdgeI = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[0]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[0]), axialEdgesSubModel[kk] ),))
			
			thetaI = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][1])
			thetaII = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][1])
			if thetaI > thetaII:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
			else:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
			
			compositeMeshEdgeIIVertice = curvedCompositeSubModelPart.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[-1]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[-1]), axialEdgesSubModel[kk] ))
			compositeMeshEdgeII = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[-1]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[-1]), axialEdgesSubModel[kk] ),))
			
			thetaI = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[0]].pointOn[0][1])
			thetaII = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[1]].pointOn[0][1])
			if thetaI > thetaII:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeII, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
			else:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeII, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
			
			for jj in range(1,len(thetaFacesSubModel)-1):
				curvedCompositeMeshEdge = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaFacesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaFacesSubModel[jj]), axialEdgesSubModel[kk] ),))
				if int(round(mTSubModel*(thetaEdgesSubModel[jj+1]-thetaEdgesSubModel[jj])/angleOpening)) < 1:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mTRatio, number=1, constraint=FINER)
				else:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=1, number=int(round(mTSubModel*(thetaEdgesSubModel[jj+1]-thetaEdgesSubModel[jj])/angleOpening)), constraint=FINER)
	
	# Vernetzung in Breitenrichtung:
	mZSubModel = int(3*mRSubModel)
	mZRatioSubModel = 10
	for ii in range(len(rInterfaceAll)):
		for jj in range(len(thetaEdgesSubModel)):
			compositeMeshEdgeIVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[0]))
			compositeMeshEdgeI = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[0]),))
			
			if compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][2] > compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][2]:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
			else:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
			
			compositeMeshEdgeIIVertice = curvedCompositeSubModelPart.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[-1] ))
			compositeMeshEdgeII = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaEdgesSubModel[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaEdgesSubModel[jj]), axialFacesSubModel[-1] ),))
			
			if compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[0]].pointOn[0][2] < compositeSubModelInstance.vertices[compositeMeshEdgeIIVertice.getVertices()[1]].pointOn[0][2]:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeII, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
			else:
				curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeII, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)

def meshCurvedCompositeCornerSubModel(quadrant=1):
	cornerAxialFace,cornerThetaFace,cornerAxialEdges,cornerThetaEdges = subModelCheckCorner(quadrant)
	# Remesh submodel:
	# Vernetzung in radialer Richtung:
	mRSubModel = 6
	mRInterfaceRatioSubModel = 2.0
	mRInterfaceRatioSubModelIsotropic = 10.0
	# mRConstSubModel = int(round(3*mRSubModel/5))
	mRConstSubModel = mRSubModel
	mRConstRatioSubModel = 1
	for ii in range (len(rFacesAll)):
		for kk in range(len(cornerAxialEdges)):
			for jj in range(len(cornerThetaEdges)):
				compositeMeshEdgeVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rFacesAll[ii], cornerThetaEdges[jj]), ct.pol2cart_y(rFacesAll[ii], cornerThetaEdges[jj]), cornerAxialEdges[kk]))
				curvedCompositeMeshEdge = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], cornerThetaEdges[jj]), ct.pol2cart_y(rFacesAll[ii], cornerThetaEdges[jj]), cornerAxialEdges[kk] ),))
				
				rI = ct.cart2pol_radius(compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[0]].pointOn[0][1])
				rII = ct.cart2pol_radius(compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[1]].pointOn[0][1])
				
				if plyAngle[int(ii/2)] == 2024:
					if rI > rII:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
					else:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModelIsotropic, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
				else:
					if rI > rII:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
					else:
						if rInterfacesBias[ii] == 1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						elif rInterfacesBias[ii] == -1:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
						else:
							curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
	
	# Vernetzung in tangentialer Richtung:
	mTSubModel = int(2.5*mRSubModel)
	mTRatioSubModel = 15
	for ii in range (len(rInterfaceAll)):
		for kk in range(len(cornerAxialEdges)):
			compositeMeshEdgeIVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], cornerThetaFace), ct.pol2cart_y(rInterfaceAll[ii], cornerThetaFace), cornerAxialEdges[kk] ))
			compositeMeshEdgeI = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], cornerThetaFace), ct.pol2cart_y(rInterfaceAll[ii], cornerThetaFace), cornerAxialEdges[kk] ),))
			
			thetaI = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][1])
			thetaII = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][1])
			if quadrant in [1,4]:
				if thetaI > thetaII:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
				else:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
			else:
				if thetaI < thetaII:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
				else:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=mTSubModel, constraint=FINER)
	
	# Vernetzung in Breitenrichtung:
	mZSubModel = int(3*mRSubModel)
	mZRatioSubModel = 10
	for ii in range(len(rInterfaceAll)):
		for jj in range(len(cornerThetaEdges)):
			compositeMeshEdgeIVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rInterfaceAll[ii], cornerThetaEdges[jj]), ct.pol2cart_y(rInterfaceAll[ii], cornerThetaEdges[jj]), cornerAxialFace))
			compositeMeshEdgeI = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], cornerThetaEdges[jj]), ct.pol2cart_y(rInterfaceAll[ii], cornerThetaEdges[jj]), cornerAxialFace),))
			if quadrant in [1,2]:
				if compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][2] > compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][2]:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
				else:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
			else:
				if compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][2] < compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][2]:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)
				else:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=mZSubModel, constraint=FINER)

def meshCurvedCompositeSurfaceSubModel():
	# Remesh submodel:
	# Vernetzung in radialer Richtung:
	mRSubModel = 6
	mRInterfaceRatioSubModel = 2.0
	mRInterfaceRatioSubModelIsotropic = 10.0
	# mRConstSubModel = int(round(3*mRSubModel/5))
	mRConstSubModel = mRSubModel
	mRConstRatioSubModel = 1
	for ii in range(len(rkSurfaceSubModelFaces)):
		for kk in range(len(axialEdges)):
			for jj in range(len(thetaEdges)):
				compositeMeshEdgeVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rkSurfaceSubModelFaces[ii], thetaEdges[jj]), ct.pol2cart_y(rkSurfaceSubModelFaces[ii], thetaEdges[jj]), axialEdges[kk]))
				curvedCompositeMeshEdge = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rkSurfaceSubModelFaces[ii], thetaEdges[jj]), ct.pol2cart_y(rkSurfaceSubModelFaces[ii], thetaEdges[jj]), axialEdges[kk] ),))
				
				rI = ct.cart2pol_radius(compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[0]].pointOn[0][1])
				rII = ct.cart2pol_radius(compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeVertice.getVertices()[1]].pointOn[0][1])
				
				if rI > rII:
					if rInterfacesSubModelSurfaceBias[ii] == 1:
						curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
					elif rInterfacesSubModelSurfaceBias[ii] == -1:
						curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
					else:
						curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
				else:
					if rInterfacesSubModelSurfaceBias[ii] == 1:
						curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
					elif rInterfacesSubModelSurfaceBias[ii] == -1:
						curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=curvedCompositeMeshEdge, ratio=mRInterfaceRatioSubModel, number=mRSubModel, constraint=FINER)
					else:
						curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRConstRatioSubModel, number=mRConstSubModel, constraint=FINER)
	
	# Vernetzung in tangentialer Richtung:
	mTSubModel = int(10*mRSubModel)
	mTRatioSubModel = 1
	for ii in range (len(rkSurfaceSubModelEdges)):
		for kk in range(len(axialEdges)):
			for jj in range(len(thetaFaces)):
				compositeMeshEdgeIVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rkSurfaceSubModelEdges[ii], thetaFaces[jj]), ct.pol2cart_y(rkSurfaceSubModelEdges[ii], thetaFaces[jj]), axialEdges[kk] ))
				compositeMeshEdgeI = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rkSurfaceSubModelEdges[ii], thetaFaces[jj]), ct.pol2cart_y(rkSurfaceSubModelEdges[ii], thetaFaces[jj]), axialEdges[kk] ),))
				
				thetaI = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][1])
				thetaII = ct.cart2pol_theta(compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][0], compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][1])
				if thetaI > thetaII:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=int(round(mTSubModel*(thetaEdges[jj+1]-thetaEdges[jj])/angleOpening)), constraint=FINER)
				else:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mTRatioSubModel, number=int(round(mTSubModel*(thetaEdges[jj+1]-thetaEdges[jj])/angleOpening)), constraint=FINER)
	# Vernetzung in Breitenrichtung:
	mZSubModel = int(10*mRSubModel)
	mZRatioSubModel = 1
	for ii in range(len(rkSurfaceSubModelEdges)):
		for jj in range(len(thetaEdges)):
			for kk in range(len(axialFaces)):
				compositeMeshEdgeIVertice = compositeSubModelInstance.edges.findAt((ct.pol2cart_x(rkSurfaceSubModelEdges[ii], thetaEdges[jj]), ct.pol2cart_y(rkSurfaceSubModelEdges[ii], thetaEdges[jj]), axialFaces[kk]))
				compositeMeshEdgeI = curvedCompositeSubModelPart.edges.findAt(((ct.pol2cart_x(rkSurfaceSubModelEdges[ii], thetaEdges[jj]), ct.pol2cart_y(rkSurfaceSubModelEdges[ii], thetaEdges[jj]), axialFaces[kk]),))
				
				if compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[0]].pointOn[0][2] > compositeSubModelInstance.vertices[compositeMeshEdgeIVertice.getVertices()[1]].pointOn[0][2]:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=int(round(mZSubModel*(axialEdges[kk+1]-axialEdges[kk])/length)), constraint=FINER)
				else:
					curvedCompositeSubModelPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=compositeMeshEdgeI, ratio=mZRatioSubModel, number=int(round(mZSubModel*(axialEdges[kk+1]-axialEdges[kk])/length)), constraint=FINER)

def meshCurvedComposite():
	# Element types: C3D8, C3D8R, C3D20, C3D20R
	compositeElementType = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
	curvedCompositeSubModelPart.setMeshControls(regions=curvedCompositeSubModelPart.cells, elemShape=HEX, technique=STRUCTURED)
	curvedCompositeSubModelPart.setElementType(regions=(curvedCompositeSubModelPart.cells[:], ), elemTypes=(compositeElementType, ))
	curvedCompositeSubModelPart.seedPart(size=0.05, deviationFactor=0.1, constraint=FINER)
	curvedCompositeSubModelPart.generateMesh()

if boolAxialSubModel:
	meshCurvedCompositeAxialSubModel()
elif boolCircSubModel and not boolCornerSubModel and not boolSurfaceSubModel:
	meshCurvedCompositeCircSubModel()
elif boolCornerSubModel and not boolSurfaceSubModel:
	meshCurvedCompositeCornerSubModel(cornerSubModelQuadrant)
elif boolSurfaceSubModel:
	meshCurvedCompositeSurfaceSubModel()
else:
	pass

meshCurvedComposite()

#---------------------------------------------------------------
# BC - Submodel:
def boundaryConditionCurvedSubModelAxialCoupling():
	curvedCompositeStringSubModelBC =''
	for ii in range (N):
		for kk in range(len(axialFaces)):
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii]+rm[ii])/2, angleSubModelStart), ct.pol2cart_y((rk[ii]+rm[ii])/2, angleSubModelStart), axialFaces[kk]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, angleSubModelStart), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, angleSubModelStart), axialFaces[kk]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii]+rm[ii])/2, angleSubModelStart+angleSubModel), ct.pol2cart_y((rk[ii]+rm[ii])/2, angleSubModelStart+angleSubModel), axialFaces[kk]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, angleSubModelStart+angleSubModel), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, angleSubModelStart+angleSubModel), axialFaces[kk]),),) + ','
	curvedCompositeFacesSubModelBCExec = 'curvedCompositeFacesSubModelBC = compositeSubModelInstance.faces.findAt(' + curvedCompositeStringSubModelBC + ')'
	exec(curvedCompositeFacesSubModelBCExec)
	curvedCompositeSubModelBCRegion = regionToolset.Region(faces=curvedCompositeFacesSubModelBC)
	compositeSubModel.SubmodelBC(name='SubModelBC', createStepName=step, region=curvedCompositeSubModelBCRegion, globalStep='1', globalIncrement=0, timeScale=OFF, dof=(1, 2, 3), globalDrivingRegion='', absoluteExteriorTolerance=None, exteriorTolerance=0.05)

def boundaryConditionCurvedSubModelCircCoupling():
	curvedCompositeStringSubModelBC =''
	for ii in range(N):
		for jj in range(len(thetaFaces)):
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2, thetaFaces[jj]), axialEdgesSubModel[0]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), axialEdgesSubModel[0]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii]+rm[ii])/2, thetaFaces[jj]), axialEdgesSubModel[-1]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, thetaFaces[jj]), axialEdgesSubModel[-1]),),) + ','
	curvedCompositeFacesSubModelBCExec = 'curvedCompositeFacesSubModelBC = compositeSubModelInstance.faces.findAt(' + curvedCompositeStringSubModelBC + ')'
	exec(curvedCompositeFacesSubModelBCExec)
	curvedCompositeSubModelBCRegion = regionToolset.Region(faces=curvedCompositeFacesSubModelBC)
	compositeSubModel.SubmodelBC(name='SubModelBC', createStepName=step, region=curvedCompositeSubModelBCRegion, globalStep='1', globalIncrement=0, timeScale=OFF, dof=(1, 2, 3), globalDrivingRegion='', absoluteExteriorTolerance=None, exteriorTolerance=0.05)

def boundaryConditionCurvedSubModelSurfaceCoupling():
	curvedCompositeStringSubModelBC =''
	for jj in range(len(thetaFaces)):
		for kk in range(len(axialFaces)):
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x(rkSurfaceSubModelEdges[0], thetaFaces[jj]), ct.pol2cart_y(rkSurfaceSubModelEdges[0], thetaFaces[jj]), axialFaces[kk]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x(rkSurfaceSubModelEdges[-1], thetaFaces[jj]), ct.pol2cart_y(rkSurfaceSubModelEdges[-1], thetaFaces[jj]), axialFaces[kk]),),) + ','
	curvedCompositeFacesSubModelBCExec = 'curvedCompositeFacesSubModelBC = compositeSubModelInstance.faces.findAt(' + curvedCompositeStringSubModelBC + ')'
	exec(curvedCompositeFacesSubModelBCExec)
	curvedCompositeSubModelBCRegion = regionToolset.Region(faces=curvedCompositeFacesSubModelBC)
	compositeSubModel.SubmodelBC(name='SubModelBC', createStepName=step, region=curvedCompositeSubModelBCRegion, globalStep='1', globalIncrement=0, timeScale=OFF, dof=(1, 2, 3), globalDrivingRegion='', absoluteExteriorTolerance=None, exteriorTolerance=0.05)

def boundaryConditionCurvedSubModelCornerCoupling(quadrant=1):
	cornerAxialFace,cornerThetaFace,cornerAxialEdges,cornerThetaEdges = subModelCheckCorner(quadrant)
	curvedCompositeStringSubModelBC =''
	for ii in range(N):
		if quadrant == 1 or quadrant == 2:
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii]+rm[ii])/2, cornerThetaFace), ct.pol2cart_y((rk[ii]+rm[ii])/2, cornerThetaFace), cornerAxialEdges[0]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, cornerThetaFace), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, cornerThetaFace), cornerAxialEdges[0]),),) + ','
		elif quadrant == 3 or quadrant == 4:
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii]+rm[ii])/2, cornerThetaFace), ct.pol2cart_y((rk[ii]+rm[ii])/2, cornerThetaFace), cornerAxialEdges[-1]),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, cornerThetaFace), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, cornerThetaFace), cornerAxialEdges[-1]),),) + ','
		if quadrant == 1 or quadrant == 4:
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii]+rm[ii])/2, cornerThetaEdges[-1]), ct.pol2cart_y((rk[ii]+rm[ii])/2, cornerThetaEdges[-1]), cornerAxialFace),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, cornerThetaEdges[-1]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, cornerThetaEdges[-1]), cornerAxialFace),),) + ','
		elif quadrant == 2 or quadrant == 3:
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii]+rm[ii])/2, cornerThetaEdges[0]), ct.pol2cart_y((rk[ii]+rm[ii])/2, cornerThetaEdges[0]), cornerAxialFace),),) + ','
			curvedCompositeStringSubModelBC = curvedCompositeStringSubModelBC + str(((ct.pol2cart_x((rk[ii+1]+rm[ii])/2, cornerThetaEdges[0]), ct.pol2cart_y((rk[ii+1]+rm[ii])/2, cornerThetaEdges[0]), cornerAxialFace),),) + ','
	curvedCompositeFacesSubModelBCExec = 'curvedCompositeFacesSubModelBC = compositeSubModelInstance.faces.findAt(' + curvedCompositeStringSubModelBC + ')'
	exec(curvedCompositeFacesSubModelBCExec)
	curvedCompositeSubModelBCRegion = regionToolset.Region(faces=curvedCompositeFacesSubModelBC)
	compositeSubModel.SubmodelBC(name='SubModelBC', createStepName=step, region=curvedCompositeSubModelBCRegion, globalStep='1', globalIncrement=0, timeScale=OFF, dof=(1, 2, 3), globalDrivingRegion='', absoluteExteriorTolerance=None, exteriorTolerance=0.05)

def loadsSubModelSurfaceLoad():
	if OuterPressure != 0.0 or cylindricalOuterPressure != 0.0:
		# Aeussere Last - Aussendruck:
		curvedCompositeStringFaces = ''
		for jj in range(len(thetaPartition)-1):
			for kk in range(len(axialFaces)):
				curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[-1], (thetaPartition[jj]+thetaPartition[jj+1])/2), ct.pol2cart_y(rk[-1], (thetaPartition[jj]+thetaPartition[jj+1])/2), axialFaces[kk]),),) + ','
		
		curvedCompositeStringFacesExec = 'curvedCompositeOuterPressureFaces = compositeSubModelInstance.faces.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		
		if boolCylindricalPressure:
			compositeSubModel.ExpressionField(name='Load_OuterPressure_CylindricalBending', localCsys=curvedCompositeCylCoordSys, description='',expression='sin(pi*Th/angleOpening)')
			compositeSubModel.Pressure(name='Load_OuterPressure_CylindricalBending', createStepName=step, region=regionToolset.Region(side1Faces=curvedCompositeOuterPressureFaces), magnitude=cylindricalOuterPressure, distributionType=FIELD, field='Load_OuterPressure_CylindricalBending', amplitude=UNSET)
		else:
			compositeSubModel.Pressure(name='Load_OuterPressure', createStepName=step, region=regionToolset.Region(side1Faces=curvedCompositeOuterPressureFaces), magnitude=OuterPressure, distributionType=UNIFORM)
	elif InnerPressure != 0.0 or cylindricalInnerPressure != 0.0:
		# Aeussere Last - Innendruck:
		curvedCompositeStringFaces = ''
		for jj in range(len(thetaPartition)-1):
			for kk in range(len(axialFaces)):
				curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[0], (thetaPartition[jj]+thetaPartition[jj+1])/2), ct.pol2cart_y(rk[0], (thetaPartition[jj]+thetaPartition[jj+1])/2), axialFaces[kk]),),) + ','
		
		curvedCompositeStringFacesExec = 'curvedCompositeInnerPressureFaces = compositeSubModelInstance.faces.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		
		if boolCylindricalPressure:
			compositeSubModel.ExpressionField(name='Load_InnerPressure_CylindricalBending', localCsys=curvedCompositeCylCoordSys, description='',expression='sin(pi*Th/angleOpening)')
			compositeSubModel.Pressure(name='Load_InnerPressure_CylindricalBending', createStepName=step, region=regionToolset.Region(side1Faces=curvedCompositeInnerPressureFaces), magnitude=cylindricalInnerPressure, distributionType=FIELD, field='Load_InnerPressure_CylindricalBending', amplitude=UNSET)
		else:
			compositeSubModel.Pressure(name='Load_InnerPressure', createStepName=step, region=regionToolset.Region(side1Faces=curvedCompositeInnerPressureFaces), magnitude=InnerPressure, distributionType=UNIFORM)

def loadsSubModelHygrothermal():
	# Definition der Ausgangstemperatur:
	compositeSubModel.Temperature(name='Temp_Initial_SubModel', createStepName=stepPrevious, region=regionToolset.Region(cells=compositeSubModelInstance.cells), distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0,))
	compositeSubModel.Temperature(name='Temperature_Difference_SubModel', createStepName=step, region=regionToolset.Region(cells=compositeSubModelInstance.cells), distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(tempDif,))

if boolAxialSubModel:
	boundaryConditionCurvedSubModelAxialCoupling()
elif boolCircSubModel and not boolCornerSubModel and not boolSurfaceSubModel:
	boundaryConditionCurvedSubModelCircCoupling()
elif boolCornerSubModel and not boolSurfaceSubModel:
	boundaryConditionCurvedSubModelCornerCoupling(cornerSubModelQuadrant)
elif boolSurfaceSubModel:
	boundaryConditionCurvedSubModelSurfaceCoupling()
else:
	pass

if any([boolPressure,boolCylindricalPressure]):
	loadsSubModelSurfaceLoad()
elif any([boolTempDif,boolMoistDif]):
	loadsSubModelHygrothermal()
else:
	pass

#---------------------------------------------------------------
def setsCurvedCompositeAxial():
	# Sets zur Verifizierung der Free Edge-Solution:
	# Sets der Spannungen in radialer Richtung an beiden Enden des Laminats:
	for ii in range(N):
		compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(((ct.pol2cart_x((rm[ii]+rk[ii])/2, angleEval), ct.pol2cart_y((rm[ii]+rk[ii])/2, angleEval), axialEdges[0]),),((ct.pol2cart_x((rm[ii]+rk[ii+1])/2, angleEval), ct.pol2cart_y((rm[ii]+rk[ii+1])/2, angleEval), axialEdges[0]),))
		
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FEFreeEdgeAxial_'+str(int(axialEdges[0]))+'_Layer'+str(ii+1))
		
		compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(((ct.pol2cart_x((rm[ii]+rk[ii])/2, angleEval), ct.pol2cart_y((rm[ii]+rk[ii])/2, angleEval), axialEdges[-1]),),((ct.pol2cart_x((rm[ii]+rk[ii+1])/2, angleEval), ct.pol2cart_y((rm[ii]+rk[ii+1])/2, angleEval), axialEdges[-1]),))
		
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FEFreeEdgeAxial_'+str(int(axialEdges[-1]))+'_Layer'+str(ii+1))
	# Sets zur Bestimmung des Abklingverhaltens der singulaeren Spannungen (in Dickenrichtung) in den Auswertungsinterfaces an beiden Enden:
	for kk in range(len(iInterfaceEval)):
		curvedCompositeStringFaces = ''
		ii = iInterfaceEval[kk]
		for jj in range(len(axialFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[ii], angleEval), ct.pol2cart_y(rk[ii], angleEval), axialFaces[jj]),),) + ','
		curvedCompositeStringFacesExec = 'compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FEInterfaceAxial'+str(iInterfaceEval[kk]))
	# Sets zur Bestimmung der Spannungsverlaeufe innerhalb jeder physikalischen Schicht in Dickenrichtung:
	curvedCompositeStringFaces = ''
	for ii in range(N):
		curvedCompositeStringFaces = ''
		for jj in range(len(axialFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rm[ii], angleEval), ct.pol2cart_y(rm[ii], angleEval), axialFaces[jj]),),) + ','
		curvedCompositeStringFacesExec = 'compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FELayerAxial'+str(ii+1))
	# Sets zur Verifizierung der Inner-Solution:
	# Sets der Spannungen in radialer Richtung in Laminatmitte:
	for ii in range(N):
		compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(((ct.pol2cart_x((rm[ii]+rk[ii])/2, angleEval), ct.pol2cart_y((rm[ii]+rk[ii])/2, angleEval), axialEval),),((ct.pol2cart_x((rm[ii]+rk[ii+1])/2, angleEval), ct.pol2cart_y((rm[ii]+rk[ii+1])/2, angleEval), axialEval),))
		
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FEInner'+str(ii+1))

def setsCurvedCompositeCirc():
	# Sets zur Bestimmung des Abklingverhaltens der singulaeren Spannungen (in Dickenrichtung) in den Auswertungsinterfaces an beiden Enden:
	for kk in range(len(iInterfaceEval)):
		curvedCompositeStringFaces = ''
		ii = iInterfaceEval[kk]
		for jj in range(len(thetaFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[ii], thetaFaces[jj]), ct.pol2cart_y(rk[ii], thetaFaces[jj]), axialEval),),) + ','
		curvedCompositeStringFacesExec = 'compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FEInterfaceCirc'+str(iInterfaceEval[kk]))
	# Sets zur Bestimmung der Spannungsverlaeufe innerhalb jeder physikalischen Schicht in Dickenrichtung:
	curvedCompositeStringFaces = ''
	for ii in range(N):
		curvedCompositeStringFaces = ''
		for jj in range(len(thetaFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rm[ii], thetaFaces[jj]), ct.pol2cart_y(rm[ii], thetaFaces[jj]), axialEval),),) + ','
		curvedCompositeStringFacesExec = 'compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FELayerCirc'+str(ii+1))
	# Sets zur Verifizierung der Free Edge-Solution:
	# Sets der Spannungen in radialer Richtung an beiden Enden des Laminats:
	for ii in range(N):
		compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(((ct.pol2cart_x((rm[ii]+rk[ii])/2, thetaEdgesSubModel[0]), ct.pol2cart_y((rm[ii]+rk[ii])/2, thetaEdgesSubModel[0]), axialEval),),((ct.pol2cart_x((rm[ii]+rk[ii+1])/2, thetaEdgesSubModel[0]), ct.pol2cart_y((rm[ii]+rk[ii+1])/2, thetaEdgesSubModel[0]), axialEval),))
		
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FEFreeEdgeCirc_'+str(int(thetaEdgesSubModel[0]*180/np.pi))+'_Layer'+str(ii+1))
		
		compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(((ct.pol2cart_x((rm[ii]+rk[ii])/2, thetaEdgesSubModel[-1]), ct.pol2cart_y((rm[ii]+rk[ii])/2, thetaEdgesSubModel[-1]), axialEval),),((ct.pol2cart_x((rm[ii]+rk[ii+1])/2, thetaEdgesSubModel[-1]), ct.pol2cart_y((rm[ii]+rk[ii+1])/2, thetaEdgesSubModel[-1]), axialEval),))
		
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FEFreeEdgeCirc_'+str(int(thetaEdgesSubModel[-1]*180/np.pi))+'_Layer'+str(ii+1))

def setsCurvedCompositeCorner(quadrant=1):
	_,_,cornerAxialEdges,cornerThetaEdges = subModelCheckCorner(quadrant)
	if quadrant == 1 or quadrant == 4:
		thetaEvalCorner = cornerThetaEdges[0]
	else:
		thetaEvalCorner = cornerThetaEdges[-1]
	if quadrant == 1 or quadrant == 2:
		axialEvalCorner = cornerAxialEdges[-1]
	else:
		axialEvalCorner = cornerAxialEdges[0]
	for ii in range(N):
		compositeSetEvalEdge = compositeSubModelInstance.edges.findAt(((ct.pol2cart_x((rm[ii]+rk[ii])/2, thetaEvalCorner), ct.pol2cart_y((rm[ii]+rk[ii])/2, thetaEvalCorner), axialEvalCorner),),((ct.pol2cart_x((rm[ii]+rk[ii+1])/2, thetaEvalCorner), ct.pol2cart_y((rm[ii]+rk[ii+1])/2, thetaEvalCorner), axialEvalCorner),))
		
		compositeSubModelAssembly.Set(edges=compositeSetEvalEdge, name='FEFreeEdgeCorner_t'+str(int(thetaEvalCorner*180/np.pi))+'_z'+str(int(axialEvalCorner))+'_Layer'+str(ii+1))

def setsCurvedCompositeSurface():
	# Sets zur Bestimmung des Abklingverhaltens der singulaeren Spannungen (in Dickenrichtung) in den Auswertungsinterfaces an beiden Enden:
	curvedCompositeStringFaces = ''
	for jj in range(len(thetaFaces)):
		for kk in range(len(axialFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rkSurfaceSubModelEdges[1], thetaFaces[jj]), ct.pol2cart_y(rkSurfaceSubModelEdges[1], thetaFaces[jj]), axialFaces[kk]),),) + ','
	curvedCompositeStringFacesExec = 'compositeSetEvalFace = compositeSubModelInstance.faces.findAt(' + curvedCompositeStringFaces + ')'
	exec(curvedCompositeStringFacesExec)
	compositeSubModelAssembly.Set(faces=compositeSetEvalFace, name='FEInterfaceSurface'+str(surfaceSubModelInterface))

if boolAxialSubModel:
	setsCurvedCompositeAxial()
elif boolCircSubModel and not boolCornerSubModel and not boolSurfaceSubModel:
	setsCurvedCompositeCirc()
elif boolCornerSubModel and not boolSurfaceSubModel:
	setsCurvedCompositeCorner(cornerSubModelQuadrant)
elif boolSurfaceSubModel:
	setsCurvedCompositeSurface()
else:
	pass

# ---------------------------------------------------------------
# Field Output:
compositeSubModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'COORD'))

#-------------------------------------------------------------------
# Job
# Erstellen des Jobs:
jobNameSubModel = modelName + '_Job_SubModel'
mdb.Job(name=jobNameSubModel, model=modelName+'_SubModel',description='Run FE-analysis', parallelizationMethodExplicit=DOMAIN, numDomains=12, numCpus=12, memory=85,echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)

# Speichern der CAE-Datei:
mdb.saveAs(pathName=analysis_newPath + '\\' + modelName + '_CAE_Model')

# Job starten:
mdb.jobs[jobNameSubModel].submit(consistencyChecking=OFF)
mdb.jobs[jobNameSubModel].waitForCompletion()

#-------------------------------------------------------------------
# Visualization:
compositeSubModelOdbPath = jobNameSubModel + '.odb'
compositeSubModelOdbObject = session.openOdb(name=compositeSubModelOdbPath)
session.viewports['Viewport: 1'].setValues(displayedObject=compositeSubModelOdbObject)
compositeSubModelViewport = session.viewports['Viewport: 1']

# Erstellen eines zylindrischen Koordinatensystems fuer das Postprocessing:
postProc_Cyl_CS = compositeSubModelOdbObject.rootAssembly.DatumCsysByThreePoints(coordSysType = CYLINDRICAL, name='Fixed cylindrical coordinate system - PostProcessing',
																				point1 = (rk[-1], 0, length/2),
																				point2 = (0, rk[-1], length/2),
																				origin= (0, 0, length/2))
compositeSubModelViewport.odbDisplay.basicOptions.setValues(transformationType=USER_SPECIFIED, datumCsys=postProc_Cyl_CS)

# Anzeigen der Radialspannungen am verformten FE-Model:
compositeSubModelViewport.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S11'), )

startTime = time.time()
#-------------------------------------------------------------------
# Postprocessing:
def postProcessingCurvedCompositeAxial():
	# Auswertung der im Preprocessing definierten Knoten-/Element-Sets mithilfe der definierten Funktionen:
	thetaInterfaceEvalDegree = 180*angleEval/np.pi
	for ii in range(len(iInterfaceEval)):
		stressesExtrapolatedAndIntegrationPointInterfaceAxial('FEInterfaceAxial'+str(iInterfaceEval[ii]),step,iInterfaceEval[ii],rInterfaceEval[ii],thetaInterfaceEvalDegree)
	
	for ii in range(N):
		stressesExtrapolatedAndIntegrationPointLayerAxial('FELayerAxial'+str(ii+1),step,rm[ii],thetaInterfaceEvalDegree)
		stressesExtrapolatedAndIntegrationPointRadial('FEInner'+str(ii+1),step,thetaInterfaceEvalDegree,axialEval)
		stressesExtrapolatedAndIntegrationPointRadial('FEFreeEdgeAxial_'+str(int(axialEdges[0]))+'_Layer'+str(ii+1),step, thetaInterfaceEvalDegree, axialEdges[0])
		stressesExtrapolatedAndIntegrationPointRadial('FEFreeEdgeAxial_'+str(int(axialEdges[-1]))+'_Layer'+str(ii+1),step, thetaInterfaceEvalDegree, axialEdges[-1])

def postProcessingCurvedCompositeCirc():
	# Auswertung der im Preprocessing definierten Knoten-/Element-Sets mithilfe der definierten Funktionen:
	for ii in range(len(iInterfaceEval)):
		stressesExtrapolatedAndIntegrationPointInterfaceCirc('FEInterfaceCirc'+str(iInterfaceEval[ii]),step,iInterfaceEval[ii],rInterfaceEval[ii],axialEval)
	
	for ii in range(N):
		stressesExtrapolatedAndIntegrationPointLayerCirc('FELayerCirc'+str(ii+1),step,rm[ii],axialEval)
		stressesExtrapolatedAndIntegrationPointRadial('FEFreeEdgeCirc_'+str(int(thetaEdgesSubModel[0]*180/np.pi))+'_Layer'+str(ii+1),step,thetaEdgesSubModel[0]*180/np.pi, axialEval)
		stressesExtrapolatedAndIntegrationPointRadial('FEFreeEdgeCirc_'+str(int(thetaEdgesSubModel[-1]*180/np.pi))+'_Layer'+str(ii+1),step,thetaEdgesSubModel[-1]*180/np.pi, axialEval)

def postProcessingCurvedCompositeCorner(quadrant=1):
	_,_,cornerAxialEdges,cornerThetaEdges = subModelCheckCorner(quadrant)
	if quadrant == 1 or quadrant == 4:
		thetaEvalCorner = cornerThetaEdges[0]
	else:
		thetaEvalCorner = cornerThetaEdges[-1]
	if quadrant == 1 or quadrant == 2:
		axialEvalCorner = cornerAxialEdges[-1]
	else:
		axialEvalCorner = cornerAxialEdges[0]
	for ii in range(N):
		stressesExtrapolatedAndIntegrationPointRadial('FEFreeEdgeCorner_t'+str(int(thetaEvalCorner*180/np.pi))+'_z'+str(int(axialEvalCorner))+'_Layer'+str(ii+1),step,thetaEvalCorner*180/np.pi, axialEvalCorner)

def postProcessingCurvedCompositeSurface():
		stressesExtrapolatedAndIntegrationPointSurface('FEInterfaceSurface'+str(surfaceSubModelInterface),step)

if boolAxialSubModel:
	postProcessingCurvedCompositeAxial()
elif boolCircSubModel and not boolCornerSubModel and not boolSurfaceSubModel:
	postProcessingCurvedCompositeCirc()
elif boolCornerSubModel and not boolSurfaceSubModel:
	postProcessingCurvedCompositeCorner(cornerSubModelQuadrant)
elif boolSurfaceSubModel:
	postProcessingCurvedCompositeSurface()
else:
	pass

print(time.time()-startTime)
