# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:01:48 2024

@author: psingh
"""

# RigidEdges.unsetPrimaryObject()
import __main__
import os
import time
import datetime as dt
import operator

from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgot
import displayGroupMdbToolset as dgmt
import connectorBehavior
import numpy as np
import math

path_modules = 'N:\\Sachdeva\\MT_Nair\\ABAQUS\\MT\\Macros'

os.chdir(path_modules)

# Further packages:
import coordinateTransformation as ct

analysis_Path = 'C:\\Users\\anair\\Documents\\01_Promotion\\01_ANA_2.1\\FE\\2D_bending_monoclinic'

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

############################################################################
# PARAMETER ZUM ANPASSEN
# Verschiebung der Bolten
disp_u1 =0.00000001 # 4  
#Auswertungwinkel
angleEval = np.pi / 4
#Einstellung, ob bei 0° ausgewertet werden soll
Null_Auswaertung = True

# Meshparameter für den gekrümmten Bereich
# Anzahl Elemente radial, tangential
mR, mT = 4,150
mFlanken = mT/3
# Ratio für gekrümmten Bereich
mRRatio, mTRatio = 1, 1

#Meshparameter für Flanken
#mFlanken = 500  # 50
#Verteilungsparameter
mBias = mTRatio*2   
rBias = 7  

#Meshparameter für Bolzen
mRigid = 0.35
meshsize = 000  # nur für das gewebte Modell wichtig, wird aber für Export gebraucht

#Pruefstandspositionierung
lt = 75
lb = 100
#Durchmesser der Bolzen
drigid = 10

#Probekörpergeometrie
#Innenradius
ri = 6.4
#Länge der Tangente
ltan=10.53
#3D tiefe
length=2

compositeMaterialName = 'cfk'  # 'cfk', 'cfkDuro', 'cfk'
AngleName = 'Test'  # 'bestSymm', 'bestUnsymm', 'Vergleichslaminat', 'Test'

if compositeMaterialName == 'cfk':
    if AngleName == 'Test':
        ric_plyAngle = '0110' 
    elif AngleName == 'Thick':
        ric_plyAngle = '0000111100001111'
    elif AngleName == 'bestSymm':
        ric_plyAngle = '00010110111010111011011101011101101000'  
    elif AngleName == 'onlyB':
        ric_plyAngle = '00011011011011011011011011011011011000'  
    elif AngleName == 'bestUnsymm':
        ric_plyAngle = '00111011101110100101101110111011001000' 
    elif AngleName == 'Vergleichslaminat':
        ric_plyAngle = '00110010011001011100101001101101001110'  
    elif AngleName == 'Foto':
        ric_plyAngle = '0010110'  
     
elif compositeMaterialName == 'gfk':
    if AngleName == 'Test':
        ric_plyAngle = '01010101010101010101010101'  
    elif AngleName == 'bestSymm':
        ric_plyAngle = '00101101101011010110110100'  
    elif AngleName == 'bestUnsymm':
        ric_plyAngle = '00110110110110110110100100' 
    elif AngleName == 'Vergleichslaminat':
        ric_plyAngle = '00110010010011001001001100'  


############################################################################
# Bestimmung des Zeitpunktes der FE-Analyse:
analysis_currentDateTime = dt.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
modelName = 'mR' + str(mR) +'_mT'+ str(mT) + '_' + '_'+ analysis_currentDateTime + '_NORMAL'
#modelName = '0_90_2_R4h_L4h_T1_AO90_AE45_CFK_PS'
# Erstellen eines neuen Ordners in dem angebeben Pfad und Definition dessen als Arbeitsverzeichnis:
analysis_newPath = analysis_Path + '\\' + modelName
os.makedirs(analysis_newPath)
os.chdir(analysis_newPath)
#------------------------------------------------

#------------------------------------------------
# Transformieren in  0 90 Schreibweise
plyAngle = [90 if char == '1' else 0 for char in ric_plyAngle]

#Anzahl und Orientierung der Schichten
N = len(plyAngle)
# Auswertungsinterfaces (Interfaces zweier aufeinanderfolgender Schichten mit unterschiedlichen Faserorientierungswinkeln):
iInterfaceEval = list(range(1,N))
# Startwinkel des Modells, Winkelbereich des Sub-Modells, Oeffnungswinkel,
angleStart,angleOpening = 0.0,np.pi/2

# Materialparameter der 0-Schicht (1 - Radialkoordinate, 2 - Tangentialkoordinate, 3 - Axialkoordinate):
if compositeMaterialName == 'cfk':
	E1, E2, E3 = 7460.0, 118148.0, 7460.0
	Nu12, Nu13, Nu23 = 0.021, 0.37, 0.34
	G12, G13, G23 = 4800.0, 2701.0, 4800.0
	alpha11,alpha22,alpha33=2.6e-5,-1.0e-6,2.6e-5
	dsingle = 0.125
elif compositeMaterialName == 'gfk':
	E1, E2, E3 = 9552.6, 39296.0, 9552.6
	Nu21, Nu13, Nu23 = 0.29, 0.38, 0.29
	Nu12=E1/E2*Nu21
	G12, G13, G23 = 3080.5, 3449.0, 3080.5
	alpha11,alpha22,alpha33=2.6e-5,8.6e-6,2.6e-5
	dsingle = 0.190
elif compositeMaterialName == 'cfkDuro':
	E1, E2, E3 = 6895.0, 172375.0, 6895.0
	Nu21, Nu13, Nu23 = 0.25, 0.25, 0.25
	Nu12=E1/E2*Nu21
	G12,G13,G23=3448.0,1379.0,3448.0
	alpha11,alpha22,alpha33=2.6e-5,-1.0e-6,2.6e-5
	dsingle = 0.125
else:
	pass

#--------------------------------------------------------------
#Berechnung verschiedener Punkte

#Gesamtdicke des Pruefkörpers
d=dsingle*N

#Umrechnung der BolzenCoords in Flankenlänge
lt_transf = lt/np.cos(np.pi/4)/2 - drigid/2
lb_transf = lb/np.cos(np.pi/4)/2 + drigid/2
#Außenradius
ra=ri+d

#Mittelpunkt der Kurvenradien
centerpoint=(ra,ra)
rc = np.sqrt(2*ra**2)
#Äußerer Halbkreis Punkt oben x y
outer_upper_arc = [0,ra]

#Äußerer Halbkreis Punkt unten
outer_lower_arc = [ra,0]

#Innerer Halbkreis punkt oben
inner_upper_arc = [d,ra]

#innerer Halbkreis Punkt unten
inner_lower_arc = [ra,d]

#Radien der Mittelflächen in der Kruemmung
rk =[ri]
for ii in range (N):
    rk.append (rk[ii]+dsingle)

rN = rk[-1]
#--------------------------------------------------------------
#Parameter
def radialGeometryParameters():
    """
    Berechnet und gibt verschiedene radiale Geometrieparameter für Laminatschichten zurück.

    Returns:
        tuple: Ein Tupel, das Folgendes enthält:
            - rm (list): Radien der Mittelflächen der Laminat-Einzelschichten.
            - rEdgesEval (list): Radien aller zu evaluierenden Interfaces.
            - rEdgesAll (list): Radien aller Partitionierungsinterfaces.
            - rFacesAll (list): Mittlere Radien aller partitionierten physikalischen Schichten.
    """
    # Radien der Mittelflächen der Laminat-Einzelschichten:
    rm = [(rk[ii] + rk[ii + 1]) / 2 for ii in range(N)]
     
    # Radien aller zu evaluierenden Interfaces:
    rEdgesEval = [rk[iInterEval] for iInterEval in iInterfaceEval]
     
    # Radien aller Interfaces:
    rInterfaceEval = [rk[iInterEval] for iInterEval in iInterfaceEval]
     
    # Radien aller Partitionierungsinterfaces:
    rEdgesAll = []
    for ii in range(N + 1):
        rEdgesAll.append(rk[ii])
     
    # Mittlere Radien aller partitionierten physikalischen Schichten:
    rFacesAll = [(rEdgesAll[ii] + rEdgesAll[ii + 1]) / 2 for ii in range(N)]
     
    return (rm, rEdgesEval, rEdgesAll, rInterfaceEval, rFacesAll)

def circGeometryParameters():
    """
    Berechnet und gibt verschiedene Umfangswinkel-Geometrieparameter zurück.

    Returns:
        tuple: Ein Tupel, das Folgendes enthält:
            - thetaPartition (list): Umfangswinkel-Partitionierungen.
            - thetaEdges (list): Umfangswinkel-Kanten.
            - thetaFaces (list): Umfangswinkel-Mittenflächen.
            - thetaPartitionFaces (list): Mittlere Umfangswinkel-Partitionierungen.
    """
    # Umfangswinkel-Kanten:
    thetaEdges = [angleStart, angleEval, angleOpening]
     
    # Umfangswinkel-Partitionierungen:
    thetaPartition = [angleEval]
     
    # Umfangswinkel-Mittenflächen:
    refinementInterval = np.pi * 10 / 180
    thetaPartiotionsAll = [
        angleStart,
        angleEval - refinementInterval,
        angleEval,
        angleEval + refinementInterval,
        angleOpening
    ]
    thetaPartitionFaces = [(thetaPartiotionsAll[ii] + thetaPartiotionsAll[ii + 1]) / 2 for ii in range(len(thetaPartiotionsAll) - 1)]
    thetaFaces = [(thetaEdges[ii] + thetaEdges[ii + 1]) / 2 for ii in range(len(thetaEdges) - 1)]
     
    return (thetaPartition, thetaEdges, thetaFaces, thetaPartitionFaces)

rm,rEdgesEval,rEdgesAll,rInterfaceEval,rFacesAll = radialGeometryParameters()
thetaPartition,thetaEdges,thetaFaces,thetaPartitionFaces = circGeometryParameters()

#--------------------------------------------------------------
#Dateien für Auswärtung
Ratio = (ri + d) / ri 

def export_matrix(matrix, filename, delimiter='\t'):
    """
    Exportiert eine Matrix in eine Datei mit einem bestimmten Trennzeichen.

    Args:
        matrix (list of lists): Die zu exportierende Matrix.
        filename (str): Der Name der Datei, in die die Matrix exportiert wird.
        delimiter (str, optional): Das Trennzeichen, das zum Trennen der Werte verwendet wird. Standard ist Tabulator ('\t').
    """
    with open(filename, 'w') as file:
        for row in matrix:
            formatted_row = delimiter.join(
                "{:.0f}".format(x) if x.is_integer() else "{:.2f}".format(x) for x in row
            )
            file.write(formatted_row + '\n')

def materialParametersExport():
    """
    Exportiert die Materialparameter als 'materialParameters'.
    """
    # Erstellen einer 3x3 Matrix für die Konfiguration des Composite-Materials
    matConfigCompositeMatrix = np.array([
        [E1, E2, E3],
        [G12, G13, G23],
        [Nu12, Nu13, Nu23]
    ])
    
    # Definieren des Trennzeichens
    delimiter = '\t'  # Tab-getrennte Werte
    
    # Export der Matrix mit benutzerdefiniertem Format
    with open('materialParameters', 'w') as file:
        for row in matConfigCompositeMatrix:
            formatted_row = delimiter.join("{:.2f}".format(x).lstrip('0') if x != 0 else '0' for x in row)
            file.write(formatted_row + '\n')

materialParametersExport()


def export_geometry(plyAngle, MDetailed, orderLagrangeDetailed, Ratio, dL, width_h, angleBegin, openingAngle, angleEval):
    """
    Exportiert die Geometrie- und Lagenwinkelinformationen in zwei separate Dateien.

    Args:
        plyAngle (list): Liste der Lagenwinkel.
        MDetailed (int): Detaillierungsgrad des Modells.
        orderLagrangeDetailed (int): Ordnung des Lagrange-Polynoms.
        Ratio (float): Verhältnis von (ri + d) zu ri.
        dL (float): Länge.
        width_h (float): Breite.
        angleBegin (float): Startwinkel in Radiant.
        openingAngle (float): Öffnungswinkel in Radiant.
        angleEval (float): Evaluierungswinkel in Radiant.
    """
    # Exportiert plyAngle zu compositeLayup
    N = len(plyAngle)
    pAExport = np.zeros((1, N))
    for ii in range(N):
        pAExport[0, ii] = plyAngle[ii]
    export_matrix(pAExport, 'compositeLayup', delimiter='\t')
    
    # Exportiert geometryCompositeExport zu compositeGeometry
    geometryCompositeExport = np.zeros((1, 8))
    geometryCompositeExport[0, 0] = MDetailed
    geometryCompositeExport[0, 1] = orderLagrangeDetailed
    geometryCompositeExport[0, 2] = Ratio
    geometryCompositeExport[0, 3] = float(dL)
    geometryCompositeExport[0, 4] = width_h
    geometryCompositeExport[0, 5] = 180 * angleBegin / np.pi
    geometryCompositeExport[0, 6] = 180 * openingAngle / np.pi
    geometryCompositeExport[0, 7] = 180 * angleEval / np.pi
    export_matrix(geometryCompositeExport, 'compositeGeometry', delimiter='\t')


# Platzhalter Werte:
MDetailed = 8.0  # Example value
orderLagrangeDetailed = 3  # Example value
Ratio = 4  # Example value
dL = dsingle  # Example value
width_h = 10  # Example value
angleBegin = 0  # Example begin angle in radians
openingAngle = np.pi *0.5  # Example opening angle in radians

export_geometry(plyAngle, MDetailed, orderLagrangeDetailed, Ratio, dL, width_h, angleBegin, openingAngle, angleEval)



#--------------------------------------------------------------
#SKIZZE ERSTELLEN
Mdb()
# Name ändern
mdb.models.changeKey(fromName='Model-1', toName= modelName)
LModel = mdb.models[modelName]

def sketchLModel():
    """
    Erstellt die Skizze des L-Profils und die IDs der äußeren Elemente zurück.

    Returns:
        tuple:
            - LModelPart: Die erstellte Part-Instanz des L-Profils.
            - OuterCirc_id: Die ID des äußeren Halbkreises.
            - LowerOuterTan_id: Die ID der unteren äußeren Tangente.
            - UpperOuterTan_id: Die ID der oberen äußeren Tangente.
    """
    LModelSketch = LModel.ConstrainedSketch(name='L-Profil-Sketch', sheetSize=(2 * (ltan + ra)))
     
    # Innerer Halbkreis
    LModelSketch.ArcByCenterEnds(
        center=(0, 0), 
        point1=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0)), 
        point2=(ct.pol2cart_x(rk[0], thetaEdges[-1]), ct.pol2cart_y(rk[0], thetaEdges[-1])), 
        direction=COUNTERCLOCKWISE)
     
    # Äußererer Halbkreis
    OuterCirc = LModelSketch.ArcByCenterEnds(
        center=(0, 0), 
        point1=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0)), 
        point2=(ct.pol2cart_x(rk[-1], thetaEdges[-1]), ct.pol2cart_y(rk[-1], thetaEdges[-1])), 
        direction=COUNTERCLOCKWISE)
    OuterCirc_id = OuterCirc.id
     
    # Obere Innere Tangente
    LModelSketch.Line(
        point1=(ct.pol2cart_x(rk[0], thetaEdges[-1]), ct.pol2cart_y(rk[0], thetaEdges[-1])), 
        point2=(ct.pol2cart_x(rk[0], thetaEdges[-1]) - ltan, ct.pol2cart_y(rk[0], thetaEdges[-1])))
     
    # Untere Innere Tangente
    LModelSketch.Line(
        point1=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0)), 
        point2=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0) - ltan))
     
    # Untere äußere Tangente
    LowerOuterTan = LModelSketch.Line(
        point1=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0)), 
        point2=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0) - ltan))
    LowerOuterTan_id = LowerOuterTan.id
     
    # Obere äußere Tangente
    UpperOuterTan = LModelSketch.Line(
        point1=(ct.pol2cart_x(rk[-1], thetaEdges[-1]), ct.pol2cart_y(rk[-1], thetaEdges[-1])), 
        point2=(ct.pol2cart_x(rk[-1], thetaEdges[-1]) - ltan, ct.pol2cart_y(rk[-1], thetaEdges[-1])))
    UpperOuterTan_id = UpperOuterTan.id
     
    # Abschließende Linie oben
    LModelSketch.Line(
        point1=(ct.pol2cart_x(rk[-1], thetaEdges[-1]) - ltan, ct.pol2cart_y(rk[-1], thetaEdges[-1])), 
        point2=(ct.pol2cart_x(rk[0], thetaEdges[-1]) - ltan, ct.pol2cart_y(rk[0], thetaEdges[-1])))
    
    # Abschließende Linie unten
    LModelSketch.Line(
        point1=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0) - ltan), 
        point2=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0) - ltan))
     
    # Abschließen der Flanken von Radius
    LModelPart = LModel.Part(name='L-Profil-Part', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    LModelPart.BaseShell(sketch=LModelSketch)
     
    return (LModelPart, OuterCirc_id, LowerOuterTan_id, UpperOuterTan_id)


def sketchRigid():
    """
    Erstellt eine Skizze des Bolzens.

    Returns:
        Part: Die erstellte Part-Instanz des Bolzens.
    """
    RigidSketch = LModel.ConstrainedSketch(name='Rigid', sheetSize=10)
     
    # Zeichnet einen Kreis mit gegebenem Mittelpunkt und Perimeterpunkt
    RigidSketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, drigid / 2))
     
    # Erstellt ein Part-Objekt für den starren Körper
    RigidPart = LModel.Part(name='Rigid-Part', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)  # ANALYTIC_RIGID_SURFACE
     
    # Erzeugt die Basis-Shell aus der Skizze
    RigidPart.BaseShell(sketch=RigidSketch)
     
    return(RigidPart)

LModelPart,OuterCirc_id,LowerOuterTan_id,UpperOuterTan_id=sketchLModel()
#RigidPart=sketchRigid()
#--------------------------------------------------------------
#Abkürzungen als Variablen speichern
LModelAssembly = LModel.rootAssembly
LModelPart = mdb.models[modelName].parts['L-Profil-Part']
LModelFaces = LModelPart.faces

#Liste mit den Kanten der Flanken
flankenEdges = [0,-lb_transf + ra -d,-lt_transf + ra, -ltan]
def instanceLModel():
    """
    Erstellt Instanzen und positioniert sie entsprechend.

    Returns:
        Instance: Die letzte erstellte Instanz des starren Körpers.
    """
    LModelinstance = LModelAssembly.Instance(name='L_Profil_Instance', part=LModelPart, dependent=ON)
    #RigidPartinstance = LModelAssembly.Instance(name='Rigid-Part-1', part=RigidPart, dependent=ON)
    #LModelAssembly.translate(instanceList=('Rigid-Part-1',), vector=(flankenEdges[1], ct.pol2cart_y(rk[0], thetaEdges[-1]) - drigid / 2, 0))
    #RigidPartinstance = LModelAssembly.Instance(name='Rigid-Part-2', part=RigidPart, dependent=ON)
    #LModelAssembly.translate(instanceList=('Rigid-Part-2',), vector=(flankenEdges[2], ct.pol2cart_y(rk[-1], thetaEdges[-1]) + drigid / 2, 0))
    #RigidPartinstance = LModelAssembly.Instance(name='Rigid-Part-3', part=RigidPart, dependent=ON)
    #LModelAssembly.translate(instanceList=('Rigid-Part-3',), vector=(ct.pol2cart_x(rk[-1], thetaEdges[0]) + drigid / 2, flankenEdges[2], 0))
    #RigidPartinstance = LModelAssembly.Instance(name='Rigid-Part-4', part=RigidPart, dependent=ON)
    #LModelAssembly.translate(instanceList=('Rigid-Part-4',), vector=(ct.pol2cart_x(rk[0], thetaEdges[0]) - drigid / 2, flankenEdges[1], 0))
    
    return (LModelinstance)


#RigidPartInstance = instanceLModel()
#LModelInstance = LModelAssembly.instances['L_Profil_Instance']
LModelInstance = instanceLModel()
#--------------------------------------------------------------
#Partionierung

def partitionComposite():
    """
    Erstellt Partitionen im L-Profil-Teil basierend auf den definierten Geometrieparametern und den Lagenabständen.
    """
    LModelPartition = mdb.models[modelName].parts['L-Profil-Part']
    partitionfaces = LModelPartition.faces
    
    # Erstellen einer Skizzentransformation für die Partitionierung
    LModelTransform = LModelPartition.MakeSketchTransform(
        sketchPlane=partitionfaces[0], 
        sketchPlaneSide=SIDE1, 
        origin=(0, 0, 0))
    
    # Erstellen einer Skizze für die Partitionierung
    LModelPartitionSketch = mdb.models[modelName].ConstrainedSketch(
        name='L-Profil-Partition-Sketch', 
        sheetSize=(2 * (ltan + ra)), 
        transform=LModelTransform)
    g = LModelPartitionSketch.geometry
    
    # Projizieren von Referenzen auf die Skizze
    LModelPartition.projectReferencesOntoSketch(sketch=LModelPartitionSketch, filter=COPLANAR_EDGES)
    
    # Erstellen von Offsets und Partitionen für jede Schicht
    for ii in range(N - 1):
        LModelPartitionSketch.offset(
            distance=((ii + 1) * dsingle), 
            objectList=(
                g.findAt((ct.pol2cart_x(rk[0], thetaFaces[0]), ct.pol2cart_y(rk[0], thetaFaces[0]))),
                g.findAt((ct.pol2cart_x(rk[0], thetaEdges[0]), 0 - ltan / 2)),
                g.findAt((ct.pol2cart_x(rk[0], thetaEdges[-1]) - ltan / 2, ct.pol2cart_y(rk[0], thetaEdges[-1])))), side=RIGHT)
        
        LModelPickedFaces = partitionfaces.findAt((ct.pol2cart_x(rm[ii], thetaFaces[0]), ct.pol2cart_y(rm[ii], thetaFaces[0]), 0),)
        LModelPartition.PartitionFaceBySketch(faces=LModelPickedFaces, sketch=LModelPartitionSketch)

partitionComposite()
#--------------------------------------------------------------
#Materialzuweisung
def materialComposite():
    """
    Definiert die Materialeigenschaften für das Composite-Material und erstellt entsprechende Materialabschnitte.
    """
    # Composite-Material für 0-Grad Orientierung
    compositeMaterial = LModel.Material(compositeMaterialName + '_0')
    compositeMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=((E1, E2, E3, Nu12, Nu13, Nu23, G23, G13, G12), ))
    compositeMaterial.Expansion(type=ORTHOTROPIC, table=((alpha11, alpha22, alpha33), ))
    LModel.HomogeneousSolidSection(name='Composite_Section' + '_0', material=compositeMaterialName + '_0', thickness=None)
    
    # Composite-Material für 90-Grad Orientierung
    compositeMaterial = LModel.Material(compositeMaterialName + '_90')
    compositeMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=((E1, E3, E2, Nu13, Nu12, Nu23 * E3 / E2, G23, G12, G13), ))
    compositeMaterial.Expansion(type=ORTHOTROPIC, table=((alpha11, alpha33, alpha22), ))
    LModel.HomogeneousSolidSection(name='Composite_Section' + '_90', material=compositeMaterialName + '_90', thickness=None)

def materialRigid():
    """
    Definiert die Materialeigenschaften für das starre Material (C_10_Stahl) und weist es dem starren Teil zu.
    """
    # Definition des starren Materials
    mdb.models[modelName].Material(name='C_10_Stahl')
    mdb.models[modelName].materials['C_10_Stahl'].Elastic(table=((210000, 0.27), ))
    mdb.models[modelName].HomogeneousSolidSection(name='rigidsection', material='C_10_Stahl', thickness=None)
    
    # Zuweisung des starren Materials zum starren Teil
    LModelPart = mdb.models[modelName].parts['Rigid-Part']
    f = LModelPart.faces
    faces = f.findAt(((0.0, drigid / 4, 0.0), ))
    region = LModelPart.Set(faces=faces, name='Original-Bolt')
    LModelPart.SectionAssignment(region=region, sectionName='rigidsection', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

materialComposite()
#materialRigid()

#--------------------------------------------------------------
#Koordinatensystem festlegen
LModelPart.DatumCsysByThreePoints(name='Cartesian-KOS', coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
cartesianCompositeCSCylId = LModelPart.features['Cartesian-KOS'].id

#Zyl Coordinatensystem mit Ursprung im Centerpoint für Punkt Berechnung
#LModelPart.DatumCsysByThreePoints(origin=(ra,ra,0), point1=LModelPart.vertices.findAt(coordinates=(0, ra, 0.0)), point2=LModelPart.vertices.findAt(coordinates=(ra, 0.0, 0.0)), name='Cylindrical-KOS', coordSysType=CYLINDRICAL)
#LModelPart.DatumCsysByThreePoints(origin=(0,0,0), point1=LModelPart.vertices.findAt(coordinates=(0, ra, 0.0)), point2=LModelPart.vertices.findAt(coordinates=(ra, 0.0, 0.0)), name='Cylindrical-KOS', coordSysType=CYLINDRICAL)
LModelPart.DatumCsysByThreePoints(origin=(rk[-1],rk[-1],0), point1=(rk[-1], -ltan, 0.0), point2=(-ltan, rk[-1], 0.0), name='Cylindrical-KOS', coordSysType=CYLINDRICAL)
#LModelPart.DatumCsysByThreePoints(name='Cylindrical-KOS', coordSysType=CARTESIAN, origin=(ra, ra, 0.0), point1=(0.0, 0.0, 0.0), line2=(-ltan, ltan , 0.0))
#LModelPart.DatumCsysByThreePoints(name='Cylindrical-KOS', coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
curvedCompositeCSCylId = LModelPart.features['Cylindrical-KOS'].id
curvedCompositeCylCoordSys = LModel.rootAssembly.instances['L_Profil_Instance'].datums[curvedCompositeCSCylId]

#Zyl Coordinatensystem in 0 0 0
LModelPart.DatumCsysByThreePoints(origin=(0,0,0), point1=(1, 0 , 0.0), point2=(0, 1.0, 0.0), name='Cylindrical-KOS-Orientation', coordSysType=CYLINDRICAL)
orientationKOSid = LModelPart.features['Cylindrical-KOS-Orientation'].id


#KOS Für Boundary Contditions
LModelAssembly.DatumCsysByThreePoints(name='BoundaryCD-Datum', coordSysType=CARTESIAN, origin=(ra, ra, 0.0), 
                                      point1=(0.0, 0.0, 0.0), line2=(-ltan, ltan , 0.0))
BoundaryCDDatum_ID = LModelAssembly.features['BoundaryCD-Datum'].id
BoundaryCDDatum = LModelAssembly.datums[BoundaryCDDatum_ID]

#--------------------------------------------------------------

def hilfslinie():
    """
    Erstellt Hilfslinien und partitioniert das Modell am Auswärtungswinkel.
    """
    compositePartitionPlaneLaufVar = 1
     
    # Partitionierung in tangentialer Richtung
    for ii in range(len(thetaPartition)):
        LModelPartitionDatumPointLaufVar = len(LModelPart.datums)
         
        # Erstellen von Datumspunkten entlang der Partitionswinkel
        LModelPart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[0], thetaPartition[ii]), ct.pol2cart_y(rk[0], thetaPartition[ii]), length))
        LModelPart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1], thetaPartition[ii]), ct.pol2cart_y(rk[-1], thetaPartition[ii]), length))
        LModelPart.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1], thetaPartition[ii]), ct.pol2cart_y(rk[-1], thetaPartition[ii]), 0))
         
        compositePartDatumsKeys = list(LModelPart.datums.keys())
        compositePartDatumsKeys.sort()
         
        compositePartitionFace = LModelPart.faces
         
        # Erstellen einer Datumsebene durch die drei definierten Punkte
        LModelPart.DatumPlaneByThreePoints(
            point1=LModelPart.datums[compositePartDatumsKeys[LModelPartitionDatumPointLaufVar]],
            point2=LModelPart.datums[compositePartDatumsKeys[LModelPartitionDatumPointLaufVar + 1]],
            point3=LModelPart.datums[compositePartDatumsKeys[LModelPartitionDatumPointLaufVar + 2]])
         
        compositePartitionPlaneId = LModelPart.features['Datum plane-' + str(compositePartitionPlaneLaufVar)].id
         
        # Partitionieren der Flächen anhand der Datumsebene
        LModelPart.PartitionFaceByDatumPlane(faces=compositePartitionFace, datumPlane=LModelPart.datums[compositePartitionPlaneId])
         
        compositePartitionPlaneLaufVar += 1

hilfslinie()

#--------------------------------------------------------------
#Partitionierung der Flanken für besseres Mesh
def partitionflanken():
    """
    Partitioniert die Flanken des L-Profils für eine bessere Vernetzung.
     
    """
    LModelPart = mdb.models[modelName].parts['L-Profil-Part']
    LModelfaces = LModelPart.faces
    LModelTransform = LModelPart.MakeSketchTransform(
        sketchPlane=LModelfaces.findAt((ct.pol2cart_x(rm[0], thetaFaces[0]), ct.pol2cart_y(rm[0], thetaFaces[0]), 0.0), normal=(0.0, 0.0, 1.0)), 
        sketchPlaneSide=SIDE1, 
        origin=(0, 0, 0.0))
    LModelSketch = mdb.models[modelName].ConstrainedSketch(name='FlankenPartition', sheetSize=200, transform=LModelTransform)
    LModelgeom = LModelSketch.geometry
        # Partitionierung der Flanken
    for jj in range(len(thetaFaces)):
        LModelSketch.Line(
            point1=(ct.pol2cart_x(rk[0], thetaEdges[-jj]), ct.pol2cart_y(rk[0], thetaEdges[-jj])), 
            point2=(ct.pol2cart_x(rk[-1], thetaEdges[-jj]), ct.pol2cart_y(rk[-1], thetaEdges[-jj])))
        for ii in range(N):
            pickedFaces = LModelfaces.findAt(((ct.pol2cart_x(rm[ii], thetaFaces[-jj]), ct.pol2cart_y(rm[ii], thetaFaces[-jj]), 0.0),))
            LModelPart.PartitionFaceBySketch(faces=pickedFaces, sketch=LModelSketch)

partitionflanken()

#--------------------------------------------------------------
#Materialorientierung
def CompositeOrientation():
    """
    Teilt die Materialorientierung zu.
     
    Raises:
        ValueError: Wenn der plyAngle weder 0 noch 90 Grad beträgt.
    """
    # In Kurve
    for ii in range(N):
        for mm in range(len(thetaFaces)):
            curvedCompositeFaces = LModelPart.faces.findAt(((ct.pol2cart_x(rFacesAll[ii], thetaFaces[mm]), ct.pol2cart_y(rFacesAll[ii], thetaFaces[mm]), 0.0),))
            curvedCompositeRegion = regionToolset.Region(faces=curvedCompositeFaces)
            if plyAngle[ii] == 0:
                LModelPart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_' + str(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
            elif plyAngle[ii] == 90:
                LModelPart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_' + str(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
            else:
                raise ValueError('Please check the laminate layup. Only cross-plies allowed.')
            LModelPart.MaterialOrientation(region=curvedCompositeRegion, orientationType=SYSTEM, localCsys=LModelPart.datums[orientationKOSid], axis=AXIS_3, additionalRotationType=ROTATION_ANGLE, angle=0, stackDirection=STACK_3)
    
    # In Flanke unten
    for ii in range(N):
        curvedCompositeFaces = LModelPart.faces.findAt(((ct.pol2cart_x(rFacesAll[ii], thetaEdges[0]), -ltan/2, 0.0),))
        curvedCompositeRegion = regionToolset.Region(faces=curvedCompositeFaces)
        if plyAngle[ii] == 0:
            LModelPart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_' + str(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
        elif plyAngle[ii] == 90:
            LModelPart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_' + str(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
        else:
            raise ValueError('Please check the laminate layup. Only cross-plies allowed.')
        LModelPart.MaterialOrientation(region=curvedCompositeRegion, orientationType=SYSTEM, localCsys=LModelPart.datums[cartesianCompositeCSCylId], axis=AXIS_3, additionalRotationType=ROTATION_NONE, angle=0, stackDirection=STACK_3)
    
    # In Flanke oben    
    for ii in range(N):
        curvedCompositeFaces = LModelPart.faces.findAt(((-ltan/2, ct.pol2cart_y(rFacesAll[ii], thetaEdges[-1]), 0.0),))
        curvedCompositeRegion = regionToolset.Region(faces=curvedCompositeFaces)
        if plyAngle[ii] == 0:
            LModelPart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_' + str(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
        elif plyAngle[ii] == 90:
            LModelPart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_' + str(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
        else:
            raise ValueError('Please check the laminate layup. Only cross-plies allowed.')
        LModelPart.MaterialOrientation(region=curvedCompositeRegion, orientationType=SYSTEM, localCsys=LModelPart.datums[cartesianCompositeCSCylId], axis=AXIS_3, additionalRotationType=ROTATION_ANGLE, angle=90, stackDirection=STACK_3)

CompositeOrientation()

#--------------------------------------------------------------
#Partitionierung an Kontaktstelle für Besseres Mesh und mesh verlauf
def partitionContact():
    """
    Partitioniert die Kontaktstellen für ein besseres Mesh und Mesh-Verlauf.
    """
    LModelTransform = LModelPart.MakeSketchTransform(
        sketchPlane=LModelFaces.findAt((ct.pol2cart_x(rm[0], thetaFaces[0]), ct.pol2cart_y(rm[0], thetaFaces[0]), 0.0), normal=(0.0, 0.0, 1.0)), 
        sketchPlaneSide=SIDE1, 
        origin=(0, 0, 0.0))
    LModelSketch = mdb.models[modelName].ConstrainedSketch(name='ContactPartition', sheetSize=200, transform=LModelTransform)
    
    # Partitionierung der Flanken
    # untere Flanke 
    LModelSketch.Line(point1=(ri, flankenEdges[1]), point2=(ra, flankenEdges[1]))
    LModelSketch.Line(point1=(ri, flankenEdges[2]), point2=(ra, flankenEdges[2]))
    
    # obere Flanke   
    LModelSketch.Line(point1=(flankenEdges[1], ri), point2=(flankenEdges[1], ra))
    LModelSketch.Line(point1=(flankenEdges[2], ri), point2=(flankenEdges[2], ra))
    
    for ii in range(N):
        pickedFaces = LModelFaces.findAt(((-ltan / 2, ct.pol2cart_y(rm[ii], thetaEdges[-1]), 0.0),))
        LModelPart.PartitionFaceBySketch(faces=pickedFaces, sketch=LModelSketch)
        pickedFaces = LModelFaces.findAt(((ct.pol2cart_x(rm[ii], thetaEdges[0]), -ltan / 2, 0.0),))
        LModelPart.PartitionFaceBySketch(faces=pickedFaces, sketch=LModelSketch)

#partitionContact()


#--------------------------------------------------------------
#Meshgenerierung
flankenFaces = []
for ii in range(1,6,2):
    flankenFaces.append(-float(ltan)*ii/6)

def meshLPart(mR, mT, mFlanken, mBias, mRRatio, mTRatio):
    """
    Erstellt ein Netz für das L-Profil unter Verwendung der angegebenen Parameter.
     
    Args:
        mR (int): Anzahl der Elemente in radialer Richtung.
        mT (int): Anzahl der Elemente in tangentialer Richtung.
        mFlanken (int): Anzahl der Elemente entlang der Flanken.
        mBias (float): Bias-Wert für das Netz.
        mRRatio (float): Verhältnis für das Bias im radialen Netz.
        mTRatio (float): Verhältnis für das Bias im tangentialen Netz.
     
    Raises:
        ValueError: Wenn eine Kante keine Richtung hat.
    """
    for ii in range(len(rFacesAll)):
        for jj in range(len(thetaEdges)):
            curvedCompositeMeshEdge = LModelPart.edges.findAt(
                ((ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj]),
                  ct.pol2cart_y(rFacesAll[ii], thetaEdges[jj]), 0.0),)
            )
            LModelPart.seedEdgeByBias(
                biasMethod=DOUBLE,
                endEdges=curvedCompositeMeshEdge,
                ratio=mRRatio,
                number=mR,
                constraint=FINER
            )
     
    for ii in range(len(rk)):
          curvedCompositeMeshEdge = LModelPart.edges.findAt(
                ((ct.pol2cart_x(rk[ii], thetaFaces[-1]),
                  ct.pol2cart_y(rk[ii], thetaFaces[-1]), 0.0),)
            )
          LModelPart.seedEdgeByBias(
            biasMethod=SINGLE,
            end2Edges=curvedCompositeMeshEdge,
            ratio=mTRatio,
            number=mT/2,
            constraint=FINER
            ) 
     
    for ii in range(len(rk)):
          curvedCompositeMeshEdge = LModelPart.edges.findAt(
                ((ct.pol2cart_x(rk[ii], thetaFaces[0]),
                  ct.pol2cart_y(rk[ii], thetaFaces[0]), 0.0),)
            )
          LModelPart.seedEdgeByBias(
            biasMethod=SINGLE,
            end1Edges=curvedCompositeMeshEdge,
            ratio=mTRatio,
            number=mT/2,
            constraint=FINER
            )   
     
    for ii in range(len(rFacesAll)):
        for jj in range(len(flankenEdges)):
            curvedCompositeMeshEdge = LModelPart.edges.findAt(
                ((ct.pol2cart_x(rm[ii], thetaEdges[0]),
                  flankenEdges[jj], 0.0),)
            )
            LModelPart.seedEdgeByBias(
                biasMethod=DOUBLE,
                endEdges=curvedCompositeMeshEdge,
                ratio=mRRatio,
                number=mR,
                constraint=FINER
            )
     
    for jj in range(len(flankenFaces)):
        for ii in range(len(rk)):
            CurvedCompositeMeshEdge = LModelPart.edges.findAt(
                ((rk[ii], flankenFaces[jj], 0.0),)
            )
            LModelPart.seedEdgeByBias(
                biasMethod=DOUBLE,
                endEdges=CurvedCompositeMeshEdge,
                ratio=mBias,
                number=mFlanken,
                constraint=FINER
            )
     
    for jj in range(len(flankenEdges)):
        for ii in range(len(rFacesAll)):
            curvedCompositeMeshEdge = LModelPart.edges.findAt(
                ((flankenEdges[jj],
                  ct.pol2cart_y(rm[ii], thetaEdges[-1]), 0.0),)
            )
            LModelPart.seedEdgeByBias(
                biasMethod=DOUBLE,
                endEdges=curvedCompositeMeshEdge,
                ratio=mRRatio,
                number=mR,
                constraint=FINER
            )
     
    for jj in range(len(flankenFaces)):
        for ii in range(len(rk)):
            CurvedCompositeMeshEdge = LModelPart.edges.findAt(
                ((flankenFaces[jj], rk[ii], 0.0),)
            )
            LModelPart.seedEdgeByBias(
                biasMethod=DOUBLE,
                endEdges=CurvedCompositeMeshEdge,
                ratio=mBias,
                number=mFlanken,
                constraint=FINER
            )
     
    compositeElementType = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
    LModelPart.setMeshControls(
        regions=LModelPart.faces,
        elemShape=HEX,
        technique=STRUCTURED
    )
    LModelPart.setElementType(
        regions=(LModelPart.faces[:],),
        elemTypes=(compositeElementType,)
    )
    LModelPart.seedPart(size=0.05, deviationFactor=0.1, constraint=FINER)
    LModelPart.generateMesh()

meshLPart(mR, mT, mFlanken, mBias, mRRatio, mTRatio)


def meshrigid(mRigid):
    """
    Erstellt die Vernetzung der Bolzen.

    Args:
        mRigid (float): Größe der Elemente im Netz.
    """
    RigidPart = mdb.models[modelName].parts['Rigid-Part']
    RigidPart.seedPart(size=mRigid, deviationFactor=0.1, minSizeFactor=0.1)
    pickedRegions = RigidPart.faces.findAt(((0.0, drigid / 2, 0.0),))
    RigidPart.setMeshControls(regions=pickedRegions, algorithm=MEDIAL_AXIS)
    RigidPart.generateMesh()

#meshrigid(mRigid)

#--------------------------------------------------------------
#Sets zur Auswärtung

if Null_Auswaertung == True:
    angleEval_FEInnner = 0
else:
    angleEval_FEInnner =angleEval

def setsCurvedCompositeCirc(rFacesAll, angleEval_FEInnner, iInterfaceEval, thetaPartitionFaces, angleEval, rk):
    """
    Erstellt Sets zur Verifizierung der Inner-Solution und zur Bestimmung des Abklingverhaltens der singulären Spannungen.

    Args:
        rFacesAll (list): Liste aller Radien der Flächen.
        angleEval_FEInnner (float): Winkel zur Auswertung der inneren Lösung.
        iInterfaceEval (list): Liste der auszuwertenden Interfaces.
        thetaPartitionFaces (list): Liste der Winkelpartitionen.
        angleEval (float): Winkel zur Auswertung.
        rk (list): Liste der Radien.
    """
    # Sets zur Verifizierung der Inner-Solution:
    # Sets zur Verifizierung der Inner-Solution:
    for ii in range(len(rFacesAll)):
        compositeSetEvalEdge = LModelInstance.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], angleEval_FEInnner), ct.pol2cart_y(rFacesAll[ii], angleEval_FEInnner), 0.0),))
        LModelAssembly.Set(edges=compositeSetEvalEdge, name='FEInner' + str(ii + 1))
    # Sets zur Bestimmung des Abklingverhaltens der singulären Spannungen (in Dickenrichtung) in den Auswertungsinterfaces an beiden Enden:
    for kk in range(len(iInterfaceEval)):
        curvedCompositeStringEdgesRight = ''
        ii = iInterfaceEval[kk]
        for jj in range(len(thetaPartitionFaces) / 2):
            curvedCompositeStringEdgesRight = curvedCompositeStringEdgesRight + str(((ct.pol2cart_x(rk[ii], angleEval_FEInnner - 0.00001), ct.pol2cart_y(rk[ii], angleEval_FEInnner - 0.00001), 0.0),),) + ','
        curvedCompositeStringEdgesRightExec = 'compositeSetEdgeRight = LModelInstance.edges.findAt(' + curvedCompositeStringEdgesRight + ')'
        exec(curvedCompositeStringEdgesRightExec)
        LModelAssembly.Set(edges=compositeSetEdgeRight, name='FEInterfaceCircRight' + str(iInterfaceEval[kk]))
    for kk in range(len(iInterfaceEval)):
        curvedCompositeStringEdgesLeft = ''
        ii = iInterfaceEval[kk]
        for jj in range(len(thetaPartitionFaces) / 2):
            curvedCompositeStringEdgesLeft = curvedCompositeStringEdgesLeft + str(((ct.pol2cart_x(rk[ii], angleEval - 0.00001), ct.pol2cart_y(rk[ii], angleEval - 0.00001), 0.0),),) + ','
        curvedCompositeStringEdgesLeftExec = 'compositeSetEdgeLeft = LModelInstance.edges.findAt(' + curvedCompositeStringEdgesLeft + ')'
        exec(curvedCompositeStringEdgesLeftExec)
        LModelAssembly.Set(edges=compositeSetEdgeLeft, name='FEInterfaceCircLeft' + str(iInterfaceEval[kk]))

setsCurvedCompositeCirc(rFacesAll, angleEval_FEInnner, iInterfaceEval, thetaPartitionFaces, angleEval, rk)


#--------------------------------------------------------------
#def stepInitialise():
    #"""
    #Erstellt einen quasi-statischen Schritt.

    #Returns:
    #    str: Der Name des erstellten Schritts.
    #"""
    #step = 'BendingMoment'
    #LModel.StaticStep(name=step, previous='Initial', description='Introduce bending moment', nlgeom=ON, initialInc=0.1)
    #return step

#step = stepInitialise()



def stepBending():
	stepPrevious,step,stepDescription = 'Initial','Mechanical_Loading','Introduce mechanical loadings'
	LModel.StaticStep(name=step, previous=stepPrevious, description=stepDescription, nlgeom=OFF, initialInc=1.0)
	return (stepPrevious,step)

bendingMoment = -1.0
stepPrevious,step = stepBending()

def loadsCreateKinCoup(jj):
	# Erstellen eines Referenzpunktes und Fixierung dessen:
	LModelAssembly.ReferencePoint(point=((
		ct.pol2cart_x(rk[N//2], thetaEdges[jj-1]) + (jj-1)*ltan, 
		ct.pol2cart_y(rk[N//2], thetaEdges[jj-1]) - (jj)*ltan, 
		0)))
	
	curvedCompositeLoadRefPoint = (
		LModelAssembly.referencePoints.findAt((
			ct.pol2cart_x(rk[N//2], thetaEdges[jj-1]) + (jj-1)*ltan, 
			ct.pol2cart_y(rk[N//2], thetaEdges[jj-1]) - (jj)*ltan, 
			0),),)
	
	curvedCompositeLoadRefPointRegion = regionToolset.Region(referencePoints=curvedCompositeLoadRefPoint)
	
	curvedCompositeFixRefPointSet = LModelAssembly.Set(
		region=curvedCompositeLoadRefPointRegion, 
		name='Set_Fix_RefPoint_Pos_' + str(jj*int(180*angleOpening/np.pi)))
	
	#LModel.DisplacementBC(name='Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeFixRefPointSet, localCsys=curvedCompositeCylCoordSys, u3=0.0, ur1=0.0, ur2=0.0)
	LModel.DisplacementBC(name='Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeFixRefPointSet, u3=0.0, ur1=0.0, ur2=0.0)
	
	#LModel.DisplacementBC(
		#name='Fix_RefPoint_Pos_' + str(jj*int(180*angleOpening/np.pi)), 
		#createStepName=step, 
		#region=curvedCompositeFixRefPointSet, 
		#localCsys=LModelPart.datums[cartesianCompositeCSCylId], 
		#u3=0.0, ur1=0.0, ur2=0.0)
    
	# Kinematic Coupling:
	curvedCompositeStringKinCoupFaces =''
	for ii in range(len(rFacesAll)):
		curvedCompositeStringKinCoupFaces = curvedCompositeStringKinCoupFaces + str(((ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj-1]) + (jj-1)*ltan, ct.pol2cart_y(rFacesAll[ii], thetaEdges[jj-1]) - (jj)*ltan, 0.0),),) + ','
	curvedCompositeKinCoupFacesExec = 'curvedCompositeKinCoupFaces = LModelInstance.edges.findAt(' + curvedCompositeStringKinCoupFaces + ')'
	exec(curvedCompositeKinCoupFacesExec)
	curvedCompositeKinCoupRegion = regionToolset.Region(edges=curvedCompositeKinCoupFaces)
	#LModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion,influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=OFF, u2=ON, ur3=ON)
	if jj == 0:
		LModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion,influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=ON, u2=OFF, ur3=ON)
	else:
		LModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion,influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=OFF, u2=ON, ur3=ON)        
	#LModel.constraints['Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi))].setValues(localCsys=curvedCompositeCylCoordSys)
	LModel.constraints['Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi))].setValues()
	return curvedCompositeLoadRefPointRegion

def loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion):
	# Aeussere Last - Moment:
	#LModel.Moment(name='Load_BendingMoment_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeLoadRefPointRegion, cm3=((-1)**(jj))*bendingMoment, localCsys=curvedCompositeCylCoordSys, distributionType=UNIFORM)
	LModel.Moment(name='Load_BendingMoment_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeLoadRefPointRegion, cm3=((-1)**(jj))*bendingMoment, distributionType=UNIFORM)

noKinCoup = 2
for jj in range(noKinCoup):
	curvedCompositeLoadRefPointRegion = loadsCreateKinCoup(jj)
	loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion)


# ---------------------------------------------------------------
# Field Output:
#LModel.FieldOutputRequest(name='F-Output-1', createStepName='Initial', variables=('S', 'U', 'COORD'))

LModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'COORD'))

#-------------------------------------------------------------------
# Speichern der CAE-Datei:
mdb.saveAs(pathName=analysis_newPath + '\\' + modelName + '_CAE_Model')

#--------------------------------------------------------------
#Job
jobName = modelName + '_Job'

def LmodelJob():
    """
    Erstellt und führt einen Job für das Modell aus.

    Returns:
        float: Simulationszeit in Sekunden.
    """
    global startTime
    mdb.Job(
        name=jobName, 
        model=modelName, 
        description='', 
        type=ANALYSIS, 
        atTime=None, 
        waitMinutes=0, 
        waitHours=0, 
        queue=None, 
        memory=90, 
        memoryUnits=PERCENTAGE, 
        getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, 
        echoPrint=OFF, 
        modelPrint=OFF, 
        contactPrint=OFF, 
        historyPrint=OFF, 
        userSubroutine='', 
        scratch='', 
        resultsFormat=ODB, 
        multiprocessingMode=DEFAULT, 
        numCpus=6, 
        numDomains=6
    )
    startTime = time.time()
    localTime = time.strftime("%H:%M", time.localtime(startTime))
    print('Start:' + localTime + ' Uhr')
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()
    simTime = time.time() - startTime
    return(simTime)

simTime = LmodelJob()

#-------------------------------------------------------------------
#+++++++++++++++++++++++++ Post- Processing ++++++++++++++++++++++++
#------------------------------------------------------------------- 

elementsLoc = [['r < r_eval','r > r_eval'],['t < t_eval','t > t_eval']]
elementNodesLoc = [['first three element nodes','last three element nodes'],['r < r_eval','r > r_eval'],['t < t_eval','t > t_eval']]



# ---------------------------------------------------------------
#---------------------------------------------------------------
# Ausgabe der extrapolierten Spannungen an den jeweiligen Elementknoten des ersten und letzten Elements des Sets (in Breitenrichtung):
def getDispComp2D():
    return('x','y')

def getStressComp2D():
    return('Sigma_rr','Sigma_tt','Sigma_zz','Tau_rt')

def getStressesCoordinatesNode(step,FieldValueObject):
    stressComp = {stressCompStr:FieldValueObject.data[ii] for ii,stressCompStr in enumerate(getStressComp2D())}
    dispComp = {dispCompStr:compositeModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].\
                instance.getNodeFromLabel(FieldValueObject.nodeLabel).coordinates[ii] \
                for ii,dispCompStr in enumerate(getDispComp2D())}
    rCoord = ct.cart2pol_radius(dispComp['x'],dispComp['y'])
    tCoord = (180*(ct.cart2pol_theta(dispComp['x'],dispComp['y']))/np.pi)
    
    return(stressComp[getStressComp2D()[0]], stressComp[getStressComp2D()[1]], stressComp[getStressComp2D()[2]], stressComp[getStressComp2D()[3]], dispComp[getDispComp2D()[0]], dispComp[getDispComp2D()[1]],rCoord, tCoord, FieldValueObject.nodeLabel, FieldValueObject.elementLabel)

# Ausgabe der Spannungen in den Integrationspunkten der jeweiligen Elemente des Sets:
def getStressesIntegrationPoint(FieldValueObject):
    stressComp = {stressCompStr:FieldValueObject.data[ii] for ii,stressCompStr in enumerate(getStressComp2D())}
    
    return (FieldValueObject.elementLabel, FieldValueObject.integrationPoint,\
            stressComp[getStressComp2D()[0]], stressComp[getStressComp2D()[1]], stressComp[getStressComp2D()[2]], stressComp[getStressComp2D()[3]])

#---------------------------------------------------------------
# Funktionen zur Anwendung von Sortieralgorithmen:
def getCoordinateR(unsortedElementNodes):
    return ct.cart2pol_radius(unsortedElementNodes[1],unsortedElementNodes[2])

def getCoordinateX(cart_unsortedElementNodes):
    return cart_unsortedElementNodes[1]

def getCoordinateT(unsortedElementNodes):
    return ct.cart2pol_theta(unsortedElementNodes[1],unsortedElementNodes[2])

def getCoordinateY(cart_unsortedElementNodes):
    return cart_unsortedElementNodes[2]

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
    compositeSetElements=compositeModelOdbObject.rootAssembly.elementSets[set]
    compositeSetElementNodes=compositeModelOdbObject.rootAssembly.nodeSets[set]
    
    # Bestimmung des letzten Frames der FE-Analyse:
    indexLastFrame = len(compositeModelOdbObject.steps[step].frames) - 1
    
    # Bestimmung der Knoten-Koordinaten/Knoten-Verschiebungen/Element-Spannungen im zylindrischen Koordinatensystem und abspeichern dieser in Arrays:
    coordsFirstFrame = compositeModelOdbObject.steps[step].frames[0].fieldOutputs['COORD'].getSubset(region=compositeSetElementNodes).values
    dispLastFrame = compositeModelOdbObject.steps[step].frames[indexLastFrame].fieldOutputs['U'].getSubset(region=compositeSetElementNodes).getTransformedField(datumCsys=postProc_Cyl_CS).values
    stressLastFrameElementNodalTemp = compositeModelOdbObject.steps[step].frames[-1].fieldOutputs['S'].getSubset(position=ELEMENT_NODAL, region=compositeSetElements).getTransformedField(datumCsys=postProc_Cyl_CS).values
    stressLastFrameGaussPoint = compositeModelOdbObject.steps[step].frames[-1].fieldOutputs['S'].getSubset(position=INTEGRATION_POINT, region=compositeSetElements).getTransformedField(datumCsys=postProc_Cyl_CS).values
    
    # Knotenkoordinaten des Knotensets in aufsteigender Reihenfolge sortieren - freier Rand (z = 0) bis entsprechend Ende (Knotenkoordinaten des letzten Elementknotens):
    unsortedElementNodes = []
    for ii in range(len(dispLastFrame)):
        unsortedElementNodes.append((dispLastFrame[ii].nodeLabel,coordsFirstFrame[ii].data[0],coordsFirstFrame[ii].data[1]))
    
    if sortDirection == 'r':
        sortedElementNodes = sorted(unsortedElementNodes, key=getCoordinateR)
    elif sortDirection == 't':
        sortedElementNodes = sorted(unsortedElementNodes, key=getCoordinateT)
    
    return(compositeSetElements,stressLastFrameElementNodalTemp,stressLastFrameGaussPoint,sortedElementNodes)

def cart_getStressesDisplacements(set,step,sortDirection):
    # Einlesen des im Preprocessing definierten Knoten-/Element-Sets:
    compositeSetElements=compositeModelOdbObject.rootAssembly.elementSets[set]
    compositeSetElementNodes=compositeModelOdbObject.rootAssembly.nodeSets[set]
    
    # Bestimmung des letzten Frames der FE-Analyse:
    indexLastFrame = len(compositeModelOdbObject.steps[step].frames) - 1
    
    # Bestimmung der Knoten-Koordinaten/Knoten-Verschiebungen/Element-Spannungen im zylindrischen Koordinatensystem und abspeichern dieser in Arrays:
    coordsFirstFrame = compositeModelOdbObject.steps[step].frames[0].fieldOutputs['COORD'].getSubset(region=compositeSetElementNodes).values
    cart_dispLastFrame = compositeModelOdbObject.steps[step].frames[indexLastFrame].fieldOutputs['U'].getSubset(region=compositeSetElementNodes).values
    cart_stressLastFrameElementNodalTemp = compositeModelOdbObject.steps[step].frames[-1].fieldOutputs['S'].getSubset(position=ELEMENT_NODAL, region=compositeSetElements).values
    stressLastFrameGaussPoint = compositeModelOdbObject.steps[step].frames[-1].fieldOutputs['S'].getSubset(position=INTEGRATION_POINT, region=compositeSetElements).getTransformedField(datumCsys=postProc_Cyl_CS).values
    
    # Knotenkoordinaten des Knotensets in aufsteigender Reihenfolge sortieren - freier Rand (z = 0) bis entsprechend Ende (Knotenkoordinaten des letzten Elementknotens):  
    cart_unsortedElementNodes = []
    for ii in range(len(cart_dispLastFrame)):
        cart_unsortedElementNodes.append((cart_dispLastFrame[ii].nodeLabel,coordsFirstFrame[ii].data[0],coordsFirstFrame[ii].data[1]))
    
    if sortDirection == 'r':
        cart_sortedElementNodes = sorted(cart_unsortedElementNodes, key=getCoordinateX)
        
    elif sortDirection == 't':
        cart_sortedElementNodes = sorted(cart_unsortedElementNodes, key=getCoordinateY)
    
    return(compositeSetElements,cart_stressLastFrameElementNodalTemp,stressLastFrameGaussPoint,cart_sortedElementNodes)

def elemLocCol(colDir):
    if colDir == 'r':
        elements = {elementsLocT:[] for elementsLocT in elementsLoc[1]}
        elementLabels = {elementsLocT:[] for elementsLocT in elementsLoc[1]}
    elif colDir == 't':
        elements = {elementsLocR:[] for elementsLocR in elementsLoc[0]}
        elementLabels = {elementsLocR:[] for elementsLocR in elementsLoc[0]}
    return(elements,elementLabels)

def elemNodeLocCol(colDir):
    if colDir == 'r':
        stressesElementNodes = {elementNodesLocR:{elementNodesLocT:[] for elementNodesLocT in elementNodesLoc[2]} for elementNodesLocR in elementNodesLoc[0]}
    elif colDir == 't':
        stressesElementNodes = {elementNodesLocT:{elementNodesLocR:[] for elementNodesLocR in elementNodesLoc[1]} for elementNodesLocT in elementNodesLoc[0]}
    return stressesElementNodes

def elemGaussPointsLocCol(colDir):
    if colDir == 'r':
        stressGaussPoints = {elementsLocT:[] for elementsLocT in elementsLoc[1]}
    elif colDir == 't':
        stressGaussPoints = {elementsLocR:[] for elementsLocR in elementsLoc[0]}
    return stressGaussPoints

def dataGaussPointsLocCol(colDir):
    if colDir == 'r':
        dataOutputGaussPoints = {elementsLocT:[] for elementsLocT in elementsLoc[1]}
    elif colDir == 't':
        dataOutputGaussPoints = {elementsLocR:[] for elementsLocR in elementsLoc[0]}
    return dataOutputGaussPoints

def stressElementNodesF3NRadialInner(stressesElementNodes,step,elements,conVar):
    try:
        for ii in range(np.array(elementsLoc).shape[1]):
            getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][conVar])
        for ii in range(np.array(elementsLoc).shape[1]):
            stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][ii]].append((getStressesCoordinatesNode(step,elements[elementsLoc[1][ii]][conVar])))
    except:
        raise KeyError
    return(stressesElementNodes)

def stressElementNodesF3N(stressesElementNodes,step,elements,conVar,colDir):
    if colDir == 'r':
        try:
            stressesElementNodes = stressElementNodesF3NRadialInner(stressesElementNodes,step,elements,conVar)
        except:
            raise KeyError("Problem in method 'stressElementNodesF3N'")
    elif colDir == 't':
        for ii in range(np.array(elementsLoc).shape[1]):
            stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[1][ii]].append((getStressesCoordinatesNode(step,elements[elementsLoc[0][ii]][conVar])))
    return(stressesElementNodes)

def elemDefLocCenGrav(rEval, tEval, stressLastFrameElementNodal,compositeSetElements):
    # Berechnung des Schwerpunktes aller Elemente des jeweiligen Elementknotens des Sets und aufgrund der Position, Abspeichern dieser in den jeweiligen Arrays:
    if rEval == 'None':
        colDir = 'r'
    else:
        colDir = 't'
    elements, elementLabels = elemLocCol(colDir)
    elementsExist = []
    if colDir == 'r':
        for ii in range(len(stressLastFrameElementNodal)):
            # Ueberpruefung, ob die Koordinaten des Elementschwerpunktes bereits berechnet wurden:
            if stressLastFrameElementNodal[ii].elementLabel in elementsExist:
                for jj in range(np.array(elementsLoc).shape[1]):
                    if stressLastFrameElementNodal[ii].elementLabel in elementLabels[elementsLoc[1][jj]]:
                        elements[elementsLoc[1][jj]].append(stressLastFrameElementNodal[ii])
            else:
                elementsExist.append(stressLastFrameElementNodal[ii].elementLabel)
                # Identifikations des Elementes, dessen Schwerpunkt berechnet werden soll:
                for jj in range(len(compositeSetElements.elements[0])):
                    if compositeSetElements.elements[0][jj].label == stressLastFrameElementNodal[ii].elementLabel:
                        # Betrachte alle Elementknoten des jeweiligen Elements, bestimme mithilfe der kartesischen Koordinaten ueber eine Mittelung die Radien/Winkel der Gauss-Punkte
                        # und speichere Elemente, deren Positionen der Gauss-Punkte gewuenschte Kriterien erfuellen, zur Weiterverarbeitung ab:
                        nodeLabelConnect = np.array(list(compositeSetElements.elements[0][jj].connectivity))
                        nodeLabelConnectInst = compositeModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].instance
                        xCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[0] for kk in nodeLabelConnect])
                        yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
                        
                        tCoord = ct.cart2pol_theta(xCoord,yCoord)
                        sumNodesPosThetaMean = 180*np.mean(tCoord)/np.pi
                        
                        # Abspeichern der Element-Labels der Elemente, welche die gewuenschten Kriterien erfuellen:
                        if sumNodesPosThetaMean > tEval:
                            elements['t > t_eval'].append(stressLastFrameElementNodal[ii])
                            elementLabels['t > t_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
                        else:
                            elements['t < t_eval'].append(stressLastFrameElementNodal[ii])
                            elementLabels['t < t_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
    elif colDir == 't':
        for ii in range(len(stressLastFrameElementNodal)):
            # Ueberpruefung, ob die Koordinaten des Elementschwerpunktes bereits berechnet wurden:
            if stressLastFrameElementNodal[ii].elementLabel in elementsExist:
                for jj in range(np.array(elementsLoc).shape[1]):
                    if stressLastFrameElementNodal[ii].elementLabel in elementLabels[elementsLoc[0][jj]]:
                        elements[elementsLoc[0][jj]].append(stressLastFrameElementNodal[ii])
            else:
                elementsExist.append(stressLastFrameElementNodal[ii].elementLabel)
                # Identifikations des Elementes, dessen Schwerpunkt berechnet werden soll:
                for jj in range(len(compositeSetElements.elements[0])):
                    if compositeSetElements.elements[0][jj].label == stressLastFrameElementNodal[ii].elementLabel:
                        # Betrachte alle Elementknoten des jeweiligen Elements, bestimme mithilfe der kartesischen Koordinaten ueber eine Mittelung die Radien/Winkel der Gauss-Punkte
                        # und speichere Elemente, deren Positionen der Gauss-Punkte gewuenschte Kriterien erfuellen, zur Weiterverarbeitung ab:
                        nodeLabelConnect = np.array(list(compositeSetElements.elements[0][jj].connectivity))
                        nodeLabelConnectInst = compositeModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].instance
                        xCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[0] for kk in nodeLabelConnect])
                        yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
                        
                        rCoord = ct.cart2pol_radius(xCoord,yCoord)
                        sumNodesPosRMean  = np.mean(rCoord)
                        
                        # Abspeichern der Element-Labels der Elemente, welche die gewuenschten Kriterien erfuellen:
                        if sumNodesPosRMean > rEval:
                            elements['r > r_eval'].append(stressLastFrameElementNodal[ii])
                            elementLabels['r > r_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
                        else:
                            elements['r < r_eval'].append(stressLastFrameElementNodal[ii])
                            elementLabels['r < r_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
    return(elements, elementLabels)

def cart_elemDefLocCenGrav(rEval, tEval, stressLastFrameElementNodal,compositeSetElements):
    # Berechnung des Schwerpunktes aller Elemente des jeweiligen Elementknotens des Sets und aufgrund der Position, Abspeichern dieser in den jeweiligen Arrays:
    if rEval == 'None':
        colDir = 'r'
    else:
        colDir = 't'
    elements, elementLabels = elemLocCol(colDir)
    elementsExist = []
    if colDir == 'r':#not adjusted, irrelevant
        for ii in range(len(stressLastFrameElementNodal)):
            # Ueberpruefung, ob die Koordinaten des Elementschwerpunktes bereits berechnet wurden:
            if stressLastFrameElementNodal[ii].elementLabel in elementsExist:
                for jj in range(np.array(elementsLoc).shape[1]):
                    if stressLastFrameElementNodal[ii].elementLabel in elementLabels[elementsLoc[1][jj]]:
                        elements[elementsLoc[1][jj]].append(stressLastFrameElementNodal[ii])
            else:
                elementsExist.append(stressLastFrameElementNodal[ii].elementLabel)
                # Identifikations des Elementes, dessen Schwerpunkt berechnet werden soll:
                for jj in range(len(compositeSetElements.elements[0])):
                    if compositeSetElements.elements[0][jj].label == stressLastFrameElementNodal[ii].elementLabel:
                        # Betrachte alle Elementknoten des jeweiligen Elements, bestimme mithilfe der kartesischen Koordinaten ueber eine Mittelung die Radien/Winkel der Gauss-Punkte
                        # und speichere Elemente, deren Positionen der Gauss-Punkte gewuenschte Kriterien erfuellen, zur Weiterverarbeitung ab:
                        nodeLabelConnect = np.array(list(compositeSetElements.elements[0][jj].connectivity))
                        nodeLabelConnectInst = compositeModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].instance
                        xCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[0] for kk in nodeLabelConnect])
                        yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
                        
                        tCoord = ct.cart2pol_theta(xCoord,yCoord)
                        sumNodesPosThetaMean = 180*np.mean(tCoord)/np.pi
                        
                        # Abspeichern der Element-Labels der Elemente, welche die gewuenschten Kriterien erfuellen:
                        if sumNodesPosThetaMean > tEval:
                            elements['t > t_eval'].append(stressLastFrameElementNodal[ii])
                            elementLabels['t > t_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
                        else:
                            elements['t < t_eval'].append(stressLastFrameElementNodal[ii])
                            elementLabels['t < t_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
    elif colDir == 't':
        for ii in range(len(stressLastFrameElementNodal)):
            # Ueberpruefung, ob die Koordinaten des Elementschwerpunktes bereits berechnet wurden:
            if stressLastFrameElementNodal[ii].elementLabel in elementsExist:
                for jj in range(np.array(elementsLoc).shape[1]):
                    if stressLastFrameElementNodal[ii].elementLabel in elementLabels[elementsLoc[0][jj]]:
                        elements[elementsLoc[0][jj]].append(stressLastFrameElementNodal[ii])
            else:
                elementsExist.append(stressLastFrameElementNodal[ii].elementLabel)
                # Identifikations des Elementes, dessen Schwerpunkt berechnet werden soll:
                for jj in range(len(compositeSetElements.elements[0])):
                    if compositeSetElements.elements[0][jj].label == stressLastFrameElementNodal[ii].elementLabel:
                        # Betrachte alle Elementknoten des jeweiligen Elements, bestimme mithilfe der kartesischen Koordinaten ueber eine Mittelung die Radien/Winkel der Gauss-Punkte
                        # und speichere Elemente, deren Positionen der Gauss-Punkte gewuenschte Kriterien erfuellen, zur Weiterverarbeitung ab:
                        nodeLabelConnect = np.array(list(compositeSetElements.elements[0][jj].connectivity))
                        nodeLabelConnectInst = compositeModelOdbObject.steps[step].frames[0].fieldOutputs['S'].values[0].instance
                        xCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[0] for kk in nodeLabelConnect])
                        yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
                        
#                        rCoord = ct.cart2pol_radius(xCoord,yCoord)
                        sumNodesPosRMean  = np.mean(xCoord)
                        
                        # Abspeichern der Element-Labels der Elemente, welche die gewuenschten Kriterien erfuellen:
                        if sumNodesPosRMean > rEval:
                            elements['r > r_eval'].append(stressLastFrameElementNodal[ii])
                            elementLabels['r > r_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
                        else:
                            elements['r < r_eval'].append(stressLastFrameElementNodal[ii])
                            elementLabels['r < r_eval'].append(stressLastFrameElementNodal[ii].elementLabel)
    return(elements, elementLabels)

def elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,boolInterface,colDir):
    # Bestimmung des Elementtyps:
    if colDir == 'r':
        try:
            elementTyp = elements[elementsLoc[1][0]][0].baseElementType
        except:
            try:
                elementTyp = elements[elementsLoc[1][0]][0].baseElementType
            except:
                try:
                    elementTyp = elements[elementsLoc[1][1]][0].baseElementType
                except:
                    elementTyp = elements[elementsLoc[1][1]][0].baseElementType
    elif colDir == 't':
        elementTyp = elements[elementsLoc[0][0]][0].baseElementType
    stressesElementNodes = elemNodeLocCol(colDir)
    laufVarElementNodesFieldOutput = 0
    for ii in range(len(sortedElementNodes)):
        if ii == 0:
                stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
                laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput + 1
        elif ii == len(sortedElementNodes)-1:
                stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
        else:
            if ii % 2 == 0:
                stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
                laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+2
            else:
                stressesElementNodes = stressElementNodesF3N(stressesElementNodes,step,elements,laufVarElementNodesFieldOutput,colDir)
                laufVarElementNodesFieldOutput = laufVarElementNodesFieldOutput+1
    return stressesElementNodes

def elemLocInterface(stressLastFrameGaussPoint,elementLabels,colDir):
    # Bestimmung der Position der Elemente und abspeichern relevater Daten:
    stressGaussPoints = elemGaussPointsLocCol(colDir)
    dataOutputGaussPoints = dataGaussPointsLocCol(colDir)
    if colDir == 'r':
        for ii in range(np.array(elementsLoc).shape[1]):
            for mm in range(len(elementLabels[elementsLoc[1][ii]])):
                stressGaussPointsTemp = elemGaussPointsLocCol(colDir)
                for nn in range(len(stressLastFrameGaussPoint)):
                    if elementLabels[elementsLoc[1][ii]][mm] == stressLastFrameGaussPoint[nn].elementLabel:
                        stressGaussPointsTemp[elementsLoc[1][ii]].append((stressLastFrameGaussPoint[nn], stressLastFrameGaussPoint[nn].elementLabel, stressLastFrameGaussPoint[nn].integrationPoint))
                stressGaussPointsTemp[elementsLoc[1][ii]]= sorted(stressGaussPointsTemp[elementsLoc[1][ii]],key=getIntegrationPoint)
                for nn in range(len(stressGaussPointsTemp[elementsLoc[1][ii]])):
                    stressGaussPoints[elementsLoc[1][ii]].append(stressGaussPointsTemp[elementsLoc[1][ii]][nn][0])
        # Abspeichern der Daten, die letztendlich ausgegeben werden sollen:
        for ii in range(np.array(elementsLoc).shape[1]):
            for kk in range(len(stressGaussPoints[elementsLoc[1][ii]])):
                dataOutputGaussPoints[elementsLoc[1][ii]].append((getStressesIntegrationPoint(stressGaussPoints[elementsLoc[1][ii]][kk])))
    elif colDir == 't':
        for ii in range(np.array(elementsLoc).shape[1]):
            for mm in range(len(elementLabels[elementsLoc[0][ii]])):
                stressGaussPointsTemp = elemGaussPointsLocCol(colDir)
                for nn in range(len(stressLastFrameGaussPoint)):
                    if elementLabels[elementsLoc[0][ii]][mm] == stressLastFrameGaussPoint[nn].elementLabel:
                        stressGaussPointsTemp[elementsLoc[0][ii]].append((stressLastFrameGaussPoint[nn], stressLastFrameGaussPoint[nn].elementLabel, stressLastFrameGaussPoint[nn].integrationPoint))
                stressGaussPointsTemp[elementsLoc[0][ii]] = sorted(stressGaussPointsTemp[elementsLoc[0][ii]],key=getIntegrationPoint)
                for nn in range(len(stressGaussPointsTemp[elementsLoc[0][ii]])):
                    stressGaussPoints[elementsLoc[0][ii]].append(stressGaussPointsTemp[elementsLoc[0][ii]][nn][0])
        # Abspeichern der Daten, die letztendlich ausgegeben werden sollen:
        for ii in range(np.array(elementsLoc).shape[1]):
            for kk in range(len(stressGaussPoints[elementsLoc[0][ii]])):
                dataOutputGaussPoints[elementsLoc[0][ii]].append((getStressesIntegrationPoint(stressGaussPoints[elementsLoc[0][ii]][kk])))
    return dataOutputGaussPoints

def dtypeOutputElementNodeFnc():
    return [('sigrr',float),('sigtt',float),('sigzz',float),('taurt',float),('xCoords',float),('yCoords',float),('r',float),('theta',float), ('labelNode',int),('labelElement',int)]

def dtypeOutputGaussPointFnc():
    return [('labelElement',int),('integrationPoint',int),('sigrr',float),('sigtt',float),('sigzz',float),('taurt',float)]

def stressesExtrapolatedAndIntegrationPointInterfaceCirc(set, step, iInterfaceEval, rEval):
    colDir = 't'
    compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(),step,colDir)
    
    # sortedElementNodeLabels = np.array([ii[0] for ii in sortedElementNodes])
    stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
    
    elements, elementLabels = elemDefLocCenGrav(rEval, 'None', stressLastFrameElementNodal,compositeSetElements)
    
    # Sortierung der Elemente Umfangsrichtung:
    for ii in range(np.array(elementsLoc).shape[1]):
        elements[elementsLoc[0][ii]] = sortElements(sortedElementNodes,elements[elementsLoc[0][ii]])
    
    stressesElementNodes = elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,False,colDir)
    
    # Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
    # 1: GreaterAngleEval
    # 2: SmallerAngleEval
    # Definition des Ausgabeformates des auszugebenden *.txt-Files:	
    dataOutputElemNodes = {elementNodesLocNo:{elementNodesLocR:[] for elementNodesLocR in elementNodesLoc[1]} for elementNodesLocNo in elementNodesLoc[0]}
    
    dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_1'
    
    dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_1'
    
    dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][0]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_2'
    
    dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][1]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_2'
    dtypeOutputElementNode = dtypeOutputElementNodeFnc()
    for ii in range(np.array(elementNodesLoc).shape[1]):
        for jj in range(np.array(elementNodesLoc).shape[1]):
            np.savetxt(dataOutputElemNodes[elementNodesLoc[0][ii]][elementNodesLoc[1][jj]]+'.txt', np.sort(np.array(stressesElementNodes[elementNodesLoc[0][ii]][elementNodesLoc[1][jj]], dtype=dtypeOutputElementNode),order='theta'))
    
    #---------------------------------------------------------------
    # Definition des Ausgabeformates des auszugebenden *.txt-Files:
    dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint,elementLabels,colDir)
    
    # Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
    dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
    dataOutputGaussPointsTxt[elementsLoc[0][1]]= set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_IntPo'
    dataOutputGaussPointsTxt[elementsLoc[0][0]]= set + '_Layer_' + str(plyAngle[iInterfaceEval-1]) + '_IntPo'
    dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
    for ii in range(np.array(elementsLoc).shape[1]):
        np.savetxt(dataOutputGaussPointsTxt[elementsLoc[0][ii]]+'.txt', np.array(dataOutputGaussPoints[elementsLoc[0][ii]], dtype=dtypeOutputGaussPoint))

def cart_stressesExtrapolatedAndIntegrationPointInterfaceCirc(set, step, iInterfaceEval, rEval):
    colDir = 't'
    compositeSetElements, cart_stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, cart_sortedElementNodes = cart_getStressesDisplacements(set.upper(),step,colDir)
    
    # sortedElementNodeLabels = np.array([ii[0] for ii in sortedElementNodes])
    stressLastFrameElementNodal = [cart_stressLastFrameElementNodalTemp[kk] for sortedElementNode in cart_sortedElementNodes for kk in range(len(cart_stressLastFrameElementNodalTemp)) if cart_stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
    
    elements, elementLabels = cart_elemDefLocCenGrav(rEval, 'None', stressLastFrameElementNodal,compositeSetElements)
    
    # Sortierung der Elemente Umfangsrichtung:
    for ii in range(np.array(elementsLoc).shape[1]):
        elements[elementsLoc[0][ii]] = sortElements(cart_sortedElementNodes,elements[elementsLoc[0][ii]])
    
    stressesElementNodes = elemNodesFirstLast3Nodes(elements,cart_sortedElementNodes,step,False,colDir)
    
    # Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
    # 1: GreaterAngleEval
    # 2: SmallerAngleEval
    # Definition des Ausgabeformates des auszugebenden *.txt-Files:	
    dataOutputElemNodes = {elementNodesLocNo:{elementNodesLocR:[] for elementNodesLocR in elementNodesLoc[1]} for elementNodesLocNo in elementNodesLoc[0]}
    
    dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][0]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_1'
    
    dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[1][1]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_1'
    
    dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][0]] = set + '_Layer_' + str(plyAngle[(iInterfaceEval-1)]) + '_RS_2'
    
    dataOutputElemNodes[elementNodesLoc[0][1]][elementNodesLoc[1][1]] = set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_RG_2'
    dtypeOutputElementNode = dtypeOutputElementNodeFnc()
    for ii in range(np.array(elementNodesLoc).shape[1]):
        for jj in range(np.array(elementNodesLoc).shape[1]):
            np.savetxt(dataOutputElemNodes[elementNodesLoc[0][ii]][elementNodesLoc[1][jj]]+'.txt', np.sort(np.array(stressesElementNodes[elementNodesLoc[0][ii]][elementNodesLoc[1][jj]], dtype=dtypeOutputElementNode),order='theta'))
    
    #---------------------------------------------------------------
    # Definition des Ausgabeformates des auszugebenden *.txt-Files:
    dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint,elementLabels,colDir)
    
    # Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
    dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
    dataOutputGaussPointsTxt[elementsLoc[0][1]]= set + '_Layer_' + str(plyAngle[iInterfaceEval]) + '_IntPo'
    dataOutputGaussPointsTxt[elementsLoc[0][0]]= set + '_Layer_' + str(plyAngle[iInterfaceEval-1]) + '_IntPo'
    dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
    for ii in range(np.array(elementsLoc).shape[1]):
        np.savetxt(dataOutputGaussPointsTxt[elementsLoc[0][ii]]+'.txt', np.array(dataOutputGaussPoints[elementsLoc[0][ii]], dtype=dtypeOutputGaussPoint))

#---------------------------------------------------------------

'''
Definierte Knoten-/Element-Sets werden mithilfe folgender Funktionen im Postprocessing verarbeitet, d.h. gewuenschte Daten, wie z. B. Spannungen o. Verschiebungen,  werden in einer *.txt-File abgespeichert.
Hinsichtlich der Spannungen werden diese sowohl an Elementknoten(durch Extrapolation) als auch in den Integrations-Punkten ausgegeben:
'''
def stressesExtrapolatedAndIntegrationPointRadial(set, step, tEval):
    colDir = 'r'
    compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(), step, colDir)
    
    stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel == sortedElementNode[0]]
    
    elements, elementLabels = elemDefLocCenGrav('None', tEval, stressLastFrameElementNodal, compositeSetElements)
    
    # Debugging information
    print "Elements: {}".format(elements)
    print "Element Labels: {}".format(elementLabels)
    
    # Radiale Sortierung der Elemente:
    for ii in range(np.array(elementsLoc).shape[1]):
        elements[elementsLoc[1][ii]] = sortElements(sortedElementNodes, elements[elementsLoc[1][ii]])
    
    stressesElementNodes = elemNodesFirstLast3Nodes(elements, sortedElementNodes, step, False, colDir)
    
    # Output file names
    dataOutputElemNodes = {elementNodesLoc[0][0]: {elementNodesLocT: [] for elementNodesLocT in elementNodesLoc[2]}}
    dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][1]] = set + '_TG'
    dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][0]] = set + '_TS'
    
    dtypeOutputElementNode = dtypeOutputElementNodeFnc()
    
    for jj in range(np.array(elementNodesLoc).shape[1]):
        # Check if there is data to write
        if bool(stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]]):
            output_filename = dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]] + '.txt'
            try:
                # Debug print before writing
                print "Writing element node data to: {}".format(output_filename)
                np.savetxt(output_filename, np.sort(np.array(stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]], dtype=dtypeOutputElementNode), order='r'))
            except Exception as e:
                print "Error writing to file {}: {}".format(output_filename, str(e))
    
    # Gauss points data
    dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint, elementLabels, colDir)
    
    dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
    dataOutputGaussPointsTxt[elementsLoc[1][1]] = set + '_TG_IntPo'
    dataOutputGaussPointsTxt[elementsLoc[1][0]] = set + '_TS_IntPo'
    
    dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
    
    for ii in range(np.array(elementsLoc).shape[1]):
        # Check if there is data to write
        if bool(dataOutputGaussPoints[elementsLoc[1][ii]]):
            output_filename = dataOutputGaussPointsTxt[elementsLoc[1][ii]] + '.txt'
            try:
                # Debug print before writing
                print "Writing Gauss points data to: {}".format(output_filename)
                np.savetxt(output_filename, np.array(dataOutputGaussPoints[elementsLoc[1][ii]], dtype=dtypeOutputGaussPoint))
            except Exception as e:
                print "Error writing to file {}: {}".format(output_filename, str(e))

#-------------------------------------------------------------------
# Visualization:

compositeModelOdbPath = jobName + '.odb'
compositeModelOdbObject = session.openOdb(name=compositeModelOdbPath)
session.viewports['Viewport: 1'].setValues(displayedObject=compositeModelOdbObject)
compositeModelViewport = session.viewports['Viewport: 1']

# Erstellen eines zylindrischen Koordinatensystems fuer das Postprocessing:

postProc_Cyl_CS = compositeModelOdbObject.rootAssembly.DatumCsysByThreePoints(coordSysType = CYLINDRICAL, name='Fixed cylindrical coordinate system - PostProcessing',
                                                                                point1 = (1, 0, 0),
                                                                                point2 = (0, 1, 0),
                                                                                origin= (0, 0, 0))
compositeModelViewport.odbDisplay.basicOptions.setValues(transformationType=USER_SPECIFIED, datumCsys=postProc_Cyl_CS)

# Anzeigen der Radialspannungen am verformten FE-Model:
compositeModelViewport.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S11'), )


nIncrements = len(compositeModelOdbObject.steps[step].frames)
def exportInfo(mR, mT, disp_u1, mFlanken, mBias, rBias, simTime, nIncrements):
    """
    Exportiert die Netzparameter und Simulationsinformationen in eine Datei.

    Args:
        mR (float): Anzahl der Elemente in radialer Richtung.
        mT (float): Anzahl der Elemente in tangentialer Richtung.
        disp_u1 (float): Verschiebung in Richtung u1.
        mFlanken (float): Anzahl der Elemente entlang der Flanken.
        mBias (float): Bias-Wert für das Netz.
        rBias (float): Radialer Bias-Wert.
        simTime (float): Simulationszeit in Sekunden.
        nIncrements (float): Anzahl der Inkremente.
    """
    delimiter = '\t'
    matrix = [
        [float(mR), float(mT), float(disp_u1), float(mFlanken), float(mBias), float(rBias), float(simTime), float(nIncrements)]
    ]
    export_matrix(matrix, 'meshParameters', delimiter)

exportInfo(mR, mT, disp_u1, mFlanken, mBias, rBias, simTime, nIncrements)

startTime = time.time()
#-------------------------------------------------------------------
def postProcessingCurvedCompositeCirc():
    thetaInterfaceEvalDegree = 180*angleEval/np.pi
    for ii in range(N):
        stressesExtrapolatedAndIntegrationPointRadial('FEInner'+str(ii+1),step,thetaInterfaceEvalDegree)
    for ii in range(len(iInterfaceEval)):
        stressesExtrapolatedAndIntegrationPointInterfaceCirc('FEInterfaceCircLeft'+str(iInterfaceEval[ii]),step,iInterfaceEval[ii],rInterfaceEval[ii])
    for ii in range(len(iInterfaceEval)):
        cart_stressesExtrapolatedAndIntegrationPointInterfaceCirc('FEInterfaceCircRight'+str(iInterfaceEval[ii]),step,iInterfaceEval[ii],rInterfaceEval[ii])

postProcessingCurvedCompositeCirc()

print('Berechnungszeit Simulation: ' + str(simTime) +  '  Berechnungszeit PostProcessing: ' + str(time.time()-startTime))
