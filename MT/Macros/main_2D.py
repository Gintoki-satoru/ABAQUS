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

path_modules = 'C:\\Users\\psachdeva\\Documents\\01_Promotion\\04_MT\\Kp\\09_Zusatzmaterial\\02_Script_2D'
os.chdir(path_modules)
# Further packages:
import coordinateTransformation as ct

#---------------------------------------------------------------
# Names des Modells:
modelName = '0_90_2_R4h_L4h_T1_AO90_AE45_CFK'

# Modellparameter:
# Schichtwinkel der einzelnen physikalischen Schichten:
plyAngle = [0,90,0,90]

# Definition des Pfades zur Sicherung der gesamten FE-Analyse: 
analysis_Path = 'C:\\Users\\psachdeva\\Documents\\01_Promotion\\01_ANA_2.1\\FE\\shell'

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

# Materialparameter der 0-Schicht (1 - Radialkoordinate, 2 - Tangentialkoordinate, 3 - Axialkoordinate):
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

# Startwinkel des Modells, Winkelbereich des Sub-Modells, Oeffnungswinkel, Auswertungswinkel 
angleStart,angleOpening = 0.0,np.pi/2

# Auswertungswinkel:
angleEval = 2*angleOpening/4

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
cylindricalOuterPressure, cylindricalInnerPressure = 0.0, 1.0
if cylindricalOuterPressure != 0.0 or cylindricalInnerPressure != 0.0:
	boolCylindricalPressure = True
else:
	boolCylindricalPressure = False

# Temperaturdifferenz:
tempDif = 0.0
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
		else:
			rInterfaceAll.append(rk[ii])
	# Mittlere Radien aller partitionierten physikalischen Schichten:
	rFacesAll = [(rInterfaceAll[ii]+rInterfaceAll[ii+1])/2 for ii in range(N)]
	# Vektor zur Definition der Interfaces, die eine Element-Konzentration erfahren:
	return(rm,rInterfaceEval,rInterfaceAll,rFacesAll)

def circGeometryParameters():
	# Umfangswinkel-Kanten:
	thetaEdges = [angleStart,angleEval,angleOpening]
	# Umfangswinkel-Partitionierungen:
	thetaPartition = [angleEval]
	# Umfangswinkel-Mittenkanten:
	thetaFaces = [(thetaEdges[ii]+thetaEdges[ii+1])/2 for ii in range(len(thetaEdges)-1)]
	return(thetaPartition,thetaEdges,thetaFaces)

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

rm,rInterfaceEval,rInterfaceAll,rFacesAll = radialGeometryParameters(N,rk,iInterfaceEval)
thetaPartition,thetaEdges,thetaFaces = circGeometryParameters()

#---------------------------------------------------------------
def sketchCurvedComposite():
	# Skizze/Koerper des zylindrisch gekruemmten Laminats:
	curvedCompositeSketch = compositeModel.ConstrainedSketch(name='Cylindrical_Composite_Sketch', sheetSize=2*rk[-1])
	curvedCompositeSketch.ArcByCenterEnds(center=(0.0, 0.0), point1=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0)), point2=(ct.pol2cart_x(rk[0],  thetaEdges[-1]), ct.pol2cart_y(rk[0], thetaEdges[-1])), direction=COUNTERCLOCKWISE)
	curvedCompositeSketch.ArcByCenterEnds(center=(0.0, 0.0), point1=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0)), point2=(ct.pol2cart_x(rk[-1],  thetaEdges[-1]), ct.pol2cart_y(rk[-1],  thetaEdges[-1])), direction=COUNTERCLOCKWISE)
	
	curvedCompositeSketch.Line(point1=(ct.pol2cart_x(rk[-1], 0), ct.pol2cart_y(rk[-1], 0)), point2=(ct.pol2cart_x(rk[0], 0), ct.pol2cart_y(rk[0], 0)))
	curvedCompositeSketch.Line(point1=(ct.pol2cart_x(rk[0],  thetaEdges[-1]), ct.pol2cart_y(rk[0],  thetaEdges[-1])), point2=(ct.pol2cart_x(rk[-1],  thetaEdges[-1]), ct.pol2cart_y(rk[-1],  thetaEdges[-1])))
	
	curvedCompositePart = compositeModel.Part(name='Cylindrically_Curved_Composite', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
	curvedCompositePart.BaseShell(sketch=curvedCompositeSketch)
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
	return (curvedCompositeCylCoordSys,curvedCompositeCSCylId)

curvedCompositeCylCoordSys,curvedCompositeCSCylId = coordinateSystemCurvedComposite()
#---------------------------------------------------------------
# Partitionierung in radialer Richtung:
# Einteilung des Laminats in physikalische Schichten:
def partitionCurvedCompositeRadial():
	for ii in range (N-1):
		curvedCompositeMakeSketchTransformPlane = curvedCompositePart.faces.findAt((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii],  thetaEdges[-1]/2), 0.0),)
		curvedCompositeMakeSketchTransformUpEdge = curvedCompositePart.edges.findAt((ct.pol2cart_x(rk[0],  thetaEdges[-1]/2), ct.pol2cart_y(rk[0],  thetaEdges[-1]/2), 0.0),)
		curvedCompositeMakeSketchTransform = curvedCompositePart.MakeSketchTransform(sketchPlane=curvedCompositeMakeSketchTransformPlane, sketchPlaneSide=SIDE1, sketchUpEdge=curvedCompositeMakeSketchTransformUpEdge, origin=(0, 0, 0))
		curvedCompositeTransformedSketch = compositeModel.ConstrainedSketch(name='Transformed_Sketch', sheetSize=2*rk[-1], transform=curvedCompositeMakeSketchTransform)
		
		curvedCompositePart.projectReferencesOntoSketch(sketch=curvedCompositeTransformedSketch, filter=COPLANAR_EDGES)
		curvedCompositeTransformedSketch.ArcByCenterEnds(center=(0.0,0.0), point1=(ct.pol2cart_x(rk[ii+1], 0), ct.pol2cart_y(rk[ii+1], 0)),point2=(ct.pol2cart_x(rk[ii+1],  thetaEdges[-1]), ct.pol2cart_y(rk[ii+1],  thetaEdges[-1])),direction=COUNTERCLOCKWISE)
		
		curvedCompositePickedFaces = curvedCompositePart.faces.findAt(((ct.pol2cart_x(rm[ii],  thetaEdges[-1]/2), ct.pol2cart_y(rm[ii], thetaEdges[-1]/2), 0.0),),)
		curvedCompositePart.PartitionFaceBySketch(faces=curvedCompositePickedFaces, sketch=curvedCompositeTransformedSketch, sketchUpEdge=curvedCompositeMakeSketchTransformUpEdge)

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
		
		compositePartitionFace = curvedCompositePart.faces
		curvedCompositePart.DatumPlaneByThreePoints(point1=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar]], 
																		point2=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar+1]],
																		point3=curvedCompositePart.datums[compositePartDatumsKeys[compositePartitionDatumPointLaufVar+2]])
		compositePartitionPlaneId = curvedCompositePart.features['Datum plane-'+ str(compositePartitionPlaneLaufVar)].id
		curvedCompositePart.PartitionFaceByDatumPlane(faces=compositePartitionFace, datumPlane=curvedCompositePart.datums[compositePartitionPlaneId])
		
		compositePartitionPlaneLaufVar += 1

partitionCurvedCompositeRadial()
partitionCurvedCompositeCirc()
#---------------------------------------------------------------
def materialMetal():
	# Materialdefinition ueber Ingenieurkonstanten:
	aluminumMaterial = compositeModel.Material(metalMaterialName) 
	aluminumMaterial.Elastic(table=((E2024, Nu2024), ))
	compositeModel.HomogeneousSolidSection(name='Metal_Section', material=metalMaterialName, thickness=None)

def materialComposite():
	compositeMaterial = compositeModel.Material(compositeMaterialName+'_0')
	compositeMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=((E1, E2, E3, Nu12, Nu13, Nu23, G23, G13, G12), ))
	compositeMaterial.Expansion(type=ORTHOTROPIC, table=((alpha11, alpha22, alpha33), ))
	compositeModel.HomogeneousSolidSection(name='Composite_Section'+'_0', material=compositeMaterialName+'_0', thickness=None)
	
	compositeMaterial = compositeModel.Material(compositeMaterialName+'_90')
	compositeMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=((E1, E3, E2, Nu13, Nu12, Nu23*E3/E2, G23, G12, G13), ))
	compositeMaterial.Expansion(type=ORTHOTROPIC, table=((alpha11, alpha33, alpha22), ))
	compositeModel.HomogeneousSolidSection(name='Composite_Section'+'_90', material=compositeMaterialName+'_90', thickness=None)

materialComposite()
if boolHybrid:
	materialMetal()

def sectionCurvedComposite():
	# Section-Zuweisung bzw. Definition der Faserorientierungswinkel der einzelnen physikalischen Schichten:
	jj = 0
	for ii in range(N):
		curvedCompositeAngleInc = 0
		for mm in range(len(thetaFaces)):
			curvedCompositeFaces = curvedCompositePart.faces.findAt(((ct.pol2cart_x(rFacesAll[jj], thetaFaces[mm]), ct.pol2cart_y(rFacesAll[jj], thetaFaces[mm]), 0.0),))
			
			curvedCompositeRegion = regionToolset.Region(faces=curvedCompositeFaces)
			if plyAngle[ii] == 2024:
				curvedCompositePart.SectionAssignment(region=curvedCompositeRegion, sectionName='Metal_Section', offsetType=MIDDLE_SURFACE)
			elif plyAngle[ii] == 0:
				curvedCompositePart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_'+str(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
			elif plyAngle[ii] == 90:
				curvedCompositePart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_'+str(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
			else:
				raise ValueError('Please check the laminate layup. Only cross-plies allowed.')
			curvedCompositePart.MaterialOrientation(region=curvedCompositeRegion, orientationType=SYSTEM, localCsys=curvedCompositePart.datums[curvedCompositeCSCylId], axis=AXIS_3, additionalRotationType=ROTATION_NONE, angle=0.0, stackDirection=STACK_3)
			curvedCompositeAngleInc = curvedCompositeAngleInc + 1
		jj = jj + 1

sectionCurvedComposite()
#---------------------------------------------------------------
# Vernetzung in radialer Richtung:
mR,mT = 5,100
mRRatio,mTRatio = 1,1
def meshCurvedCompositeGlobal():
	for ii in range(len(rFacesAll)):
		for jj in range(len(thetaEdges)):
			curvedCompositeMeshEdge = curvedCompositePart.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj]), ct.pol2cart_y(rFacesAll[ii], thetaEdges[jj]),0.0),))
			curvedCompositePart.seedEdgeByBias(biasMethod=DOUBLE, endEdges=curvedCompositeMeshEdge, ratio=mRRatio, number=mR, constraint=FINER)
	# Vernetzung in tangentialer Richtung:
	for ii in range (len(rInterfaceAll)):
		for jj in range(len(thetaFaces)):
			curvedCompositeMeshEdge = curvedCompositePart.edges.findAt(((ct.pol2cart_x(rInterfaceAll[ii], thetaFaces[jj]), ct.pol2cart_y(rInterfaceAll[ii], thetaFaces[jj]), 0.0),))
			if int(round(mT*(thetaEdges[jj+1]-thetaEdges[jj])/thetaEdges[-1])) < 1:
				curvedCompositePart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mTRatio, number=1, constraint=FINER)
			else:
				curvedCompositePart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=curvedCompositeMeshEdge, ratio=mTRatio, number=int(round(mT*(thetaEdges[jj+1]-thetaEdges[jj])/thetaEdges[-1])), constraint=FINER)
	# Element types: C3D8, C3D8R, C3D20, C3D20R
	compositeElementType = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
	curvedCompositePart.setMeshControls(regions=curvedCompositePart.faces, elemShape=HEX, technique=STRUCTURED)
	curvedCompositePart.setElementType(regions=(curvedCompositePart.faces[:], ), elemTypes=(compositeElementType, ))
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
	curvedCompositeFixRadialEdge = curvedCompositeInstance.vertices.findAt(((ct.pol2cart_x(r0,  thetaEdges[-1]), ct.pol2cart_y(r0, thetaEdges[-1]), 0.0),))
	curvedCompositeFixRadialSet = compositeAssembly.Set(vertices=curvedCompositeFixRadialEdge, name='Set_Fix_U_Interface_'+ str(0))
	compositeModel.DisplacementBC(name='Fix_U_Interface_'+ str(0), createStepName=stepPrevious, region=curvedCompositeFixRadialSet, localCsys=curvedCompositeCylCoordSys, u1=0.0)
	# Sperren der tangentialen Bewegung
	for ii in range(len(rFacesAll)):
		curvedCompositeStringFixCircFaces = ''
		curvedCompositeStringFixCircFaces = curvedCompositeStringFixCircFaces + str(((ct.pol2cart_x(rFacesAll[ii], thetaEdges[-1]), ct.pol2cart_y(rFacesAll[ii], thetaEdges[-1]), 0.0),),) + ','
		curvedCompositeFixCircFacesExec = 'curvedCompositeFixCircFaces = curvedCompositeInstance.edges.findAt(' + curvedCompositeStringFixCircFaces + ')'
		exec(curvedCompositeFixCircFacesExec)
		
		curvedCompositeFixCircSet = compositeAssembly.Set(edges=curvedCompositeFixCircFaces, name='Set_Fix_V_Layer_'+ str(ii+1))
		compositeModel.DisplacementBC(name='Fix_V_Layer_'+ str(ii+1), createStepName=stepPrevious, region=curvedCompositeFixCircSet, localCsys=curvedCompositeCylCoordSys, u2 = 0.0)

def boundaryConditionSurfaceLoad():
	for jj,ll in enumerate([thetaEdges[0],thetaEdges[-1]]):
		curvedCompositeStringFixCircFaces = ''
		for ii in range(len(rFacesAll)):
			curvedCompositeStringFixCircFaces = curvedCompositeStringFixCircFaces + str(((ct.pol2cart_x(rFacesAll[ii], ll), ct.pol2cart_y(rFacesAll[ii], ll), 0.0),),) + ','
			curvedCompositeFixCircFacesExec = 'curvedCompositeFixCircFaces = curvedCompositeInstance.edges.findAt(' + curvedCompositeStringFixCircFaces + ')'
			exec(curvedCompositeFixCircFacesExec)
		curvedCompositeFixCircSet = compositeAssembly.Set(edges=curvedCompositeFixCircFaces, name='Set_Fix_U_SS_'+ str(jj))
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
	compositeAssembly.ReferencePoint(point=(ct.pol2cart_x((rk[0]+rk[-1])/2, jj*angleOpening), ct.pol2cart_y((rk[0]+rk[-1])/2, jj*angleOpening), 0.0))
	curvedCompositeLoadRefPoint = (compositeAssembly.referencePoints.findAt((ct.pol2cart_x((rk[0]+rk[-1])/2, jj*angleOpening), ct.pol2cart_y((rk[0]+rk[-1])/2, jj*angleOpening), 0.0),),)
	curvedCompositeLoadRefPointRegion = regionToolset.Region(referencePoints=curvedCompositeLoadRefPoint)
	curvedCompositeFixRefPointSet = compositeAssembly.Set(region=curvedCompositeLoadRefPointRegion, name='Set_Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)))
	compositeModel.DisplacementBC(name='Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeFixRefPointSet, localCsys=curvedCompositeCylCoordSys, u3=0.0, ur1=0.0, ur2=0.0)
	# Kinematic Coupling:
	curvedCompositeStringKinCoupFaces =''
	for ii in range(len(rFacesAll)):
		curvedCompositeStringKinCoupFaces = curvedCompositeStringKinCoupFaces + str(((ct.pol2cart_x(rFacesAll[ii], jj*angleOpening), ct.pol2cart_y(rFacesAll[ii], jj*angleOpening), 0.0),),) + ','
	curvedCompositeKinCoupFacesExec = 'curvedCompositeKinCoupFaces = curvedCompositeInstance.edges.findAt(' + curvedCompositeStringKinCoupFaces + ')'
	exec(curvedCompositeKinCoupFacesExec)
	curvedCompositeKinCoupRegion = regionToolset.Region(edges=curvedCompositeKinCoupFaces)
	if boolBendingMoment:
		compositeModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion,influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=OFF, u2=ON, ur3=ON)
	else:
		compositeModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion,influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=ON, u2=ON, ur3=ON)
	compositeModel.constraints['Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi))].setValues(localCsys=curvedCompositeCylCoordSys)
	return curvedCompositeLoadRefPointRegion

def loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion):
	# Aeussere Last - Moment:
	compositeModel.Moment(name='Load_BendingMoment_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeLoadRefPointRegion, cm3=((-1)**(jj))*bendingMoment, localCsys=curvedCompositeCylCoordSys, distributionType=UNIFORM)

def loadsRadialForce():
	# Aeussere Last - Querkraft:
	compositeModel.ConcentratedForce(name='Load_RadialForce', createStepName=step, region=curvedCompositeLoadRefPointRegion, cf1=radialForce, localCsys=curvedCompositeCylCoordSys, distributionType=UNIFORM)

def loadsCircForce():
	# Aeussere Last - Normalkraft:
	compositeModel.ConcentratedForce(name='Load_CircForce', createStepName=step, region=curvedCompositeLoadRefPointRegion, cf2=-circForce, localCsys=curvedCompositeCylCoordSys, distributionType=UNIFORM)

def loadsSurfaceLoad():
	if OuterPressure != 0.0 or cylindricalOuterPressure != 0.0:
		# Aeussere Last - Aussendruck:
		curvedCompositeStringFaces = ''
		for jj in range(len(thetaFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[-1], thetaFaces[jj]), ct.pol2cart_y(rk[-1], thetaFaces[jj]), 0.0),),) + ','
		
		curvedCompositeStringFacesExec = 'curvedCompositeOuterPressureFaces = curvedCompositeInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		
		if boolCylindricalPressure:
			compositeModel.ExpressionField(name='Load_OuterPressure_CylindricalBending', localCsys=curvedCompositeCylCoordSys, description='',expression='sin(pi*Th/angleOpening)')
			compositeModel.Pressure(name='Load_OuterPressure_CylindricalBending', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeOuterPressureFaces), magnitude=cylindricalOuterPressure, distributionType=FIELD, field='Load_OuterPressure_CylindricalBending', amplitude=UNSET)
		else:
			compositeModel.Pressure(name='Load_OuterPressure', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeOuterPressureFaces), magnitude=OuterPressure, distributionType=UNIFORM)
	if InnerPressure != 0.0 or cylindricalInnerPressure != 0.0:
		# Aeussere Last - Innendruck:
		curvedCompositeStringFaces = ''
		for jj in range(len(thetaFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[0], thetaFaces[jj]), ct.pol2cart_y(rk[0], thetaFaces[jj]), 0.0),),) + ','
		
		curvedCompositeStringFacesExec = 'curvedCompositeInnerPressureFaces = curvedCompositeInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		
		if boolCylindricalPressure:
			compositeModel.ExpressionField(name='Load_InnerPressure_CylindricalBending', localCsys=curvedCompositeCylCoordSys, description='',expression='sin(pi*Th/angleOpening)')
			compositeModel.Pressure(name='Load_InnerPressure_CylindricalBending', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeInnerPressureFaces), magnitude=cylindricalInnerPressure, distributionType=FIELD, field='Load_InnerPressure_CylindricalBending', amplitude=UNSET)
		else:
			compositeModel.Pressure(name='Load_InnerPressure', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeInnerPressureFaces), magnitude=InnerPressure, distributionType=UNIFORM)

def loadsHygrothermal():
	# Definition der Ausgangstemperatur:
	compositeModel.Temperature(name='Temp_Initial', createStepName=stepPrevious, region=regionToolset.Region(faces=curvedCompositeInstance.faces), distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0,))
	compositeModel.Temperature(name='Temperature_Difference', createStepName=step, region=regionToolset.Region(faces=curvedCompositeInstance.faces), distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(tempDif,))

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
	loadsSurfaceLoad()
elif any([boolTempDif,boolMoistDif]):
	loadsHygrothermal()
else:
	pass

#---------------------------------------------------------------
def setsCurvedCompositeCirc():
	# Sets zur Verifizierung der Inner-Solution:
	# Sets der Spannungen in radialer Richtung in Laminatmitte:
	for ii in range(len(rFacesAll)):
		compositeSetEvalEdge = curvedCompositeInstance.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], angleEval), ct.pol2cart_y(rFacesAll[ii], angleEval), 0.0),))
		
		compositeAssembly.Set(edges=compositeSetEvalEdge, name='FEInner'+str(ii+1))
	# Sets zur Bestimmung des Abklingverhaltens der singulaeren Spannungen (in Dickenrichtung) in den Auswertungsinterfaces an beiden Enden:
	for kk in range(len(iInterfaceEval)):
		curvedCompositeStringFaces = ''
		ii = iInterfaceEval[kk]
		for jj in range(len(thetaFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[ii], thetaFaces[jj]), ct.pol2cart_y(rk[ii], thetaFaces[jj]), 0.0),),) + ','
		curvedCompositeStringFacesExec = 'compositeSetEvalEdge = curvedCompositeInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		compositeAssembly.Set(edges=compositeSetEvalEdge, name='FEInterfaceCirc'+str(iInterfaceEval[kk]))

setsCurvedCompositeCirc()
# ---------------------------------------------------------------
# Field Output:
compositeModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'E', 'COORD'))

# ---------------------------------------------------------------
# History Output

# ---------------------------------------------------------------
# Speichern der CAE-Datei:
mdb.saveAs(pathName=analysis_newPath + '\\' + modelName + '_CAE_Model')

#-------------------------------------------------------------------
# Job
# Erstellen des Jobs:
jobName = modelName + '_Job'
mdb.Job(name=jobName, model=modelName,
		description='Run FE-analysis', parallelizationMethodExplicit=DOMAIN, numDomains=1, numCpus=1, memory=80,echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)

# Job starten:
mdb.jobs[jobName].submit(consistencyChecking=OFF)
mdb.jobs[jobName].waitForCompletion()

elementsLoc = [['r < r_eval','r > r_eval'],['t < t_eval','t > t_eval']]
elementNodesLoc = [['first three element nodes','last three element nodes'],['r < r_eval','r > r_eval'],['t < t_eval','t > t_eval']]

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

def getCoordinateT(unsortedElementNodes):
	return ct.cart2pol_theta(unsortedElementNodes[1],unsortedElementNodes[2])

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

#---------------------------------------------------------------
'''
Definierte Knoten-/Element-Sets werden mithilfe folgender Funktionen im Postprocessing verarbeitet, d.h. gewuenschte Daten, wie z. B. Spannungen o. Verschiebungen,  werden in einer *.txt-File abgespeichert.
Hinsichtlich der Spannungen werden diese sowohl an Elementknoten(durch Extrapolation) als auch in den Integrations-Punkten ausgegeben:
'''
def stressesExtrapolatedAndIntegrationPointRadial(set, step, tEval):
	colDir = 'r'
	compositeSetElements, stressLastFrameElementNodalTemp, stressLastFrameGaussPoint, sortedElementNodes = getStressesDisplacements(set.upper(),step,colDir)
	
	stressLastFrameElementNodal = [stressLastFrameElementNodalTemp[kk] for sortedElementNode in sortedElementNodes for kk in range(len(stressLastFrameElementNodalTemp)) if stressLastFrameElementNodalTemp[kk].nodeLabel  == sortedElementNode[0]]
	
	elements, elementLabels = elemDefLocCenGrav('None', tEval, stressLastFrameElementNodal,compositeSetElements)
	
	# Radiale Sortierung der Elemente:
	for ii in range(np.array(elementsLoc).shape[1]):
		elements[elementsLoc[1][ii]] = sortElements(sortedElementNodes,elements[elementsLoc[1][ii]])
	
	stressesElementNodes = elemNodesFirstLast3Nodes(elements,sortedElementNodes,step,False,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	# 1: GreaterAngleEval
	# 2: SmallerAngleEval
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:	
	dataOutputElemNodes = {elementNodesLoc[0][0]:{elementNodesLocT:[] for elementNodesLocT in elementNodesLoc[2]}}
	
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][1]] = set + '_TG'
	dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][0]] = set + '_TS'
	
	dtypeOutputElementNode = dtypeOutputElementNodeFnc()
	for jj in range(np.array(elementNodesLoc).shape[1]):
		if bool(stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]]):
			np.savetxt(dataOutputElemNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]]+'.txt', np.sort(np.array(stressesElementNodes[elementNodesLoc[0][0]][elementNodesLoc[2][jj]], dtype=dtypeOutputElementNode),order='r'))
		
	#---------------------------------------------------------------
	# Definition des Ausgabeformates des auszugebenden *.txt-Files:
	dataOutputGaussPoints = elemLocInterface(stressLastFrameGaussPoint,elementLabels,colDir)
	
	# Ausgabe der zuvor definierten Daten in einem *.txt-File, geordnet nach der Dickenkoordinate des zylindrisch gekruemmten Laminats und speichern dieser im Arbeitsverzeichnis:
	dataOutputGaussPointsTxt = dataGaussPointsLocCol(colDir)
	dataOutputGaussPointsTxt[elementsLoc[1][1]] = set + '_TG_IntPo'
	dataOutputGaussPointsTxt[elementsLoc[1][0]] = set + '_TS_IntPo'
	dtypeOutputGaussPoint = dtypeOutputGaussPointFnc()
	for ii in range(np.array(elementsLoc).shape[1]):
		if bool(dataOutputGaussPoints[elementsLoc[1][ii]]):
			np.savetxt(dataOutputGaussPointsTxt[elementsLoc[1][ii]]+'.txt', np.array(dataOutputGaussPoints[elementsLoc[1][ii]], dtype=dtypeOutputGaussPoint))

#-------------------------------------------------------------------
# Visualization:
compositeModelOdbPath = jobName + '.odb'
compositeModelOdbObject = session.openOdb(name=compositeModelOdbPath)
session.viewports['Viewport: 1'].setValues(displayedObject=compositeModelOdbObject)
compositeModelViewport = session.viewports['Viewport: 1']

# Erstellen eines zylindrischen Koordinatensystems fuer das Postprocessing:
postProc_Cyl_CS = compositeModelOdbObject.rootAssembly.DatumCsysByThreePoints(coordSysType = CYLINDRICAL, name='Fixed cylindrical coordinate system - PostProcessing',
																				point1 = (rk[-1], 0, length/2),
																				point2 = (0, rk[-1], length/2),
																				origin= (0, 0, length/2))
compositeModelViewport.odbDisplay.basicOptions.setValues(transformationType=USER_SPECIFIED, datumCsys=postProc_Cyl_CS)

# Anzeigen der Radialspannungen am verformten FE-Model:
compositeModelViewport.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S11'), )

startTime = time.time()
#-------------------------------------------------------------------
# Postprocessing:
def postProcessingCurvedCompositeCirc():
	thetaInterfaceEvalDegree = 180*angleEval/np.pi
	# Auswertung der im Preprocessing definierten Knoten-/Element-Sets mithilfe der definierten Funktionen:
	
	for ii in range(N):
		stressesExtrapolatedAndIntegrationPointRadial('FEInner'+str(ii+1),step,thetaInterfaceEvalDegree)

postProcessingCurvedCompositeCirc()

print(time.time()-startTime)
