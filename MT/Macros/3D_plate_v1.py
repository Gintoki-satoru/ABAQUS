''' File name: Script_FlatPlateAnalysis.py
    Date created: 2025-07-09
'''

#---------------------------------------------------------------
# Python-Packages:
import __main__
import os
import datetime as dt
import numpy as np

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

path_modules = 'C:\\Users\\psachdeva\\Documents\\01_Promotion\\01_ANA_2.1\\FE'
os.chdir(path_modules)

#---------------------------------------------------------------
# Model Parameters:
modelName = 'Flat_Plate_Model'
plyAngle = [0, 90, 90, 0]
analysis_Path = 'C:\\Users\\psachdeva\\Documents\\01_Promotion\\01_ANA_2.1\\FE\\flat_plate'
analysis_currentDateTime = dt.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
analysis_newPath = analysis_Path + '\\' + modelName + '_' + analysis_currentDateTime
os.makedirs(analysis_newPath)
os.chdir(analysis_newPath)
N = len(plyAngle)
dL = 0.125
mX, mY,mZ = 2, 500, 500
plateLength = 100.0
plateWidth = 100.0
plateThickness = N * dL
bendingMoment = plateLength

compositeMaterialName = 'cfk'
if compositeMaterialName == 'cfk':
    E1, E2, E3 = 7460.0, 118148.0, 7460.0
    Nu12, Nu13, Nu23 = 0.021, 0.37, 0.34
    G12, G13, G23 = 4800.0, 2701.0, 4800.0
    alpha11, alpha22, alpha33 = 2.6e-5, -1.0e-6, 2.6e-5
elif compositeMaterialName == 'gfk':
    E1, E2, E3 = 9552.6, 39296.0, 9552.6
    Nu21, Nu13, Nu23 = 0.29, 0.38, 0.29
    Nu12 = E1 / E2 * Nu21
    G12, G13, G23 = 3080.5, 3449.0, 3080.5
    alpha11, alpha22, alpha33 = 2.6e-5, 8.6e-6, 2.6e-5
elif compositeMaterialName == 'cfkDuro':
    E1, E2, E3 = 6895.0, 172375.0, 6895.0
    Nu21, Nu13, Nu23 = 0.25, 0.25, 0.25
    Nu12 = E1 / E2 * Nu21
    G12, G13, G23 = 3448.0, 1379.0, 3448.0
    alpha11, alpha22, alpha33 = 2.6e-5, -1.0e-6, 2.6e-5

#---------------------------------------------------------------
# Model Definition
def export_matrix(matrix, filename, delimiter='\t'):
    with open(filename, 'w') as file:
        for row in matrix:
            formatted_row = delimiter.join("{:.0f}".format(x) if isinstance(x, (int, float)) and x == int(x) else "{:.2f}".format(x) for x in row)
            file.write(formatted_row + '\n')

def export_composite_layup(plyAngle):
    N = len(plyAngle)
    pAExport = np.zeros((1, N))
    for ii in range(N):
        pAExport[0, ii] = plyAngle[ii]
    export_matrix(pAExport, 'compositeLayup', delimiter='\t')

export_composite_layup(plyAngle)

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
    # Expor der Matrix mit benutzerdefiniertem Format
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
dL = dL # Example value
width_h = plateLength  # Example value
angleBegin = 0  # Example begin angle in radians
openingAngle = plateWidth  # Example opening angle in radians
angleEval = plateWidth/2
export_geometry(plyAngle, MDetailed, orderLagrangeDetailed, Ratio, dL, width_h, angleBegin, openingAngle, angleEval)

Mdb()
session.viewports['Viewport: 1'].setValues(displayedObject=None)
mdb.models.changeKey(fromName='Model-1', toName=modelName)
compositeModel = mdb.models[modelName]

def sketchFlatPlate():
    flatPlatePart = compositeModel.Part(name='Flat_Plate', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    sketchTransform = flatPlatePart.MakeSketchTransform(
        sketchPlane=flatPlatePart.datums[flatPlatePart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0).id],
        sketchPlaneSide=SIDE1,
        origin=(0.0, 0.0, 0.0)
    )
    flatPlateSketch = compositeModel.ConstrainedSketch(name='Flat_Plate_Sketch', sheetSize=10.0, transform=sketchTransform)
    flatPlateSketch.rectangle(point1=(0.0, 0.0), point2=(plateLength, plateWidth))
    flatPlatePart.BaseSolidExtrude(sketch=flatPlateSketch, depth=plateThickness)
    return flatPlatePart

flatPlatePart = sketchFlatPlate()

compositeAssembly = compositeModel.rootAssembly
def instanceFlatPlate():
    flatPlateInstance = compositeAssembly.Instance(name='Flat_Plate_Instance', part=flatPlatePart, dependent=ON)
    return flatPlateInstance

flatPlateInstance = instanceFlatPlate()

def partitionFlatPlateThickness():
    for i in range(1, N):
        partitionPlane = i * dL
        datumPlaneFeature = flatPlatePart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=partitionPlane)
        datumPlane = flatPlatePart.datums[datumPlaneFeature.id]
        flatPlatePart.PartitionCellByDatumPlane(cells=flatPlatePart.cells, datumPlane=datumPlane)

partitionFlatPlateThickness()

def materialComposite():
    compositeMaterial = compositeModel.Material(compositeMaterialName)
    compositeMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=((E1, E2, E3, Nu12, Nu13, Nu23, G23, G13, G12), ))
    compositeMaterial.Expansion(type=ORTHOTROPIC, table=((alpha11, alpha22, alpha33), ))
    compositeModel.HomogeneousSolidSection(name='Composite_Section', material=compositeMaterialName, thickness=None)

materialComposite()

def sectionFlatPlate():
    for i in range(N):
        x_center = (i + 0.5) * dL
        ply_cells = flatPlatePart.cells.findAt(((x_center, plateWidth / 2, -plateLength / 2),))
        region = regionToolset.Region(cells=ply_cells)
        flatPlatePart.SectionAssignment(region=region, sectionName='Composite_Section', offsetType=MIDDLE_SURFACE)
        flatPlatePart.MaterialOrientation(
            region=region,
            orientationType=SYSTEM,
            axis=AXIS_1,
            additionalRotationType=ROTATION_ANGLE,
            angle=plyAngle[i],
            stackDirection=STACK_1
        )

sectionFlatPlate()

def partitionFlatPlateInterfaces():
    datumPlaneFeatureXZ = flatPlatePart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=plateWidth / 2)
    datumPlaneXZ = flatPlatePart.datums[datumPlaneFeatureXZ.id]
    flatPlatePart.PartitionCellByDatumPlane(cells=flatPlatePart.cells, datumPlane=datumPlaneXZ)
    datumPlaneFeatureXY = flatPlatePart.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=-plateLength / 2)
    datumPlaneXY = flatPlatePart.datums[datumPlaneFeatureXY.id]
    flatPlatePart.PartitionCellByDatumPlane(cells=flatPlatePart.cells, datumPlane=datumPlaneXY)

partitionFlatPlateInterfaces()

def meshFlatPlate(mX, mY, mZ):
    # Thickness direction (x)
    for i in range(N):
        midPlane = (i + 0.5) * dL
        edge = flatPlatePart.edges.findAt(((midPlane, plateWidth / 2, -plateLength / 2),))
        flatPlatePart.seedEdgeByNumber(edges=edge, number=mX, constraint=FINER)
    # Width direction (y): seed edges at constant x and z
    for i in range(N+1):
        xx = i * dL
        for zz in [-plateLength, -plateLength/2, 0]:
            # Lower y edge (near y=0)
            edge1 = flatPlatePart.edges.findAt(((xx, 0.0001, zz),))
            flatPlatePart.seedEdgeByNumber(edges=edge1, number=mY, constraint=FINER)
            # Upper y edge (near y=plateWidth)
            edge2 = flatPlatePart.edges.findAt(((xx, plateWidth - 0.0001, zz),))
            flatPlatePart.seedEdgeByNumber(edges=edge2, number=mY, constraint=FINER)
    # Length direction (z): seed edges at constant x and y
    for i in range(N+1):
        xx = i * dL
        for yy in [0, plateWidth/2,plateWidth]:
            edge1 = flatPlatePart.edges.findAt(((xx, yy, -0.0001),))
            flatPlatePart.seedEdgeByNumber(edges=edge1, number=mZ, constraint=FINER)
            edge2 = flatPlatePart.edges.findAt(((xx, yy, -plateLength+0.0001),))
            flatPlatePart.seedEdgeByNumber(edges=edge2, number=mZ, constraint=FINER)
    flatPlatePart.setMeshControls(regions=flatPlatePart.cells, elemShape=HEX, technique=STRUCTURED)
    flatPlatePart.setElementType(regions=(flatPlatePart.cells,), elemTypes=(mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD),))
    flatPlatePart.generateMesh()

meshFlatPlate(mX, mY, mZ)

def defineRegionsOfInterest():
    for ii in range(N):
        compositeSetEvalEdge = flatPlateInstance.edges.findAt(
            (((ii+0.5) * dL, plateWidth / 2, -plateLength / 2),)
        )
        compositeAssembly.Set(edges=compositeSetEvalEdge, name='FEInner' + str(ii + 1))
    for kk in range(N - 1):
        compositeSetEvalEdges = []
        for yy in np.linspace(0.01, plateWidth-0.01, num=10):
            compositeSetEvalEdge = flatPlateInstance.edges.findAt(
                (((kk + 1) * dL, yy, -plateLength / 2),)
            )
            compositeSetEvalEdges.append(compositeSetEvalEdge)
        compositeAssembly.Set(edges=tuple(compositeSetEvalEdges), name='FEInterfaceCircRight' + str(kk + 1))

defineRegionsOfInterest()

def stepMechanicalLoads():
    stepPrevious, step, stepDescription = 'Initial', 'Mechanical_Loading', 'Introduce mechanical loadings'
    compositeModel.StaticStep(name=step, previous=stepPrevious, description=stepDescription, nlgeom=OFF, initialInc=1.0)
    return stepPrevious, step

stepPrevious, step = stepMechanicalLoads()

def boundaryConditionAndLoads():
    refPoint1Feature = compositeAssembly.ReferencePoint(point=(plateThickness / 2, 0.0, -plateLength / 2))
    refPoint2Feature = compositeAssembly.ReferencePoint(point=(plateThickness / 2, plateWidth, -plateLength / 2))
    refPoint1 = compositeAssembly.referencePoints[refPoint1Feature.id]
    refPoint2 = compositeAssembly.referencePoints[refPoint2Feature.id]
    refRegion1 = compositeAssembly.Set(referencePoints=(refPoint1,), name='Ref_Point_1')
    refRegion2 = compositeAssembly.Set(referencePoints=(refPoint2,), name='Ref_Point_2')
    compositeModel.DisplacementBC(
        name='Fix_Ref_Point_1',
        createStepName=step,
        region=refRegion1,
        u3=0.0, ur1=0.0, ur2=0.0
    )
    compositeModel.DisplacementBC(
        name='Fix_Ref_Point_2',
        createStepName=step,
        region=refRegion2,
        u3=0.0, ur1=0.0, ur2=0.0
    )
    surfaceRegion1Faces = []
    surfaceRegion2Faces = []
    for ii in range(N):
        for zz in np.linspace(-0.00001, -plateLength+0.00001, num=10):
            surfaceRegion1Faces.append(flatPlateInstance.faces.findAt((((ii+0.5) * dL, 0.0, zz),)))
            surfaceRegion2Faces.append(flatPlateInstance.faces.findAt((((ii+0.5) * dL, plateWidth, zz),)))
    surfaceRegion1 = compositeAssembly.Set(faces=tuple(surfaceRegion1Faces), name='Surface_Region_1')
    surfaceRegion2 = compositeAssembly.Set(faces=tuple(surfaceRegion2Faces), name='Surface_Region_2')
    compositeModel.Coupling(
        name='Kinematic_Coupling_1',
        controlPoint=refRegion1,
        surface=surfaceRegion1,
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        u1=OFF, u2=ON, u3=OFF,
        ur1=OFF, ur2=OFF, ur3=ON
    )
    compositeModel.Coupling(
        name='Kinematic_Coupling_2',
        controlPoint=refRegion2,
        surface=surfaceRegion2,
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        u1=OFF, u2=ON, u3=OFF,
        ur1=OFF, ur2=OFF, ur3=ON
    )
    compositeModel.Moment(
        name='Bending_Moment_1',
        createStepName=step,
        region=refRegion1,
        cm3=bendingMoment,
        distributionType=UNIFORM
    )
    compositeModel.Moment(
        name='Bending_Moment_2',
        createStepName=step,
        region=refRegion2,
        cm3=-bendingMoment,
        distributionType=UNIFORM
    )

boundaryConditionAndLoads()

compositeModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'COORD'))

jobName = modelName + '_Job'
mdb.Job(
    name=jobName,
    model=modelName,
    description='Run FE-analysis',
    parallelizationMethodExplicit=DOMAIN,
    numDomains=12,
    numCpus=6,
    memory=85,
    echoPrint=OFF,
    modelPrint=OFF,
    contactPrint=OFF,
    historyPrint=OFF
)
mdb.jobs[jobName].submit(consistencyChecking=OFF)
mdb.jobs[jobName].waitForCompletion()

odbPath = jobName + '.odb'

#---------------------------------------------------------------
# Post-processing

def getStressComponents():
    return ['Sigma_xx', 'Sigma_yy', 'Sigma_zz', 'Tau_xy', 'Tau_xz', 'Tau_yz']

def saveNodalStress(data, filename):
    with open(filename, 'w') as file:
        for nodeLabel, elemLabel, x, y, z, stress in data:
            file.write("{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:d}\t{:d}\n".format(
                stress['Sigma_xx'], stress['Sigma_yy'], stress['Sigma_zz'],
                stress['Tau_xy'], stress['Tau_xz'], stress['Tau_yz'],
                x, y, z, nodeLabel, elemLabel
            ))

def cartesian_elemDefLocCenGrav_TS_TG_ZS_ZG(yEval, zEval, stressLastFrameElementNodal, compositeSetElements):
    """
    Assign elements to TS_ZS, TS_ZG, TG_ZS, TG_ZG by centroid y and z (with tolerance).
    Returns: elements (dict), elementLabels (dict)
    """
    elements = {'TS_ZS': [], 'TS_ZG': [], 'TG_ZS': [], 'TG_ZG': []}
    elementLabels = {'TS_ZS': [], 'TS_ZG': [], 'TG_ZS': [], 'TG_ZG': []}
    elementsExist = []
    for ii in range(len(stressLastFrameElementNodal)):
        if stressLastFrameElementNodal[ii].elementLabel in elementsExist:
            continue  # already processed
        elementsExist.append(stressLastFrameElementNodal[ii].elementLabel)
        for jj in range(len(compositeSetElements.elements[0])):
            if compositeSetElements.elements[0][jj].label == stressLastFrameElementNodal[ii].elementLabel:
                nodeLabelConnect = np.array(list(compositeSetElements.elements[0][jj].connectivity))
                nodeLabelConnectInst = stressLastFrameElementNodal[ii].instance
                yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
                zCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[2] for kk in nodeLabelConnect])
                yCentroid = np.mean(yCoord)
                zCentroid = np.mean(zCoord)
                # TS/TG split by y, ZS/ZG split by z
                if yCentroid <= yEval:
                    if zCentroid <= zEval:
                        elements['TS_ZS'].append(stressLastFrameElementNodal[ii])
                        elementLabels['TS_ZS'].append(stressLastFrameElementNodal[ii].elementLabel)
                    elif zCentroid >= zEval:
                        elements['TS_ZG'].append(stressLastFrameElementNodal[ii])
                        elementLabels['TS_ZG'].append(stressLastFrameElementNodal[ii].elementLabel)
                elif yCentroid >= yEval:
                    if zCentroid <= zEval:
                        elements['TG_ZS'].append(stressLastFrameElementNodal[ii])
                        elementLabels['TG_ZS'].append(stressLastFrameElementNodal[ii].elementLabel)
                    elif zCentroid >= zEval:
                        elements['TG_ZG'].append(stressLastFrameElementNodal[ii])
                        elementLabels['TG_ZG'].append(stressLastFrameElementNodal[ii].elementLabel)
    return elements, elementLabels

def postProcessFEInner_TS_TG_ZS_ZG(odbPath, stepName, tol_y, tol_z):
    """
    For each FEInner set, group elements by element centroid location (using cartesian_elemDefLocCenGrav_TS_TG_ZS_ZG),
    and for each group (TS_ZS, ...), output unique nodes with averaged nodal stress for that group.
    """
    odbObject = session.openOdb(name=odbPath)
    yEval = plateWidth / 2
    zEval = -plateLength / 2
    stressComps = getStressComponents()
    for layerIndex in range(N):
        setName = "FEInner{}".format(layerIndex + 1)
        elementSet = odbObject.rootAssembly.elementSets[setName.upper()]
        lastFrame = odbObject.steps[stepName].frames[-1]
        # ELEMENT_NODAL stress (list of FieldValue objects)
        stressField = lastFrame.fieldOutputs['S'].getSubset(region=elementSet, position=ELEMENT_NODAL)
        coordField  = lastFrame.fieldOutputs['COORD'].getSubset(region=elementSet, position=ELEMENT_NODAL)
        # Prepare for cartesian_elemDefLocCenGrav_TS_TG_ZS_ZG
        # This function expects stressLastFrameElementNodal (list of FieldValue) and compositeSetElements (elementSet)
        elements_dict, elementLabels_dict = cartesian_elemDefLocCenGrav_TS_TG_ZS_ZG(
            yEval, zEval, stressField.values, elementSet
        )
        # Build a mapping from (nodeLabel, elementLabel) -> coordinates and stress
        nodeElem_to_coords = {}
        nodeElem_to_stress = {}
        for v in coordField.values:
            nodeElem_to_coords[(v.nodeLabel, v.elementLabel)] = v.data
        for s in stressField.values:
            nodeElem_to_stress[(s.nodeLabel, s.elementLabel)] = s.data
        # For each group, build nodeLabel -> [stress, coords] list (collect all for nodes shared by multiple elements)
        for group in ['TS_ZS', 'TS_ZG', 'TG_ZS', 'TG_ZG']:
            node_data = {}
            elemLabelsInGroup = set([fv.elementLabel for fv in elements_dict[group]])
            for s in stressField.values:
                if s.elementLabel not in elemLabelsInGroup:
                    continue
                nodeLabel = s.nodeLabel
                elemLabel = s.elementLabel
                # Only nodes at (x, yEval, zEval) within tolerance
                coords = nodeElem_to_coords.get((nodeLabel, elemLabel), None)
                if coords is None:
                    continue
                x, y, z = coords[0], coords[1], coords[2]
                if abs(y - yEval) < tol_y and abs(z - zEval) < tol_z:
                    if nodeLabel not in node_data:
                        node_data[nodeLabel] = []
                    node_data[nodeLabel].append( (x, y, z, s.data, elemLabel) )
            # Output one line per node (average if node appears multiple times in group)
            nodeData_sorted = []
            for nodeLabel, vals in node_data.items():
                arr = np.array([v[3] for v in vals])  # stress vectors
                avgStress = arr.mean(axis=0)
                # Use the coordinates from the first occurrence (all should be identical for (x, y, z))
                x, y, z = vals[0][0], vals[0][1], vals[0][2]
                # For debugging, you may want to save elemLabels as well (for traceability)
                elemLabel = vals[0][4]
                nodeData_sorted.append( (elemLabel, nodeLabel, x, y, z, avgStress) )
            # Sort by x (thickness)
            nodeData_sorted.sort(key=lambda row: row[2])
            # Write file if there are any nodes in group
            if nodeData_sorted:
                filename = setName + '_' + group + '.txt'
                with open(filename, 'w') as file:
                    for row in nodeData_sorted:
                        elemLabel, nodeLabel, x, y, z, avgStress = row
                        file.write("{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:d}\t{:d}\n".format(
                            avgStress[0], avgStress[1], avgStress[2],
                            avgStress[3], avgStress[4], avgStress[5],
                            x, y, z, nodeLabel, elemLabel
                        ))
    odbObject.close()


# Example usage
stepName = 'Mechanical_Loading'
#postProcessFEInner_TS_TG_centroid(odbPath, stepName)
tol_y = 1e-6  # Tolerance for y-coordinate
tol_z = 1e-6  # Tolerance for z-coordinate
postProcessFEInner_TS_TG_ZS_ZG(odbPath, stepName,tol_y,tol_z)
