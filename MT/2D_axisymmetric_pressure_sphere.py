# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 2024

@author: psingh

Axisymmetric Curved Composite Section
======================================
This macro creates an axisymmetric FE model of a curved quarter-circle 
composite laminate section under internal pressure.

The sketch creates a closed curved cross-section:
- Inner radius: ri
- Outer radius: ra = ri + d (total laminate thickness)
- Angular extent: 0° to 90° (quarter circle)
- Geometry: Two concentric quarter-circle arcs connected by two radial lines

Model Type: AXISYMMETRIC (2D)
Element Type: CAX8R (8-node quadrilateral with reduced integration)
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

# path_modules = 'C:\\Users\\ps25kuhi\\Documents\\01_Promotion\\01_ANA_2.1\\FE'

# os.chdir(path_modules)

# =============================================================================
# COORDINATE TRANSFORMATION FUNCTIONS
# =============================================================================

def pol2cart_x(r, theta):
    """Convert polar to Cartesian X-coordinate."""
    return r * math.cos(theta)

def pol2cart_y(r, theta):
    """Convert polar to Cartesian Y-coordinate."""
    return r * math.sin(theta)

def cart2pol_radius(x, y):
    """Convert Cartesian to polar radius."""
    try:
        return np.sqrt(x*x + y*y)
    except:
        return math.sqrt(x*x + y*y)

def cart2pol_theta(x, y):
    """Convert Cartesian to polar angle."""
    try:
        return np.arctan2(y, x)
    except:
        return math.atan2(y, x)

# Create coordinate transformation module class
class CoordinateTransformation:
    """Coordinate transformation module class."""
    
    @staticmethod
    def pol2cart_x(r, theta):
        return r * math.cos(theta)
    
    @staticmethod  
    def pol2cart_y(r, theta):
        return r * math.sin(theta)
    
    @staticmethod
    def cart2pol_radius(x, y):
        try:
            return np.sqrt(x*x + y*y)
        except:
            return math.sqrt(x*x + y*y)
    
    @staticmethod
    def cart2pol_theta(x, y):
        try:
            return np.arctan2(y, x)
        except:
            return math.atan2(y, x)

ct = CoordinateTransformation()

# analysis_Path = 'C:\\Users\\ps25kuhi\\Documents\\01_Promotion\\01_ANA_2.1\\FE\\pressure_axisymmetric_curved_section'

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

############################################################################
# PARAMETER ZUM ANPASSEN

# Geometrie-Endwinkel (für Curved Section: π/2 = 90°)
angleEnd = np.pi / 2
# Auswertungswinkel für Post-Processing
angleEval = 0.0
Null_Auswaertung = True

# Meshparameter für den gekrümmten Bereich
mR, mT = 3, 250  # Radial elements, Tangential elements
mRRatio, mTRatio = 1, 6

# Probekörpergeometrie (Curved section only)
ri = 180.0  # Innenradius [mm]

compositeMaterialName = 'cfk'
AngleName = 'Hon25'

if compositeMaterialName == 'cfk':
    if AngleName == 'Test':
        ric_plyAngle = '90 90 90 90'.split()
    elif AngleName == 'Hon25':
        ric_plyAngle = '90 90 90 90 15 -15 -15 15 75 -75 -75 75 90 90 90 90 20 -20 -20 20 65 -65 -65 65 90 90 90 90 10 -10 -10 10 90 90 90 90 35 -35 -35 35 35 -35 -35 35 90 90 90 90 10 -10 -10 10 90 90 90 90 65 -65 -65 65 20 -20 -20 20 90 90 90 90 75 -75 -75 75 15 -15 -15 15 90 90 90 90'.split()  
    elif AngleName == 'bestSymm':
        ric_plyAngle = '90 15 75 90 20 65 90 10 90 35 35 90 10 90 65 20 90 75 15 90'.split()

############################################################################
# Bestimmung des Zeitpunktes der FE-Analyse:
analysis_currentDateTime = dt.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
modelName = 'pressure_sphere_' + '_mR' + str(mR) + '_mT' + str(mT) + '_' + analysis_currentDateTime

# Erstellen eines neuen Ordners
# analysis_newPath = analysis_Path + '\\' + modelName
# os.makedirs(analysis_newPath)
# os.chdir(analysis_newPath)

# Transformieren in 0 90 Schreibweise
plyAngle = [float(a) for a in ric_plyAngle if -180 <= float(a) <= 180]
N = len(plyAngle)
iInterfaceEval = list(range(1, N))

# Startwinkel, Öffnungswinkel
angleStart, angleOpening = 0.0, np.pi / 2

if compositeMaterialName == 'cfk':
	E1, E2, E3 = 10000.0, 147000.0, 10000.0
	# Poisson's ratios for transversely isotropic material (fiber in direction 2):
	Nu21 = 0.27  # Radial-Fiber coupling (radial contracts when fiber loaded)
	Nu31 = 0.35  # Radial-Hoop coupling (in-plane transverse, E1=E3)
	Nu23 = 0.27  # Fiber-Hoop coupling (hoop contracts when fiber loaded)
	# Reciprocal ratios (calculated by Abaqus):
	Nu12 = Nu21 * E1/E2 
	Nu13 = Nu31 * E1/E3
	# Nu32 = Nu23 * E3/E2 = 0.27 * 10000/147000 ≈ 0.018 (small - hoop stretches little when fiber loaded)
	G12, G13, G23 = 7000.0, 3700.0, 7000.0  # G12=radial-fiber, G13=radial-hoop, G23=fiber-hoop
	alpha11,alpha22,alpha33=2.6e-5,-1.0e-6,2.6e-5  # alpha11=radial, alpha22=FIBER/axial, alpha33=hoop
	dsingle = 0.35
elif compositeMaterialName == 'gfk':
	E1, E2, E3 = 9552.6, 39296.0, 9552.6
	Nu21, Nu13, Nu23 = 0.29, 0.38, 0.29
	Nu12=E1/E2*Nu21
	G12, G13, G23 = 3080.5, 3449.0, 3080.5
	alpha11,alpha22,alpha33=2.6e-5,8.6e-6,2.6e-5
	dsingle = 0.190
elif compositeMaterialName == 'cfknew':
	# Material coordinate system REDEFINED for STACK_3 compatibility:
	# Global X = radial (r, thickness/stacking) → Material 3 (STACK_3)
	# Global Y = axial (z, FIBER at 0°) → Material 1 
	# Global Z = hoop (θ, circumferential) → Material 2
	# 
	# For 0° ply: fibers align with Global Y (axial) = Material axis 1
	# Axis 1 = axial/FIBER direction (E1 = 147000 MPa) ← STRONG
	# Axis 2 = hoop (circumferential, E2 = 10000 MPa)
	# Axis 3 = radial (through-thickness, E3 = 10000 MPa) ← STACK direction
	E1, E2, E3 = 147000.0, 10000.0, 10000.0
	
	# Poisson's ratios - MUST satisfy reciprocity: ν_ij/E_i = ν_ji/E_j
	# For transversely isotropic material (E2 = E3, fiber in direction 1):
	Nu23 = 0.35  # In-plane transverse Poisson's ratio (hoop-radial), symmetric since E2=E3
	Nu12 = 0.27  # Hoop contraction when loaded in fiber direction
	Nu13 = 0.27  # Radial contraction when loaded in fiber direction (same as Nu12 due to E2=E3)
	
	# Reciprocal Poisson's ratios (calculated by Abaqus):
	# Nu21 = Nu12 * E2/E1 = 0.27 * 10000/147000 = 0.0184 (transverse contracts little when fiber loaded)
	# Nu31 = Nu13 * E3/E1 = 0.27 * 10000/147000 = 0.0184 (same as Nu21)
	# Nu32 = Nu23 * E3/E2 = 0.35 * 10000/10000 = 0.35 (symmetric, E2=E3)
	
	G12, G13, G23 = 7000.0, 7000.0, 3700.0  # G12, G13 = fiber-direction shear; G23 = in-plane shear
	alpha11,alpha22,alpha33=-1.0e-6,2.6e-5,2.6e-5  # α11=FIBER/axial, α22=hoop, α33=radial
	dsingle = 0.7
else:
	pass

# Berechnung verschiedener Punkte
d = dsingle * N
ra = ri + d

# Radien der Mittelflächen
rk = [ri]
for ii in range(N):
    rk.append(rk[ii] + dsingle)

print("rk array:", rk)
print("len(rk):", len(rk))
print("N:", N)

rN = rk[-1]

# Betrag des angreifenden Innen-/Aussendruckes:
OuterPressure, InnerPressure = 0.0, 70.0
if OuterPressure != 0.0 or InnerPressure != 0.0:
    boolPressure = True
else:
    boolPressure = False

# Temperaturdifferenz:
tempDif = 0.0
if tempDif != 0.0:
    boolTempDif = True
else:
    boolTempDif = False

# Parameter functions
def radialGeometryParameters():
    """Berechnet radiale Geometrieparameter für Laminatschichten."""
    rm = [(rk[ii] + rk[ii + 1]) / 2 for ii in range(N)]
    
    print("Creating rEdgesEval...")
    print("iInterfaceEval:", iInterfaceEval)
    
    try:
        rEdgesEval = [rk[iInterEval] for iInterEval in iInterfaceEval]
        print("rEdgesEval created successfully:", rEdgesEval)
    except IndexError as e:
        print("ERROR creating rEdgesEval:", e)
        rEdgesEval = []
    
    try:
        rInterfaceEval = [rk[iInterEval] for iInterEval in iInterfaceEval]
        print("rInterfaceEval created successfully:", rInterfaceEval)
    except IndexError as e:
        print("ERROR creating rInterfaceEval:", e)
        rInterfaceEval = []
    
    rEdgesAll = []
    for ii in range(N + 1):
        rEdgesAll.append(rk[ii])
    
    rFacesAll = [(rEdgesAll[ii] + rEdgesAll[ii + 1]) / 2 for ii in range(N)]
    
    return (rm, rEdgesEval, rEdgesAll, rInterfaceEval, rFacesAll)

def circGeometryParameters():
    """Berechnet Umfangswinkel-Geometrieparameter."""
    thetaEdges = [angleStart, angleEnd]
    
    thetaPartition = [angleEnd]
    
    refinementInterval = np.pi * 10 / 180
    thetaPartiotionsAll = [
        angleStart,
        angleEnd - refinementInterval,
        angleEnd]
    thetaPartitionFaces = [(thetaPartiotionsAll[ii] + thetaPartiotionsAll[ii + 1]) / 2 for ii in range(len(thetaPartiotionsAll) - 1)]
    thetaFaces = [(thetaEdges[ii] + thetaEdges[ii + 1]) / 2 for ii in range(len(thetaEdges) - 1)]
    
    return (thetaPartition, thetaEdges, thetaFaces, thetaPartitionFaces)

rm, rEdgesEval, rEdgesAll, rInterfaceEval, rFacesAll = radialGeometryParameters()
thetaPartition, thetaEdges, thetaFaces, thetaPartitionFaces = circGeometryParameters()

# Dateien für Auswärtung
Ratio = (ri + d) / ri

def export_matrix(matrix, filename, delimiter='\t'):
    """Exportiert eine Matrix in eine Datei."""
    with open(filename, 'w') as file:
        for row in matrix:
            formatted_row = delimiter.join(
                "{:.0f}".format(x) if x.is_integer() else "{:.2f}".format(x) for x in row
            )
            file.write(formatted_row + '\n')

def materialParametersExport():
    """Exportiert die Materialparameter."""
    matConfigCompositeMatrix = np.array([
        [E1, E2, E3],
        [G12, G13, G23],
        [Nu12, Nu13, Nu23]
    ])
    
    delimiter = '\t'
    
    with open('materialParameters', 'w') as file:
        for row in matConfigCompositeMatrix:
            formatted_row = delimiter.join("{:.2f}".format(x).lstrip('0') if x != 0 else '0' for x in row)
            file.write(formatted_row + '\n')

materialParametersExport()

def export_geometry(plyAngle, MDetailed, orderLagrangeDetailed, Ratio, dL, width_h, angleBegin, openingAngle, angleEnd):
    """Exports geometry and ply angle information for curved section.""",
    N = len(plyAngle)
    pAExport = np.zeros((1, N))
    for ii in range(N):
        pAExport[0, ii] = plyAngle[ii]
    export_matrix(pAExport, 'compositeLayup', delimiter='\t')
    
    geometryCompositeExport = np.zeros((1, 8))
    geometryCompositeExport[0, 0] = MDetailed
    geometryCompositeExport[0, 1] = orderLagrangeDetailed
    geometryCompositeExport[0, 2] = Ratio
    geometryCompositeExport[0, 3] = float(dL)
    geometryCompositeExport[0, 4] = width_h
    geometryCompositeExport[0, 5] = 180 * angleBegin / np.pi
    geometryCompositeExport[0, 6] = 180 * openingAngle / np.pi
    geometryCompositeExport[0, 7] = 180 * angleEnd / np.pi
    export_matrix(geometryCompositeExport, 'compositeGeometry', delimiter='\t')

MDetailed = 8.0
orderLagrangeDetailed = 3
dL = dsingle
width_h = 10
angleBegin = 0

export_geometry(plyAngle, MDetailed, orderLagrangeDetailed, Ratio, dL, width_h, angleBegin, angleOpening, angleEnd)

#--------------------------------------------------------------
# SKIZZE ERSTELLEN - NUR GEKRÜMMTER BEREICH
Mdb()
mdb.models.changeKey(fromName='Model-1', toName=modelName)
mdb.models[modelName].rootAssembly.clearGeometryCache()
mdb.models[modelName].rootAssembly.regenerate()
LModel = mdb.models[modelName]

def create_closed_curved_section_sketch(LModel, rk, thetaEdges, sketch_name='Curved-Section-Closed'):
    """
    Creates a closed curved sketch for the quarter-circle composite section.
    
    This function creates a quarter-circle annular sector with layered
    composite structure. The sketch is a closed path with 2 concentric
    arcs and 2 radial straight lines.
    
    For axisymmetric models, this represents a quarter-circle cross-section
    from 0° to 90° (pi/2 radians).
    
    Args:
        LModel: Abaqus model object
        rk (list): Interface radii [inner_radius, ..., outer_radius]
        thetaEdges (list): Theta angles for the arcs [angleStart, angleEnd]
        sketch_name (str): Name for the sketch
        
    Returns:
        tuple: (LModelSketch, OuterCirc_id, InnerCirc_id)
    """
    
    print("Creating closed curved section sketch...")
    
    r_inner = rk[0]   # Inner radius
    r_outer = rk[-1]  # Outer radius
    theta_start = thetaEdges[0]  # Start angle (0)
    theta_end = thetaEdges[1]    # End angle (pi/2)
    
    print("Geometry parameters: r_inner={}, r_outer={}, theta_start={}, theta_end={}".format(
        r_inner, r_outer, theta_start, theta_end))
    
    # Create the sketch
    LModelSketch = LModel.ConstrainedSketch(name=sketch_name, sheetSize=(2 * r_outer))
    
    # Add construction line for Y-axis (axis of symmetry for axisymmetric analysis)
    print("Adding construction line for Y-axis (axis of symmetry)...")
    LModelSketch.ConstructionLine(point1=(0.0, -r_outer - 10.0), point2=(0.0, r_outer + 10.0))
    print("SUCCESS: Construction line added for axis of symmetry")
    
    # Create CLOSED curved path for quarter-circle annular sector
    print("Creating complete closed curved profile...")
    
    # THE SKETCH: Two concentric quarter-circle arcs connected by two radial lines
    # Starting from inner radius at theta=0, going counterclockwise:
    
    # 1. Radial line at theta=0 (from inner to outer radius along X-axis)
    print("1. Creating first radial line at theta=0")
    LModelSketch.Line(point1=(r_inner, 0.0), point2=(r_outer, 0.0))
    
    # 2. Outer arc: quarter-circle from (r_outer, 0) to (0, r_outer)
    outer_start = (r_outer, 0.0)
    outer_end = (ct.pol2cart_x(r_outer, theta_end), ct.pol2cart_y(r_outer, theta_end))
    
    print("2. Creating outer arc: start={}, end={}".format(outer_start, outer_end))
    OuterCirc = LModelSketch.ArcByCenterEnds(
        center=(0.0, 0.0), 
        point1=outer_start, 
        point2=outer_end, 
        direction=COUNTERCLOCKWISE)
    OuterCirc_id = OuterCirc.id
    
    # 3. Radial line at theta=theta_end (from outer to inner radius along Y-axis)
    print("3. Creating second radial line at theta={}".format(theta_end))
    LModelSketch.Line(
        point1=(ct.pol2cart_x(r_outer, theta_end), ct.pol2cart_y(r_outer, theta_end)), 
        point2=(ct.pol2cart_x(r_inner, theta_end), ct.pol2cart_y(r_inner, theta_end)))
    
    # 4. Inner arc: quarter-circle from (0, r_inner) to (r_inner, 0) - CLOCKWISE to close
    inner_start = (ct.pol2cart_x(r_inner, theta_end), ct.pol2cart_y(r_inner, theta_end))
    inner_end = (r_inner, 0.0)
    
    print("4. Creating inner arc: start={}, end={}".format(inner_start, inner_end))
    InnerCirc = LModelSketch.ArcByCenterEnds(
        center=(0.0, 0.0), 
        point1=inner_start, 
        point2=inner_end, 
        direction=CLOCKWISE)
    InnerCirc_id = InnerCirc.id
    
    print("SUCCESS: Created CLOSED curved section (quarter-circle annular sector)")
    print("  - 2 quarter-circle arcs (concentric, centered at origin)")
    print("  - 2 radial connecting lines")
    print("  - Total: 4 geometric entities forming ONE CLOSED REGION")
    print("  - Angular extent: {} to {} radians ({} to {} degrees)".format(
        theta_start, theta_end, 
        math.degrees(theta_start), math.degrees(theta_end)))
    
    return LModelSketch, OuterCirc_id, InnerCirc_id

def sketchLModel():
    """
    Erstellt die Skizze des gekrümmten Bereichs.
    
    Returns:
        tuple:
            - LModelPart: Die erstellte Part-Instanz
            - OuterCirc_id: ID des äußeren Bogens
            - InnerCirc_id: ID des inneren Bogens
    """
    
    LModelSketch, OuterCirc_id, InnerCirc_id = create_closed_curved_section_sketch(
        LModel, rk, thetaEdges, 'Curved-Section-Sketch')
    
    # Create the part using AXISYMMETRIC
    LModelPart = LModel.Part(name='Curved-Section-Part', dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY)
    
    try:
        LModelPart.BaseShell(sketch=LModelSketch)
        print("SUCCESS: BaseShell operation completed successfully!")
    except Exception as e:
        print("BaseShell failed: {}".format(str(e)))
        raise e
    
    return (LModelPart, OuterCirc_id, InnerCirc_id)

LModelPart, OuterCirc_id, InnerCirc_id = sketchLModel()

#--------------------------------------------------------------
# Abkürzungen als Variablen speichern
LModelAssembly = LModel.rootAssembly
LModelPart = mdb.models[modelName].parts['Curved-Section-Part']
LModelFaces = LModelPart.faces

def instanceLModel():
    """Erstellt Instanzen."""
    LModelinstance = LModelAssembly.Instance(name='Curved_Section_Instance', part=LModelPart, dependent=ON)
    return LModelinstance

LModelInstance = instanceLModel()

#--------------------------------------------------------------
# Partionierung

def partitionComposite():
    """Erstellt Partitionen basierend auf Lagendicke."""
    print("\n=== Starting partitionComposite() for curved section ===")
    
    try:
        LModelPartition = mdb.models[modelName].parts['Curved-Section-Part']
        print("Successfully retrieved part: Curved-Section-Part")
    except Exception as e:
        print("ERROR: Failed to retrieve part - {}".format(str(e)))
        return
    
    try:
        partitionfaces = LModelPartition.faces
        print("Number of faces in part: {}".format(len(partitionfaces)))
        if len(partitionfaces) > 0:
            print("Face[0] coordinates: {}".format(partitionfaces[0]))
        else:
            print("ERROR: No faces found in part!")
            return
    except Exception as e:
        print("ERROR: Failed to get faces - {}".format(str(e)))
        return
    
    try:
        print("Creating sketch transform...")
        LModelTransform = LModelPartition.MakeSketchTransform(
            sketchPlane=partitionfaces[0],
            sketchPlaneSide=SIDE1,
            origin=(0, 0, 0))
        print("SUCCESS: Sketch transform created")
    except Exception as e:
        print("ERROR: Failed to create sketch transform - {}".format(str(e)))
        return
    
    try:
        print("Creating partition sketch...")
        print("  Sheet size: {}".format(2 * ra))
        LModelPartitionSketch = mdb.models[modelName].ConstrainedSketch(
            name='Curved-Partition-Sketch',
            sheetSize=(2 * ra),
            transform=LModelTransform)
        print("SUCCESS: Partition sketch created")
    except Exception as e:
        print("ERROR: Failed to create partition sketch - {}".format(str(e)))
        return
    
    try:
        g = LModelPartitionSketch.geometry
        print("Sketch geometry retrieved, number of entities: {}".format(len(g)))
    except Exception as e:
        print("ERROR: Failed to get sketch geometry - {}".format(str(e)))
        return
    
    try:
        print("Projecting references onto sketch...")
        LModelPartition.projectReferencesOntoSketch(sketch=LModelPartitionSketch, filter=COPLANAR_EDGES)
        print("SUCCESS: References projected")
        print("Geometry entities after projection: {}".format(len(g)))
    except Exception as e:
        print("ERROR: Failed to project references - {}".format(str(e)))
        return
    
    # Create offsets for each layer - FIXED APPROACH
    # Each partition needs its own sketch, we cannot reuse the same sketch
    print("\nCreating layer partitions...")
    print("Number of layers (N): {}".format(N))
    print("Number of partitions to create (N-1): {}".format(N - 1))
    print("rk values: {}".format(rk))
    print("rm values: {}".format(rm))
    print("dsingle: {}".format(dsingle))
    
    for ii in range(N - 1):
        print("\n--- Layer {} partition (iteration {}/{}) ---".format(ii + 1, ii + 1, N - 1))
        offset_distance = (ii + 1) * dsingle
        print("  Offset distance: {}".format(offset_distance))
        
        # CRITICAL FIX: Create a NEW sketch for each partition iteration
        try:
            print("  Creating fresh partition sketch for this layer...")
            
            # Get updated face list (important after each partition)
            partitionfaces = LModelPartition.faces
            print("  Current number of faces: {}".format(len(partitionfaces)))
            
            # CRITICAL FIX: Find the face that CONTAINS the next partition line radius
            # The partition line will be at rk[ii+1], so we need a face that spans across this radius
            target_face = None
            if len(partitionfaces) == 1:
                target_face = partitionfaces[0]
                print("  Using the only available face")
            else:
                # Find the face that contains a point at radius rk[ii+1]
                # We search for a face that contains the partition line location
                partition_radius = rk[ii + 1]
                test_point = (ct.pol2cart_x(partition_radius, thetaFaces[0]), 
                              ct.pol2cart_y(partition_radius, thetaFaces[0]), 0)
                print("  Searching for face containing partition radius {} at point {}".format(partition_radius, test_point))
                
                try:
                    target_face = partitionfaces.findAt(test_point)
                    print("  Found target face containing radius {}".format(partition_radius))
                except Exception as e:
                    print("  ERROR: Could not find face at partition location - {}".format(str(e)))
                    print("  Trying to find largest face (unpartitioned region)...")
                    # Fallback: find the face with largest radius range (the unpartitioned outer region)
                    for face in partitionfaces:
                        face_point = face.pointOn[0]
                        face_radius = ct.cart2pol_radius(face_point[0], face_point[1])
                        print("    Face has pointOn radius: {}".format(face_radius))
                        if face_radius > rk[ii]:  # Face must be beyond already-partitioned layers
                            target_face = face
                            print("  Selected face with radius {} (beyond rk[{}]={})".format(face_radius, ii, rk[ii]))
                            break
                
                if target_face is None:
                    target_face = partitionfaces[0]
                    print("  WARNING: Could not identify correct face, using partitionfaces[0]")
            
            # Create new transform using the target face
            LModelTransform_new = LModelPartition.MakeSketchTransform(
                sketchPlane=target_face,
                sketchPlaneSide=SIDE1,
                origin=(0, 0, 0))
            
            # Create new sketch with unique name
            sketch_name = 'Partition-Sketch-Layer-{}'.format(ii + 1)
            LModelPartitionSketch_new = mdb.models[modelName].ConstrainedSketch(
                name=sketch_name,
                sheetSize=(2 * ra),
                transform=LModelTransform_new)
            
            g_new = LModelPartitionSketch_new.geometry
            
            # Project references onto this new sketch
            LModelPartition.projectReferencesOntoSketch(sketch=LModelPartitionSketch_new, filter=COPLANAR_EDGES)
            print("  Fresh sketch created with {} geometry entities".format(len(g_new)))
            
        except Exception as e:
            print("  ERROR: Failed to create fresh sketch - {}".format(str(e)))
            import traceback
            print("  Traceback: {}".format(traceback.format_exc()))
            continue
        
        # ALTERNATIVE APPROACH: Draw a concentric arc directly instead of offset
        # This is more reliable for creating partition arcs in curved geometry
        try:
            print("  Drawing partition arc at radius {}...".format(rk[ii + 1]))
            
            # Calculate start and end points for the partition arc
            partition_radius = rk[ii + 1]
            arc_start = (partition_radius, 0.0)  # At theta=0
            arc_end = (ct.pol2cart_x(partition_radius, angleEnd), 
                      ct.pol2cart_y(partition_radius, angleEnd))  # At theta=angleEnd
            
            print("  Arc endpoints: start={}, end={}".format(arc_start, arc_end))
            
            # Draw concentric arc from theta=0 to theta=angleEnd at radius rk[ii+1]
            LModelPartitionSketch_new.ArcByCenterEnds(
                center=(0.0, 0.0),
                point1=arc_start,
                point2=arc_end,
                direction=COUNTERCLOCKWISE)
            
            print("  SUCCESS: Partition arc drawn")
        except Exception as e:
            print("  ERROR: Failed to draw partition arc - {}".format(str(e)))
            import traceback
            print("  Traceback: {}".format(traceback.format_exc()))
            continue
        
        # CRITICAL FIX: Use the SAME target_face that we used for sketch creation
        # This face already contains the partition arc we just created
        print("  Using target_face for partition (same face used for sketch creation)")
        LModelPickedFaces = target_face
        
        # Apply partition immediately with this fresh sketch
        try:
            print("  Applying partition by sketch...")
            LModelPartition.PartitionFaceBySketch(faces=LModelPickedFaces, sketch=LModelPartitionSketch_new)
            print("  SUCCESS: Partition applied for layer {}".format(ii + 1))
        except Exception as e:
            print("  ERROR: Failed to partition face - {}".format(str(e)))
            import traceback
            print("  Traceback: {}".format(traceback.format_exc()))
            continue
    
    print("\n=== Radial partitions completed ===\n")
    
    # Add theta partitions for refined analysis
    angleLimit = 30.0 * np.pi / 180.0
    if Null_Auswaertung:
        angleEval_partition = 0
    else:
        angleEval_partition = angleEval
    
    thetaPartitionStart = max(0.0, angleEval_partition - angleLimit)
    thetaPartitionEnd = min(angleEnd, angleEval_partition + angleLimit)
    
    try:
        for theta in [thetaPartitionStart, thetaPartitionEnd]:
            if theta != 0.0 and theta != angleEnd:
                LModelPartition.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[0], theta), ct.pol2cart_y(rk[0], theta), 0.0))
                LModelPartition.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1], theta), ct.pol2cart_y(rk[-1], theta), 0.0))
                
                datumKeys = list(LModelPartition.datums.keys())
                datumKeys.sort()
                
                LModelPartition.DatumPlaneByThreePoints(
                    point1=LModelPartition.datums[datumKeys[-2]],
                    point2=LModelPartition.datums[datumKeys[-1]], 
                    point3=(0.0, 0.0, 0.0))
                
                planeKeys = [key for key in LModelPartition.features.keys() if key.startswith('Datum plane-')]
                if planeKeys:
                    planeKeys.sort()
                    latestPlaneKey = planeKeys[-1]
                    planeId = LModelPartition.features[latestPlaneKey].id
                    LModelPartition.PartitionFaceByDatumPlane(faces=partitionfaces, datumPlane=LModelPartition.datums[planeId])
    except:
        print("Warning: Could not create theta partitions")

partitionComposite()

#--------------------------------------------------------------
def materialComposite():
    """Definiert Materialeigenschaften für Composite-Materialien."""
    def transfoRAxis_R(phi):
        """Transformation matrix R(phi)."""
        cos_phi = math.cos(math.radians(phi))
        sin_phi = math.sin(math.radians(phi))
        return np.array([
            [1, 0, 0, 0, 0, 0],
            [0, cos_phi**2, sin_phi**2, 2*cos_phi*sin_phi, 0, 0],
            [0, sin_phi**2, cos_phi**2, -2*cos_phi*sin_phi, 0, 0],
            [0, -cos_phi*sin_phi, cos_phi*sin_phi, (cos_phi**2)-(sin_phi**2), 0, 0],
            [0, 0, 0, 0, cos_phi, -sin_phi],
            [0, 0, 0, 0, sin_phi, cos_phi]
        ])
    
    S0 = np.zeros((6, 6))
    S0[0, 0] = 1 / E1
    S0[1, 1] = 1 / E2
    S0[2, 2] = 1 / E3
    S0[3, 3] = 1 / G23
    S0[4, 4] = 1 / G13
    S0[5, 5] = 1 / G12
    S0[0, 1] = S0[1, 0] = -Nu12 / E1
    S0[0, 2] = S0[2, 0] = -Nu13 / E1
    S0[1, 2] = S0[2, 1] = -Nu23 / E2
    
    for phi in range(-180, 181):
        R_phi = transfoRAxis_R(phi)
        Sr_phi = np.dot(np.dot(R_phi.T, S0), R_phi)
        
        material_name = "{}_{}".format(compositeMaterialName, phi)
        compositeMaterial = LModel.Material(material_name)
        compositeMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=(
            (1 / Sr_phi[0, 0], 1 / Sr_phi[1, 1], 1 / Sr_phi[2, 2],
             -Sr_phi[0, 1] / Sr_phi[0, 0], -Sr_phi[0, 2] / Sr_phi[0, 0], -Sr_phi[1, 2] / Sr_phi[1, 1],
             1 / Sr_phi[3, 3], 1 / Sr_phi[4, 4], 1 / Sr_phi[5, 5]),))
        compositeMaterial.Expansion(type=ORTHOTROPIC, table=((alpha11, alpha22, alpha33),))
        
        section_name = "Composite_Section_{:.1f}".format(phi)
        LModel.HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)

materialComposite()

#--------------------------------------------------------------
# Koordinatensystem festlegen
# For curved section: SPHERICAL coordinate system
LModelPart.DatumCsysByThreePoints(origin=(0, 0, 0), point1=(1, 0, 0.0), point2=(0, 1.0, 0.0), name='Cylindrical-KOS-Orientation', coordSysType=SPHERICAL)
orientationKOSid = LModelPart.features['Cylindrical-KOS-Orientation'].id

# KOS Für Boundary Conditions
LModelAssembly.DatumCsysByThreePoints(name='BoundaryCD-Datum', coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), 
                                      point1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
BoundaryCDDatum_ID = LModelAssembly.features['BoundaryCD-Datum'].id
BoundaryCDDatum = LModelAssembly.datums[BoundaryCDDatum_ID]

#--------------------------------------------------------------
def CompositeOrientation():
    """Weist Materialorientierung basierend auf plyAngle zu."""
    for ii in range(N):
        for mm in range(len(thetaFaces)):
            curvedCompositeFaces = LModelPart.faces.findAt(((ct.pol2cart_x(rFacesAll[ii], thetaFaces[mm]), ct.pol2cart_y(rFacesAll[ii], thetaFaces[mm]), 0.0),))
            curvedCompositeRegion = regionToolset.Region(faces=curvedCompositeFaces)
            LModelPart.SectionAssignment(region=curvedCompositeRegion, sectionName='Composite_Section_{:.1f}'.format(plyAngle[ii]), offsetType=MIDDLE_SURFACE)
            LModelPart.MaterialOrientation(region=curvedCompositeRegion, orientationType=SYSTEM, localCsys=LModelPart.datums[orientationKOSid], axis=AXIS_3, stackDirection=STACK_3)

CompositeOrientation()

#--------------------------------------------------------------
# Meshgenerierung

def meshLPart(mR, mT, mRRatio, mTRatio):
    """Erstellt ein Netz für den gekrümmten Bereich."""
    # Radial edges (at constant theta)
    for ii in range(len(rFacesAll)):
        for jj in range(len(thetaEdges)):
            curvedCompositeMeshEdge = LModelPart.edges.findAt(
                ((ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj]),
                  ct.pol2cart_y(rFacesAll[ii], thetaEdges[jj]), 0.0),)
            )
            LModelPart.seedEdgeByNumber(
                edges=curvedCompositeMeshEdge,
                number=mR,
                constraint=FINER
            )
    
    # Tangential edges (arcs at constant radius)
    for ii in range(len(rk)):
        curvedCompositeMeshEdge = LModelPart.edges.findAt(
            ((ct.pol2cart_x(rk[ii], thetaFaces[-1]),
              ct.pol2cart_y(rk[ii], thetaFaces[-1]), 0.0),)
        )
        LModelPart.seedEdgeByNumber(
            edges=curvedCompositeMeshEdge,
            number=mT * 2,
            constraint=FINER
        )
    
    compositeElementType = mesh.ElemType(elemCode=CAX8R, elemLibrary=STANDARD)
    LModelPart.setMeshControls(
        regions=LModelPart.faces,
        elemShape=QUAD,
        technique=STRUCTURED
    )
    LModelPart.setElementType(
        regions=(LModelPart.faces[:],),
        elemTypes=(compositeElementType,)
    )
    LModelPart.seedPart(size=0.4, deviationFactor=0.1, constraint=FINER)
    LModelPart.generateMesh()

meshLPart(mR, mT, mRRatio, mTRatio)

#--------------------------------------------------------------
# Sets zur Auswärtung

if Null_Auswaertung == True:
    angleEval_FEInnner = 0
else:
    angleEval_FEInnner = angleEval

def setsCurvedCompositeCirc(rFacesAll, angleEval_FEInnner, iInterfaceEval, thetaPartitionFaces, angleEnd, rk):
    """Erstellt Sets zur Auswertung."""
    # Sets zur Verifizierung der Inner-Solution:
    for ii in range(len(rFacesAll)):
        try:
            compositeSetEvalEdge = LModelInstance.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], angleEval_FEInnner), ct.pol2cart_y(rFacesAll[ii], angleEval_FEInnner), 0.0),))
            LModelAssembly.Set(edges=compositeSetEvalEdge, name='FEInner' + str(ii + 1))
        except:
            print("Warning: Could not create FEInner set for layer", ii + 1)
    
    # Sets für Interface-Analyse
    for kk in range(len(iInterfaceEval)):
        ii = iInterfaceEval[kk]
        try:
            midAngle = angleEval_FEInnner
            curvedCompositeEdgesLeft = LModelInstance.edges.findAt(((ct.pol2cart_x(rk[ii], midAngle + 0.00001), ct.pol2cart_y(rk[ii], midAngle + 0.00001), 0.0),))
            LModelAssembly.Set(edges=curvedCompositeEdgesLeft, name='FEInterfaceCircCurved' + str(iInterfaceEval[kk]))
        except:
            print("Warning: Could not find curved edge for interface", iInterfaceEval[kk])

setsCurvedCompositeCirc(rFacesAll, angleEval_FEInnner, iInterfaceEval, thetaPartitionFaces, angleEnd, rk)

#--------------------------------------------------------------
def stepMechanicalLoads():
    stepPrevious, step, stepDescription = 'Initial', 'Mechanical_Loading', 'Introduce mechanical loadings'
    LModel.StaticStep(name=step, previous=stepPrevious, description=stepDescription, nlgeom=OFF, initialInc=1.0)
    return (stepPrevious, step)

stepPrevious, step = stepMechanicalLoads()

#---------------------------------------------------------------
# Randbedingungen:
def boundaryConditionSurfaceLoad():
    # Fix radial displacement at theta=0
    edge_coords = []
    for ii in range(len(rFacesAll)):
        coord = (rFacesAll[ii], 0.0, 0.0)
        edge_coords.append((coord,))
    
    curvedCompositeFixCircFaces = LModelInstance.edges.findAt(*edge_coords)
    curvedCompositeFixCircSet = LModelAssembly.Set(edges=curvedCompositeFixCircFaces, name='Set_Fix_U_SS_1')
    LModel.DisplacementBC(name='Fix_U_SS_1', createStepName=stepPrevious, region=curvedCompositeFixCircSet, u2=0.0)
    
    # Fix radial displacement at theta=pi/2
    edge_coords = []
    for ii in range(len(rFacesAll)):
        coord = (ct.pol2cart_x(rFacesAll[ii], thetaEdges[1]), ct.pol2cart_y(rFacesAll[ii], thetaEdges[1]), 0.0)
        edge_coords.append((coord,))
    
    curvedCompositeFixCircFaces = LModelInstance.edges.findAt(*edge_coords)
    curvedCompositeFixCircSet = LModelAssembly.Set(edges=curvedCompositeFixCircFaces, name='Set_Fix_U_SS_2')
    LModel.DisplacementBC(name='Fix_U_SS_2', createStepName=stepPrevious, region=curvedCompositeFixCircSet, u1=0.0)

if boolPressure:
    boundaryConditionSurfaceLoad()

def loadsSurfaceLoad():
    if OuterPressure != 0.0:
        # Outer pressure
        curvedCompositeStringFaces = ''
        for jj in range(len(thetaFaces)):
            curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[-1], thetaFaces[jj]), ct.pol2cart_y(rk[-1], thetaFaces[jj]), 0.0),),) + ','
        
        curvedCompositeStringFacesExec = 'curvedCompositeOuterPressureFaces = LModelInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
        exec(curvedCompositeStringFacesExec)
        
        LModel.Pressure(name='Load_OuterPressure', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeOuterPressureFaces), magnitude=OuterPressure, distributionType=UNIFORM)
    
    if InnerPressure != 0.0:
        # Inner pressure
        angleLimit = 30.0 * np.pi / 180.0
        
        inner_edge_coords = []
        
        # Edge below angleLimit (in the small angle region)
        inner_edge_coords.append(((ct.pol2cart_x(rk[0], angleLimit / 2), 
                                   ct.pol2cart_y(rk[0], angleLimit / 2), 
                                   0.0),))
        
        # Edges above angleLimit (in the main curved region)
        for jj in range(len(thetaFaces)):
            if thetaFaces[jj] > angleLimit:
                inner_edge_coords.append(((ct.pol2cart_x(rk[0], thetaFaces[jj]),
                                           ct.pol2cart_y(rk[0], thetaFaces[jj]),
                                           0.0),))
        
        curvedCompositeInnerPressureFaces = LModelInstance.edges.findAt(*inner_edge_coords)
        LModel.Pressure(name='Load_InnerPressure', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeInnerPressureFaces), magnitude=InnerPressure, distributionType=UNIFORM)

if boolPressure:
    loadsSurfaceLoad()

#---------------------------------------------------------------
# Field Output:
LModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'COORD'))

#-------------------------------------------------------------------
# Speichern der CAE-Datei:
mdb.saveAs(pathName=analysis_newPath + '\\' + modelName + '_CAE_Model')

#--------------------------------------------------------------
# Job
jobName = modelName + '_Job'

def LmodelJob():
    """Erstellt und führt einen Job aus."""
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
        memory=99, 
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
    return simTime

simTime = LmodelJob()

print("=== SIMULATION COMPLETE ===")
print("Model: Curved Quarter-Circle Composite Section")
print("Simulation time: {:.2f} seconds".format(simTime))
compositeModelOdbPath = jobName + '.odb'
compositeModelOdbObject = session.openOdb(name=compositeModelOdbPath)
session.viewports['Viewport: 1'].setValues(displayedObject=compositeModelOdbObject)
compositeModelViewport = session.viewports['Viewport: 1']

# Create spherical coordinate system for postprocessing:
# For axisymmetric curved section: SPHERICAL coordinates
# r-direction: (1,0,0), phi-direction (azimuthal): (0,1,0)
postProc_Spherical_CS = compositeModelOdbObject.rootAssembly.DatumCsysByThreePoints(
    coordSysType=SPHERICAL, 
    name='Spherical_CS_PostProcessing',
    point1=(1, 0, 0),
    point2=(0, 1, 0),
    origin=(0, 0, 0))

# Apply spherical coordinate system to viewport for stress evaluation
compositeModelViewport.odbDisplay.basicOptions.setValues(
    transformationType=USER_SPECIFIED, 
    datumCsys=postProc_Spherical_CS)
