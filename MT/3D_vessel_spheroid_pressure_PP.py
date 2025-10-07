# -*- coding: utf-8 -*-
#---------------------------------------------------------------
# Python-Packages:
import __main__
import os
import time
import datetime as dt
import numpy as np
import operator
import math
# Abaqus-Packages:
from abaqus import *
from abaqusConstants import *
from abaqus import mdb, session, Mdb
from abaqusConstants import XZPLANE, YZPLANE, YAXIS, SWEEP, FINER, C3D8, STANDARD
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

# Ensure constants are available
from abaqusConstants import ON, CYLINDRICAL, CARTESIAN, SPHERICAL, SIDE1, COPLANAR_EDGES, REVERSE, THREE_D, DEFORMABLE_BODY, COUNTERCLOCKWISE, DELETE, GEOMETRY
from abaqusConstants import ENGINEERING_CONSTANTS, SOLID, FROM_SECTION, SPECIFY_ORIENT, ROTATION_NONE, AXIS_1, AXIS_2, AXIS_3, SPECIFY_THICKNESS, STACK_1, STACK_2, STACK_3, ROTATION_ANGLE, SYSTEM
from abaqusConstants import XYPLANE, XAXIS, UNIFORM, UNSET, OFF, GENERAL
from abaqusConstants import SHEAR, NORMAL, TRACTION
from abaqusConstants import DOMAIN  # For job parallelization
from abaqusConstants import ELEMENT_NODAL, INTEGRATION_POINT  # For post-processing
from abaqusConstants import STRUCTURED, FREE  # Additional mesh techniques

path_modules = 'D:\\psingh\\MT\\ABAQUS\\MT'
os.chdir(path_modules)

# Further packages:
import coordinateTransformation_ellipse as ct
#---------------------------------------------------------------
# Names des Modells:
modelName = '0_90_S'

# Modellparameter:
# Schichtwinkel der einzelnen physikalischen Schichten:
plyAngle = [0,90,90,0]

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
r0 = 7.5#rt*h-h/2
eccentricity = 1.0
b0 = r0 * eccentricity
phiStart,phiOpening = 0.0,np.pi/2
phiEdges = [phiStart, phiOpening]
height = phiOpening*r0/2
thetaStart = -np.pi/4  # Starting angle for revolution
thetaEnd = np.pi/4  # Total angle of revolution (previously thetaEnd)
thetaOpening = thetaEnd - thetaStart
angleEval = np.pi/2 + thetaStart + thetaOpening / 2.0

mRadial=4
mAxial=16
mMeridional=int(eccentricity*mAxial)
mCirc=12


# Partitioning parameters (number of partitions in each direction)
nCirc = 4    # Number of circumferential partitions has to be odd since cannot coincide with a pre defined plane
nAxial = 2   # Number of axial/height partitions
nMeridional = 6#4*nAxial # Number of shell/radial partitions
# User-defined bias parameters for partitioning directions
theta_bias = 1.0  # Circumferential bias (1.0 = uniform, >1.0 = bias toward thetaEnd)
axial_bias = 1.0  # Axial/height bias (1.0 = uniform, >1.0 = bias toward y=0)
meridional_bias = 1.0  # Shell/radial bias (1.0 = uniform, >1.0 = bias toward x=bk[-1])




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
OuterPressure, InnerPressure = 0.0, 1.0
if OuterPressure != 0.0 or InnerPressure != 0.0:
    boolPressure = True
else:
    boolPressure = False

# Start der FE-(Konvergenz-)Analyse:
#---------------------------------------------------------------
def export_matrix(matrix, filename, delimiter='\t'):
    """Export a matrix to a file with specified delimiter"""
    with open(filename, 'w') as file:
        for row in matrix:
            formatted_row = delimiter.join(
                "{:.0f}".format(x) if isinstance(x, (int, float)) and x == int(x) 
                else "{:.2f}".format(x) for x in row
            )
            file.write(formatted_row + '\n')

def export_composite_layup(plyAngle):
    """Export the composite layup configuration"""
    N = len(plyAngle)
    pAExport = np.zeros((1, N))
    for ii in range(N):
        pAExport[0, ii] = plyAngle[ii]
    export_matrix(pAExport, 'compositeLayup', delimiter='\t')

# Definiere eine neue Model-Database:
Mdb()
session.viewports['Viewport: 1'].setValues(displayedObject=None)

# Definiere das Abaqus-Model:
mdb.models.changeKey(fromName='Model-1', toName=modelName)
compositeModel = mdb.models[modelName]

def radialGeometryParameters(N,rk,bk,iInterfaceEval):
    # Radien der Mittelflaechen der Laminat-Einzelschichten:
    rm = [(rk[ii]+rk[ii+1])/2 for ii in range(N)]
    bm = [(bk[ii]+bk[ii+1])/2 for ii in range(N)]
    # Radien aller Interfaces:
    rInterfaceEval =[rk[iInterEval] for iInterEval in iInterfaceEval]
    bInterfaceEval =[bk[iInterEval] for iInterEval in iInterfaceEval]
    # Radien aller Partitionierungsinterfaces:
    rInterfaceAll =[]
    for ii in range(N+1):
        if ii < N:
            rInterfaceAll.append(rk[ii])
            rInterfaceAll.append(rm[ii])
        else:
            rInterfaceAll.append(rk[ii])
    bInterfaceAll =[]
    for ii in range(N+1):
        if ii < N:
            bInterfaceAll.append(bk[ii])
            bInterfaceAll.append(bm[ii])
        else:
            bInterfaceAll.append(bk[ii])
    # Mittlere Radien aller partitionierten physikalischen Schichten:
    rFacesAll = [(rInterfaceAll[ii]+rInterfaceAll[ii+1])/2 for ii in range(2*N)]
    bFacesAll = [(bInterfaceAll[ii]+bInterfaceAll[ii+1])/2 for ii in range(2*N)]
    # Vektor zur Definition der Interfaces, die eine Element-Konzentration erfahren: NOADAPT
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
    return(rm,bm,rInterfaceEval,bInterfaceEval,rInterfaceAll,bInterfaceAll,rFacesAll,bFacesAll,rInterfacesBias)


def ellipse_arc_points(a, b, phi_end, num_pts):
    phis = np.linspace(0.0, phi_end, num_pts)
    x_pts = [ct.pol2cart_x(a, phi) for phi in phis]
    y_pts = [ct.pol2cart_y(b, phi) for phi in phis]
    return list(zip(x_pts, y_pts))

rk =[r0]
bk =[b0]
for ii in range (N):
    rk.append (rk[ii]+dL)
    bk.append (bk[ii]+dL)

rN = rk[-1]
bN = bk[-1]

rm,bm,rInterfaceEval,bInterfaceEval,rInterfaceAll,bInterfaceAll,rFacesAll,bFacesAll,rInterfacesBias = radialGeometryParameters(N,rk,bk,iInterfaceEval)

# Define thetaEdges globally for use in partitioning, with bias
#thetaEdges = np.array([thetaEnd * (float(i) / nCirc) ** theta_bias for i in range(nCirc+1)])
# Define thetaEdges_mid as the midpoints between thetaEdges
#thetaEdges_mid = np.array([(thetaEdges[i] + thetaEdges[i+1]) / 2.0 for i in range(len(thetaEdges)-1)])
# Define yOffsets (axial/height) globally for use in partitioning, with bias
yOffsets = np.array([0.0 - (0.0 - (-height)) * (1 - float(i) / nAxial) ** axial_bias for i in range(nAxial+1)])
# Define yOffsets_mid as the midpoints between yOffsets
yOffsets_mid = np.array([(yOffsets[i] + yOffsets[i+1]) / 2.0 for i in range(len(yOffsets)-1)])
# Define meridional_Offsets (shell/radial) globally for use in partitioning, with bias
meridional_Offsets = np.array([phiOpening * (float(i) / nMeridional) ** meridional_bias for i in range(nMeridional+1)])
#np.array([0.0 + (bk[-1] - 0.0) * (float(i) / nMeridional) ** meridional_bias for i in range(nMeridional+1)])
# Define meridional_Offsets_mid as the midpoints between meridional_Offsets
meridional_Offsets_mid = np.array([(meridional_Offsets[i] + meridional_Offsets[i+1]) / 2.0 for i in range(len(meridional_Offsets)-1)])
# Add the mid value between last meridionalOffset and phiOpening
#additional_mid_value = (meridional_Offsets[-1] + phiOpening) / 2.0
#meridional_Offsets_mid = np.append(meridional_Offsets_mid, additional_mid_value)


# Add these parameters near the geometry parameters section:

# Replace the existing thetaEdges definition with:
thetaEdges = np.array([thetaStart + thetaOpening * (float(i) / nCirc) for i in range(nCirc+1)])
thetaEdges_mid = np.array([(thetaEdges[i] + thetaEdges[i+1]) / 2.0 for i in range(len(thetaEdges)-1)])


def create_shells_and_flanges(compositeModel, N, rk, bk, phiEdges, height):
    # Open a log file to write debugging information
    log_file_path = os.path.join(os.getcwd(), "debug_log.txt")
    with open(log_file_path, "w") as log_file:
        shells = []
        for k in range(N):
            log_file.write("Creating shell for layer: {}\n".format(k))
            log_file.write("rk[k]: {}, rk[k+1]: {}\n".format(rk[k], rk[k+1]))
            log_file.write("bk[k]: {}, bk[k+1]: {}\n".format(bk[k], bk[k+1]))
            # --- Create Rotated Datum Plane ---
            datum_axis = compositeModel.rootAssembly.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
            log_file.write("Created datum axis along Y-axis.\n")
            yz_plane = compositeModel.rootAssembly.DatumPlaneByPrincipalPlane(
                principalPlane=XYPLANE, offset=0.0
            )
            log_file.write("Created base XY plane.\n")
            rotated_plane = compositeModel.rootAssembly.DatumPlaneByRotation(
                plane=compositeModel.rootAssembly.datums[yz_plane.id],
                axis=compositeModel.rootAssembly.datums[datum_axis.id],
                angle=math.degrees(thetaStart)
            )
            log_file.write("Created rotated datum plane for layer: {}\n".format(k))
            # --- Shell ---
            arc_pts_inner = ellipse_arc_points(rk[k], bk[k], phiEdges[-1], 1000)
            arc_pts_outer = ellipse_arc_points(rk[k+1], bk[k+1], phiEdges[-1], 1000)
            log_file.write("Inner arc points (first 5): {}\n".format(arc_pts_inner[:5]))
            log_file.write("Outer arc points (first 5): {}\n".format(arc_pts_outer[:5]))
            # Create the shell sketch on the rotated plane
            shell_sketch = compositeModel.ConstrainedSketch(
                name='Shell_Sketch_%d' % k,
                sheetSize=2 * rk[-1],
                transform=compositeModel.rootAssembly.MakeSketchTransform(
                    sketchPlane=compositeModel.rootAssembly.datums[rotated_plane.id],
                    sketchPlaneSide=SIDE1,
                    origin=(0.0, 0.0, 0.0)
                )
            )
            log_file.write("Created shell sketch for layer: {}\n".format(k))
            # Add construction line for rotation axis
            shell_sketch.ConstructionLine(point1=(0.0, 0.0), point2=(0.0, height))
            # Create the base geometry
            shell_sketch.Spline(points=arc_pts_inner)
            shell_sketch.Spline(points=arc_pts_outer)
            shell_sketch.Line(point1=(rk[k+1], 0.0), point2=(rk[k+1], -height))
            shell_sketch.Line(point1=(rk[k+1], -height), point2=(rk[k], -height))
            shell_sketch.Line(point1=(rk[k], -height), point2=(rk[k], 0.0))
            shell_sketch.Line(point1=(ct.pol2cart_x(bk[k], phiEdges[-1]), ct.pol2cart_y(bk[k], phiEdges[-1])),
                              point2=(ct.pol2cart_x(bk[k+1], phiEdges[-1]), ct.pol2cart_y(bk[k+1], phiEdges[-1])))
            log_file.write("Shell sketch geometry keys: {}\n".format(list(shell_sketch.geometry.keys())))
            # Create the shell part by revolving the sketch
            shell_part = compositeModel.Part(name='Shell_%d' % k, dimensionality=THREE_D, type=DEFORMABLE_BODY)
            shell_part.BaseSolidRevolve(
                sketch=shell_sketch,
                angle=math.degrees(thetaOpening),
                pitch=0.0
            )
            log_file.write("Created shell part for layer: {}\n".format(k))
            shells.append(shell_part)
    return shells, []

def merge_shells_and_flanges(compositeModel, shells, flanges):
    compositeAssembly = compositeModel.rootAssembly
    instances = []
    for part in shells + flanges:
        inst = compositeAssembly.Instance(name='I_%s' % part.name, part=part, dependent=ON)
        instances.append(inst)
    compositeAssembly.InstanceFromBooleanMerge(name='Cylindrically_Curved_Composite',
                                              instances=tuple(instances),
                                              keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)
    merged_part = compositeModel.parts['Cylindrically_Curved_Composite']
    return merged_part

def _create_individual_sections(compositeModel, plyAngle, dL):
    for k, angle in enumerate(plyAngle):
        shell_section_name = 'Section_Shell_%d' % k
        compositeModel.HomogeneousSolidSection(
            name=shell_section_name,
            material='CFRP',
            thickness=None
        )
        
        flange_section_name = 'Section_Flange_%d' % k
        compositeModel.HomogeneousSolidSection(
            name=flange_section_name,
            material='CFRP',
            thickness=None
        )
    
    compositeModel.HomogeneousSolidSection(
        name='General_CFRP_Section',
        material='CFRP',
        thickness=None
    )

shells, flanges = create_shells_and_flanges(compositeModel, N, rk, bk, phiEdges, height)
print("=== GEOMETRY CREATION DEBUG ===")
print("Created %d shells and %d flanges" % (len(shells), len(flanges)))

curvedCompositePart = merge_shells_and_flanges(compositeModel, shells, flanges)
print("=== MERGED GEOMETRY DEBUG ===")
print("Merged part name: %s" % curvedCompositePart.name)
print("Number of cells: %d" % len(curvedCompositePart.cells))
print("Number of faces: %d" % len(curvedCompositePart.faces))
print("Number of edges: %d" % len(curvedCompositePart.edges))
print("Number of vertices: %d" % len(curvedCompositePart.vertices))

# Check cell validity
for i, cell in enumerate(curvedCompositePart.cells[:3]):  # Check first 3 cells
    try:
        cell_point = cell.pointOn[0]
        print("Cell %d point: x=%.3f, y=%.3f, z=%.3f" % (i, cell_point[0], cell_point[1], cell_point[2]))
    except Exception as e:
        print("Cell %d error: %s" % (i, str(e)))

compositeAssembly = compositeModel.rootAssembly
curvedCompositeInstance = compositeAssembly.instances['Cylindrically_Curved_Composite-1']

def coordinateSystemCylinder():
    curvedCompositePart.DatumCsysByThreePoints(coordSysType = CYLINDRICAL, name='Cylindrical_CoordinateSystem',
                                        point1 = (rk[-1], 0, 0), 
                                        point2 = (0, 0, rk[-1]),
                                        origin= (0, 0, 0))
    curvedCompositeCSCylId = curvedCompositePart.features['Cylindrical_CoordinateSystem'].id
    curvedCompositeCylCoordSys = compositeModel.rootAssembly.instances['Cylindrically_Curved_Composite-1'].datums[curvedCompositeCSCylId]
    return curvedCompositeCylCoordSys

curvedCompositeCylCoordSys = coordinateSystemCylinder()

def partition_curved_composite_axial_circ(curvedCompositePart):
    # Axial partitions remain unchanged
    for y_offset in yOffsets[1:]:
        datum_plane = curvedCompositePart.DatumPlaneByPrincipalPlane(
            principalPlane=XZPLANE, 
            offset=y_offset
        )
        curvedCompositePart.PartitionCellByDatumPlane(
            datumPlane=curvedCompositePart.datums[datum_plane.id], 
            cells=curvedCompositePart.cells
        )
    
    # Modified circumferential partitions to account for thetaStart
    for angle in thetaEdges[1:-1]:
        datum_axis = curvedCompositePart.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
        yz_plane = curvedCompositePart.DatumPlaneByPrincipalPlane(
            principalPlane=XYPLANE, 
            offset=0.0
        )
        # Add thetaStart to the rotation angle
        total_angle = np.degrees(np.pi/2-angle)
        datum_plane = curvedCompositePart.DatumPlaneByRotation(
            plane=curvedCompositePart.datums[yz_plane.id],
            axis=curvedCompositePart.datums[datum_axis.id], 
            angle=total_angle
        )
        curvedCompositePart.PartitionCellByDatumPlane(
            datumPlane=curvedCompositePart.datums[datum_plane.id], 
            cells=curvedCompositePart.cells
        )
    
    # Meridional partitions
    for angle in -meridional_Offsets[1:-1]:
        datum_axis = curvedCompositePart.DatumAxisByPrincipalAxis(principalAxis=XAXIS)
        xz_plane = curvedCompositePart.DatumPlaneByPrincipalPlane(
            principalPlane=XZPLANE, 
            offset=0.0
        )
        # Add thetaStart to account for the rotated initial position
        rotated_angle = np.degrees(angle)
        datum_plane = curvedCompositePart.DatumPlaneByRotation(
            plane=curvedCompositePart.datums[xz_plane.id],
            axis=curvedCompositePart.datums[datum_axis.id], 
            angle=rotated_angle
        )
        curvedCompositePart.PartitionCellByDatumPlane(
            datumPlane=curvedCompositePart.datums[datum_plane.id], 
            cells=curvedCompositePart.cells
        )

partition_curved_composite_axial_circ(curvedCompositePart)


def mesh_curved_composite(curvedCompositePart, mRadial, mAxial, mMeridional, mCirc):
    radial_midpoints = []
    for r, b in zip(rm, bm):
        for theta in thetaEdges:
           for y in yOffsets:
               x = ct.cyl2cart_x(r, np.pi/2+theta)
               z = ct.cyl2cart_z(r, np.pi/2+theta)
               radial_midpoints.append((x, y, z))
        for t in thetaEdges:
           for m in meridional_Offsets[1:-1]:
                theta = np.pi/2+t
                x = ct.pol2cart3D_x_intersect(r, b, m, theta)
                y = ct.pol2cart3D_y_intersect(r, b, m, theta)
                z = ct.pol2cart3D_z_intersect(r, b, m, theta)
                radial_midpoints.append((x, y, z))
        if phiOpening == np.pi/2:
            radial_midpoints.append((0, b, 0))
    circ_midpoints = []
    for r, b in zip(rk, bk):
        for theta in thetaEdges_mid:
            for y in yOffsets:
                x = ct.cyl2cart_x(r, np.pi/2+theta)
                z = ct.cyl2cart_z(r, np.pi/2+theta)
                circ_midpoints.append((x, y, z))
            for m in meridional_Offsets[1:-1]:
                x = ct.pol2cart3D_x_intersect(r, b, m, np.pi/2+theta)
                y = ct.pol2cart3D_y_intersect(r, b, m, np.pi/2+theta)
                z = ct.pol2cart3D_z_intersect(r, b, m, np.pi/2+theta)
                circ_midpoints.append((x, y, z))
    
    axial_midpoints = []
    for r, b in zip(rk, bk):
        for theta in thetaEdges:
            for y in yOffsets_mid:
                x = ct.cyl2cart_x(r, np.pi/2+theta)
                z = ct.cyl2cart_z(r, np.pi/2+theta)
                axial_midpoints.append((x, y, z))
    
    print("=== EDGE SEEDING ===")
    radial_edges_found = []
    circ_edges_found = []
    axial_edges_found = []
    
    print("Seeding radial edges...")
    for pt in radial_midpoints:
        edge_found = curvedCompositePart.edges.findAt((pt,))
        radial_edges_found.extend(edge_found)
    if radial_edges_found:
        try:
            curvedCompositePart.seedEdgeByNumber(edges=radial_edges_found, number=mRadial, constraint=FINER)
            print("✓ Seeded %d radial edges with %d elements" % (len(radial_edges_found), mRadial))
        except Exception as e:
            print("✗ Error seeding radial edges: %s" % str(e))
    else:
        print("⚠ No radial edges found!")
    
    print("Seeding circumferential edges...")
    for pt in circ_midpoints:
        edge_found = curvedCompositePart.edges.findAt((pt,))
        circ_edges_found.extend(edge_found)
    if circ_edges_found:
        try:
            curvedCompositePart.seedEdgeByNumber(edges=circ_edges_found, number=mCirc, constraint=FINER)
            print("✓ Seeded %d circumferential edges with %d elements" % (len(circ_edges_found), mCirc))
        except Exception as e:
            print("✗ Error seeding circumferential edges: %s" % str(e))
    else:
        print("⚠ No circumferential edges found!")
    
    print("Seeding axial edges...")
    for pt in axial_midpoints:
        edge_found = curvedCompositePart.edges.findAt((pt,))
        axial_edges_found.extend(edge_found)
    if axial_edges_found:
        try:
            curvedCompositePart.seedEdgeByNumber(edges=axial_edges_found, number=mAxial, constraint=FINER)
            print("✓ Seeded %d axial edges with %d elements" % (len(axial_edges_found), mAxial))
        except Exception as e:
            print("✗ Error seeding axial edges: %s" % str(e))
    else:
        print("⚠ No axial edges found!")
    
    # Seed meridional edges
    print("Seeding meridional edges...")
    meridional_edges_found = []
    meridional_midpoints = []
    for r, b in zip(rk, bk):
        for theta in thetaEdges:
            for m in meridional_Offsets_mid:
                x = ct.pol2cart3D_x_intersect(r, b, m, np.pi/2+theta)
                y = ct.pol2cart3D_y_intersect(r, b, m, np.pi/2+theta)
                z = ct.pol2cart3D_z_intersect(r, b, m, np.pi/2+theta)
                meridional_midpoints.append((x, y, z))
    
    for pt in meridional_midpoints:
        edge_found = curvedCompositePart.edges.findAt((pt,))
        meridional_edges_found.extend(edge_found)
    if meridional_edges_found:
        try:
            curvedCompositePart.seedEdgeByNumber(edges=meridional_edges_found, number=mMeridional, constraint=FINER)
            print("✓ Seeded %d meridional edges with %d elements" % (len(meridional_edges_found), mMeridional))
        except Exception as e:
            print("✗ Error seeding meridional edges: %s" % str(e))
    else:
        print("⚠ No meridional edges found!")
    
    # Store edge collections for potential fallback approaches
    edge_collections = {
        'radial': radial_edges_found,
        'circumferential': circ_edges_found, 
        'axial': axial_edges_found,
        'meridional': meridional_edges_found
    }
    
    print("=== MESH GENERATION ===")
    print("Setting mesh controls...")
    try:
        curvedCompositePart.setMeshControls(regions=curvedCompositePart.cells, technique=SWEEP)
        print("✓ Mesh controls set successfully (SWEEP technique)")
    except Exception as e:
        print("✗ Error setting SWEEP mesh controls: %s" % str(e))
        print("Trying alternative mesh technique...")
        try:
            # Try structured meshing as alternative
            curvedCompositePart.setMeshControls(regions=curvedCompositePart.cells, technique=STRUCTURED)
            print("✓ Mesh controls set successfully (STRUCTURED technique)")
        except Exception as e2:
            print("✗ Error setting STRUCTURED mesh controls: %s" % str(e2))
            print("Trying free meshing as fallback...")
            try:
                # Try free meshing as last resort
                curvedCompositePart.setMeshControls(regions=curvedCompositePart.cells, technique=FREE)
                print("✓ Mesh controls set successfully (FREE technique)")
            except Exception as e3:
                print("✗ Error setting FREE mesh controls: %s" % str(e3))
        
    print("Setting element type...")
    try:
        curvedCompositePart.setElementType(regions=(curvedCompositePart.cells,),
            elemTypes=(mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD),))
        print("✓ Element type set successfully")
    except Exception as e:
        print("✗ Error setting element type: %s" % str(e))
        
    print("Generating mesh...")
    
    # First, let's check the geometry is valid for meshing
    print("Pre-mesh validation...")
    try:
        all_cells = curvedCompositePart.cells
        print("Number of cells to mesh: %d" % len(all_cells))
        
        # Check cell geometry validity
        valid_cells = []
        invalid_cells = []
        for i, cell in enumerate(all_cells):
            try:
                cell_point = cell.pointOn[0]
                valid_cells.append(cell)
                if i < 5:  # Show first 5 cells
                    print("Valid cell %d point: x=%.3f, y=%.3f, z=%.3f" % (i, cell_point[0], cell_point[1], cell_point[2]))
            except Exception as cell_e:
                print("Invalid cell %d: %s" % (i, str(cell_e)))
                invalid_cells.append(cell)
        
        print("Valid cells: %d, Invalid cells: %d" % (len(valid_cells), len(invalid_cells)))
        
        if len(invalid_cells) > 0:
            print("⚠ Some cells are invalid - attempting to mesh only valid cells")
            mesh_cells = tuple(valid_cells)
        else:
            mesh_cells = curvedCompositePart.cells
            
    except Exception as validation_e:
        print("Error during pre-mesh validation: %s" % str(validation_e))
        mesh_cells = curvedCompositePart.cells
    
    # Try multiple mesh generation approaches
    mesh_success = False
    
    # Approach 1: Standard sweep mesh with all cells
    if not mesh_success:
        try:
            print("Attempting mesh generation (Approach 1: Standard SWEEP)...")
            curvedCompositePart.generateMesh()
            
            # Immediate validation
            nodes = curvedCompositePart.nodes
            elements = curvedCompositePart.elements
            print("Number of nodes: %d" % len(nodes))
            print("Number of elements: %d" % len(elements))
            
            if len(nodes) > 0 and len(elements) > 0:
                print("✓ Mesh generated successfully with Approach 1")
                mesh_success = True
            else:
                print("✗ Approach 1 failed - no nodes/elements created")
                # Clear any partial mesh
                try:
                    curvedCompositePart.deleteMesh()
                except:
                    pass
                
        except Exception as e1:
            print("✗ Approach 1 failed: %s" % str(e1))
            try:
                curvedCompositePart.deleteMesh()
            except:
                pass
    
    # Approach 2: Try with coarser seeding if first approach fails
    # if not mesh_success:
    #     try:
    #         print("Attempting mesh generation (Approach 2: Coarser seeding)...")
            
    #         # Reduce mesh density
    #         coarse_radial = max(1, mRadial // 2)
    #         coarse_circ = max(1, mCirc // 2)
    #         coarse_axial = max(1, mAxial // 2)
    #         coarse_meridional = max(1, mMeridional // 2)
            
    #         print("Using coarser seeding: radial=%d, circ=%d, axial=%d, meridional=%d" % (coarse_radial, coarse_circ, coarse_axial, coarse_meridional))
            
    #         # Re-seed with coarser mesh
    #         if edge_collections['radial']:
    #             curvedCompositePart.seedEdgeByNumber(edges=edge_collections['radial'], number=coarse_radial, constraint=FINER)
    #         if edge_collections['circumferential']:
    #             curvedCompositePart.seedEdgeByNumber(edges=edge_collections['circumferential'], number=coarse_circ, constraint=FINER)
    #         if edge_collections['axial']:
    #             curvedCompositePart.seedEdgeByNumber(edges=edge_collections['axial'], number=coarse_axial, constraint=FINER)
    #         if edge_collections['meridional']:
    #             curvedCompositePart.seedEdgeByNumber(edges=edge_collections['meridional'], number=coarse_meridional, constraint=FINER)
            
    #         curvedCompositePart.generateMesh()
            
    #         nodes = curvedCompositePart.nodes
    #         elements = curvedCompositePart.elements
    #         print("Number of nodes: %d" % len(nodes))
    #         print("Number of elements: %d" % len(elements))
            
    #         if len(nodes) > 0 and len(elements) > 0:
    #             print("✓ Mesh generated successfully with Approach 2 (coarser mesh)")
    #             mesh_success = True
    #         else:
    #             print("✗ Approach 2 failed - no nodes/elements created")
    #             try:
    #                 curvedCompositePart.deleteMesh()
    #             except:
    #                 pass
                
    #     except Exception as e2:
    #         print("✗ Approach 2 failed: %s" % str(e2))
    #         try:
    #             curvedCompositePart.deleteMesh()
    #         except:
    #             pass
    
    # Approach 3: Try structured mesh as fallback
    if not mesh_success:
        try:
            print("Attempting mesh generation (Approach 3: STRUCTURED technique)...")
            
            # Change mesh controls to structured
            curvedCompositePart.setMeshControls(regions=curvedCompositePart.cells, technique=STRUCTURED)
            curvedCompositePart.generateMesh()
            
            nodes = curvedCompositePart.nodes
            elements = curvedCompositePart.elements
            print("Number of nodes: %d" % len(nodes))
            print("Number of elements: %d" % len(elements))
            
            if len(nodes) > 0 and len(elements) > 0:
                print("✓ Mesh generated successfully with Approach 3 (STRUCTURED)")
                mesh_success = True
            else:
                print("✗ Approach 3 failed - no nodes/elements created")
                try:
                    curvedCompositePart.deleteMesh()
                except:
                    pass
                
        except Exception as e3:
            print("✗ Approach 3 failed: %s" % str(e3))
            try:
                curvedCompositePart.deleteMesh()
            except:
                pass
    
    # Approach 4: Try free mesh as last resort
    if not mesh_success:
        try:
            print("Attempting mesh generation (Approach 4: FREE technique)...")
            
            # Change mesh controls to free
            curvedCompositePart.setMeshControls(regions=curvedCompositePart.cells, technique=FREE)
            curvedCompositePart.generateMesh()
            
            nodes = curvedCompositePart.nodes
            elements = curvedCompositePart.elements
            print("Number of nodes: %d" % len(nodes))
            print("Number of elements: %d" % len(elements))
            
            if len(nodes) > 0 and len(elements) > 0:
                print("✓ Mesh generated successfully with Approach 4 (FREE)")
                mesh_success = True
            else:
                print("✗ Approach 4 failed - no nodes/elements created")
                
        except Exception as e4:
            print("✗ Approach 4 failed: %s" % str(e4))
    
    # Final mesh status report
    if mesh_success:
        print("=== FINAL MESH STATUS ===")
        try:
            final_nodes = curvedCompositePart.nodes
            final_elements = curvedCompositePart.elements
            print("✓ MESH GENERATION SUCCESSFUL")
            print("Final node count: %d" % len(final_nodes))
            print("Final element count: %d" % len(final_elements))
            
            # Additional mesh quality checks
            if len(final_elements) > 0:
                print("Mesh density per cell: %.1f elements" % (float(len(final_elements)) / len(all_cells)))
                
        except Exception as final_e:
            print("Error in final mesh check: %s" % str(final_e))
    else:
        print("=== MESH GENERATION FAILED ===")
        print("✗ All mesh generation approaches failed")
        print("This may be due to:")
        print("  - Complex 4-layer geometry causing mesh conflicts")
        print("  - Edge seeding incompatibilities") 
        print("  - Invalid cell geometries from boolean operations")
        print("  - Insufficient geometric tolerance")
        
        # Diagnostic information
        try:
            print("=== DIAGNOSTIC INFORMATION ===")
            print("Total edges: %d" % len(curvedCompositePart.edges))
            print("Total faces: %d" % len(curvedCompositePart.faces))
            print("Total cells: %d" % len(curvedCompositePart.cells))
            print("Total vertices: %d" % len(curvedCompositePart.vertices))
        except:
            print("Could not get diagnostic information")
    
    # Create edge sets AFTER mesh generation using Assembly-level approach
    # if mesh_success:
    #     print("=== CREATING EDGE SETS USING ASSEMBLY APPROACH ===")
    #     try:
    #         # Get the assembly and instance for set creation
    #         compositeAssembly = compositeModel.rootAssembly
    #         curvedCompositeInstance = compositeAssembly.instances['Cylindrically_Curved_Composite-1']
            
    #         # 1. Create RADIAL edge sets
    #         print("Creating radial edge sets...")
    #         curvedCompositeStringEdges = ''
    #         for r, b in zip(rm, bm):
    #             for theta in thetaEdges:
    #                 for y in yOffsets:
    #                     x = ct.cyl2cart_x(r, np.pi/2+theta)
    #                     z = ct.cyl2cart_z(r, np.pi/2+theta)
    #                     curvedCompositeStringEdges += str(((x, y, z),),) + ','
    #                 for m in meridional_Offsets[1:-1]:
    #                     x = ct.pol2cart3D_x_intersect(r, b, m, np.pi/2+theta)
    #                     y = ct.pol2cart3D_y_intersect(r, b, m, np.pi/2+theta)
    #                     z = ct.pol2cart3D_z_intersect(r, b, m, np.pi/2+theta)
    #                     curvedCompositeStringEdges += str(((x, y, z),),) + ','
    #             if phiOpening == np.pi/2:
    #                 curvedCompositeStringEdges += str(((0, b, 0),),) + ','
            
    #         if curvedCompositeStringEdges.endswith(','):
    #             curvedCompositeStringEdges = curvedCompositeStringEdges[:-1]
            
    #         curvedCompositeStringEdgesExec = ('radialEdgesForSet = ' +
    #                                         'curvedCompositeInstance.edges.findAt(' + 
    #                                         curvedCompositeStringEdges + ')')
    #         local_vars = {'curvedCompositeInstance': curvedCompositeInstance}
    #         exec(curvedCompositeStringEdgesExec, globals(), local_vars)
    #         radialEdgesForSet = local_vars['radialEdgesForSet']
            
    #         if radialEdgesForSet:
    #             compositeAssembly.Set(edges=radialEdgesForSet, name='RadialEdges')
    #             print("✓ Created RadialEdges set with %d edges" % len(radialEdgesForSet))
            
    #         # 2. Create CIRCUMFERENTIAL edge sets
    #         print("Creating circumferential edge sets...")
    #         curvedCompositeStringEdges = ''
    #         for r, b in zip(rk, bk):
    #             for theta in thetaEdges_mid:
    #                 for y in yOffsets:
    #                     x = ct.cyl2cart_x(r, np.pi/2+theta)
    #                     z = ct.cyl2cart_z(r, np.pi/2+theta)
    #                     curvedCompositeStringEdges += str(((x, y, z),),) + ','
    #                 for m in meridional_Offsets[1:-1]:
    #                     x = ct.pol2cart3D_x_intersect(r, b, m, np.pi/2+theta)
    #                     y = ct.pol2cart3D_y_intersect(r, b, m, np.pi/2+theta)
    #                     z = ct.pol2cart3D_z_intersect(r, b, m, np.pi/2+theta)
    #                     curvedCompositeStringEdges += str(((x, y, z),),) + ','
            
    #         if curvedCompositeStringEdges.endswith(','):
    #             curvedCompositeStringEdges = curvedCompositeStringEdges[:-1]
            
    #         curvedCompositeStringEdgesExec = ('circumferentialEdgesForSet = ' +
    #                                         'curvedCompositeInstance.edges.findAt(' + 
    #                                         curvedCompositeStringEdges + ')')
    #         local_vars = {'curvedCompositeInstance': curvedCompositeInstance}
    #         exec(curvedCompositeStringEdgesExec, globals(), local_vars)
    #         circumferentialEdgesForSet = local_vars['circumferentialEdgesForSet']
            
    #         if circumferentialEdgesForSet:
    #             compositeAssembly.Set(edges=circumferentialEdgesForSet, name='CircumferentialEdges')
    #             print("✓ Created CircumferentialEdges set with %d edges" % len(circumferentialEdgesForSet))
            
    #         # 3. Create AXIAL edge sets
    #         print("Creating axial edge sets...")
    #         curvedCompositeStringEdges = ''
    #         for r, b in zip(rk, bk):
    #             for theta in thetaEdges:
    #                 for y in yOffsets_mid:
    #                     x = ct.cyl2cart_x(r, np.pi/2+theta)
    #                     z = ct.cyl2cart_z(r, np.pi/2+theta)
    #                     curvedCompositeStringEdges += str(((x, y, z),),) + ','
            
    #         if curvedCompositeStringEdges.endswith(','):
    #             curvedCompositeStringEdges = curvedCompositeStringEdges[:-1]
            
    #         curvedCompositeStringEdgesExec = ('axialEdgesForSet = ' +
    #                                         'curvedCompositeInstance.edges.findAt(' + 
    #                                         curvedCompositeStringEdges + ')')
    #         local_vars = {'curvedCompositeInstance': curvedCompositeInstance}
    #         exec(curvedCompositeStringEdgesExec, globals(), local_vars)
    #         axialEdgesForSet = local_vars['axialEdgesForSet']
            
    #         if axialEdgesForSet:
    #             compositeAssembly.Set(edges=axialEdgesForSet, name='AxialEdges')
    #             print("✓ Created AxialEdges set with %d edges" % len(axialEdgesForSet))
            
    #         # 4. Create MERIDIONAL edge sets  
    #         print("Creating meridional edge sets...")
    #         curvedCompositeStringEdges = ''
    #         for r, b in zip(rk, bk):
    #             for theta in thetaEdges:
    #                 for m in meridional_Offsets_mid:
    #                     x = ct.pol2cart3D_x_intersect(r, b, m, np.pi/2+theta)
    #                     y = ct.pol2cart3D_y_intersect(r, b, m, np.pi/2+theta)
    #                     z = ct.pol2cart3D_z_intersect(r, b, m, np.pi/2+theta)
    #                     curvedCompositeStringEdges += str(((x, y, z),),) + ','
            
    #         if curvedCompositeStringEdges.endswith(','):
    #             curvedCompositeStringEdges = curvedCompositeStringEdges[:-1]
            
    #         curvedCompositeStringEdgesExec = ('meridionalEdgesForSet = ' +
    #                                         'curvedCompositeInstance.edges.findAt(' + 
    #                                         curvedCompositeStringEdges + ')')
    #         local_vars = {'curvedCompositeInstance': curvedCompositeInstance}
    #         exec(curvedCompositeStringEdgesExec, globals(), local_vars)
    #         meridionalEdgesForSet = local_vars['meridionalEdgesForSet']
            
    #         if meridionalEdgesForSet:
    #             compositeAssembly.Set(edges=meridionalEdgesForSet, name='MeridionalEdges')
    #             print("✓ Created MeridionalEdges set with %d edges" % len(meridionalEdgesForSet))
            
    #         print("✓ All edge sets created successfully using Assembly approach")
                
    #     except Exception as e:
    #         print("✗ Error creating edge sets using Assembly approach: %s" % str(e))
    # else:
    #     print("⚠ Skipping edge set creation due to mesh generation failure")

mesh_curved_composite(curvedCompositePart, mRadial, mAxial, mMeridional, mCirc)

def define_materials_sections_merged(compositeModel, plyAngle, dL, E1, E2, E3, Nu12, Nu13, Nu23, G12, G13, G23):
    compositeModel.Material(name='CFRP')
    compositeModel.materials['CFRP'].Elastic(type=ENGINEERING_CONSTANTS,
                                table=((E1, E2, E3,
                                        Nu12, Nu13, Nu23,
                                        G12, G13, G23),))
    
    merged_part = compositeModel.parts['Cylindrically_Curved_Composite']
    total_cells = merged_part.cells
    num_layers = len(plyAngle)
    
    _create_individual_sections(compositeModel, plyAngle, dL)
    
    layup = merged_part.CompositeLayup(
        name='Layup_Merged',
        description='Composite layup after Boolean merge',
        elementType=SOLID,
        symmetric=False,
        thicknessAssignment=FROM_SECTION
    )
    
    shell_cells_by_layer = [[] for _ in range(num_layers)]
    flange_cells_by_layer = [[] for _ in range(num_layers)]
    
    shell_cells = []
    flange_cells = []
    
    for cell in total_cells:
        cell_point = cell.pointOn[0]
        if cell_point[1] >= 0:  # +Y region = shell
            shell_cells.append(cell)
        else:  # -Y region = flange
            flange_cells.append(cell)
    
    for k in range(num_layers):
        r_layer = rm[k]  # Middle radius for this layer
        b_layer = bm[k]  # Middle b-parameter for this layer
        
        # Calculate edge coordinates for this specific layer using the same method as meshing
        layer_edge_coords = []
        for theta in thetaEdges[1:-1]:
            for m in meridional_Offsets[1:-1]:
                x = ct.pol2cart3D_x_intersect(r_layer, b_layer, m, np.pi/2+theta)
                y = ct.pol2cart3D_y_intersect(r_layer, b_layer, m, np.pi/2+theta)
                z = ct.pol2cart3D_z_intersect(r_layer, b_layer, m, np.pi/2+theta)
                layer_edge_coords.append((x, y, z))
        
        # Find edges at these coordinates
        found_edges = []
        for coord in layer_edge_coords:
            # Use findAt with a small tolerance to find edges near these coordinates
            edges_at_coord = merged_part.edges.findAt((coord,), )
            if edges_at_coord:
                found_edges.extend(edges_at_coord)
        
        # Find all shell cells that contain any of these layer-specific edges
        layer_shell_cells = []
        for cell in shell_cells:
            # Check if this cell contains any of the layer-specific edges
            cell_edges = cell.getEdges()
            for edge_index in cell_edges:
                edge = merged_part.edges[edge_index]
                if edge in found_edges:
                    if cell not in layer_shell_cells:
                        layer_shell_cells.append(cell)
                    break
        
        shell_cells_by_layer[k] = layer_shell_cells
    
    flange_cell_distances = []
    for cell in flange_cells:
        cell_point = cell.pointOn[0]
        x, y, z = cell_point
        radial_dist = np.sqrt(x**2 + z**2)  # Distance from Y-axis
        flange_cell_distances.append((cell, radial_dist, cell_point))
    
    # Sort flange cells by radial distance
    flange_cell_distances.sort(key=lambda x: x[1])
    
    # Assign flange cells to layers based on rm values
    for cell, radial_dist, point in flange_cell_distances:
        # Find closest rm
        best_layer = 0
        min_distance = abs(radial_dist - rm[0])
        
        for k in range(num_layers):
            r_mid = rm[k] if k < len(rm) else rm[-1]
            distance_to_layer = abs(radial_dist - r_mid)
            if distance_to_layer < min_distance:
                min_distance = distance_to_layer
                best_layer = k
        
        flange_cells_by_layer[best_layer].append(cell)
    
    # Flatten for overall statistics and compatibility with existing code
    shell_cells = [cell for layer in shell_cells_by_layer for cell in layer]
    flange_cells = [cell for layer in flange_cells_by_layer for cell in layer]
    
    # NOTE: organize_cells_by_layer functions are no longer needed
    # since we now use edge-based identification for direct layer assignment
    
    # Legacy functions kept for reference (now unused):
    # def organize_cells_by_layer(...)
    # def organize_shell_cells_by_edges(...)  
    # def organize_flange_cells_by_radial_distance(...)
    # Layers are already organized by the edge-based identification above
    shell_layers = shell_cells_by_layer
    flange_layers = flange_cells_by_layer
    
    region_all_shells = None
    region_all_flanges = None
    
    all_shell_cells = []
    for layer_cells in shell_layers:
        all_shell_cells.extend(layer_cells)
        
    all_flange_cells = []
    for layer_cells in flange_layers:
        all_flange_cells.extend(layer_cells)
    
    if all_shell_cells:
        shell_cells_list = []
        for i, cell in enumerate(total_cells):
            if cell in all_shell_cells:
                shell_cells_list.append(total_cells[i:i+1])
        
        if shell_cells_list:
            combined_shell_cells = shell_cells_list[0]
            for cell_seq in shell_cells_list[1:]:
                combined_shell_cells = combined_shell_cells + cell_seq
            region_all_shells = regionToolset.Region(cells=combined_shell_cells)
    
    if all_flange_cells:
        flange_cells_list = []
        for i, cell in enumerate(total_cells):
            if cell in all_flange_cells:
                flange_cells_list.append(total_cells[i:i+1])
        
        if flange_cells_list:
            combined_flange_cells = flange_cells_list[0]
            for cell_seq in flange_cells_list[1:]:
                combined_flange_cells = combined_flange_cells + cell_seq
            region_all_flanges = regionToolset.Region(cells=combined_flange_cells)
    
    for k, angle in enumerate(plyAngle):
        layer_shell_cells = shell_layers[k] if k < len(shell_layers) else []
        layer_flange_cells = flange_layers[k] if k < len(flange_layers) else []
        
        layer_all_cells = layer_shell_cells + layer_flange_cells
        
        if not layer_all_cells:
            continue
        
        layer_cells_list = []
        for i, cell in enumerate(total_cells):
            if cell in layer_all_cells:
                layer_cells_list.append(total_cells[i:i+1])
        
        if layer_cells_list:
            combined_layer_cells = layer_cells_list[0]
            for cell_seq in layer_cells_list[1:]:
                combined_layer_cells = combined_layer_cells + cell_seq
            region_layer_cells = regionToolset.Region(cells=combined_layer_cells)
        
            csys_local = merged_part.DatumCsysByThreePoints(
                coordSysType=CYLINDRICAL,
                name='OrientationCSYS_Layer_%d' % k,
                origin=(0.0, 0.0, 0.0),
                point1=(1.0, 0.0, 0.0),
                point2=(0.0, 1.0, 0.0)
            )
            
            if angle == 0.0:
                orientation_type = SPECIFY_ORIENT
                orientation_value = 0.0
            elif angle == 90.0:
                orientation_type = SPECIFY_ORIENT
                orientation_value = 90.0
            else:
                orientation_type = SPECIFY_ORIENT
                orientation_value = angle
            
            layup.CompositePly(
                plyName='Ply_%d' % k,
                region=region_layer_cells,
                material='CFRP',
                thickness=dL,
                orientationType=orientation_type,
                orientationValue=orientation_value,
                axis=AXIS_1,
                additionalRotationType=ROTATION_NONE,
                numIntPoints=3,
                thicknessType=SPECIFY_THICKNESS
            )
    
    merged_part.SectionAssignment(region=region_all_shells, sectionName='Section_Shell_0', offset=0.0)
    merged_part.SectionAssignment(region=region_all_flanges, sectionName='Section_Flange_0', offset=0.0)
    
    csys_global_orientation = merged_part.DatumCsysByThreePoints(
        coordSysType=CYLINDRICAL,
        name='GlobalOrientationCSYS',
        origin=(0.0, 0.0, 0.0),
        point1=(1.0, 0.0, 0.0),
        point2=(0.0, 0.0, 1.0)
    )
    
    for k, angle in enumerate(plyAngle):
        layer_shell_cells = shell_layers[k] if k < len(shell_layers) else []
        layer_flange_cells = flange_layers[k] if k < len(flange_layers) else []
        
        if not layer_shell_cells and not layer_flange_cells:
            continue
        
        if layer_shell_cells:
            shell_cells_list = []
            for i, cell in enumerate(total_cells):
                if cell in layer_shell_cells:
                    shell_cells_list.append(total_cells[i:i+1])
            
            if shell_cells_list:
                combined_shell_layer = shell_cells_list[0]
                for cell_seq in shell_cells_list[1:]:
                    combined_shell_layer = combined_shell_layer + cell_seq
                region_shell_layer = regionToolset.Region(cells=combined_shell_layer)
                
                csys_shell_spherical = merged_part.DatumCsysByThreePoints(
                    coordSysType=SPHERICAL,
                    name='ShellSphericalCSYS_Layer_%d' % k,
                    origin=(0.0, 0.0, 0.0),
                    point1=(1.0, 0.0, 0.0),
                    point2=(0.0, 0.0, 1.0)
                )
                
                merged_part.MaterialOrientation(
                    region=region_shell_layer,
                    orientationType=SYSTEM,
                    localCsys=merged_part.datums[csys_shell_spherical.id],
                    axis=AXIS_1,
                    angle=angle + 90.0,
                    additionalRotationType=ROTATION_ANGLE,
                    stackDirection=STACK_1
                )
        
        if layer_flange_cells:
            flange_cells_list = []
            for i, cell in enumerate(total_cells):
                if cell in layer_flange_cells:
                    flange_cells_list.append(total_cells[i:i+1])
            
            if flange_cells_list:
                combined_flange_layer = flange_cells_list[0]
                for cell_seq in flange_cells_list[1:]:
                    combined_flange_layer = combined_flange_layer + cell_seq
                region_flange_layer = regionToolset.Region(cells=combined_flange_layer)
                
                merged_part.MaterialOrientation(
                    region=region_flange_layer,
                    orientationType=SYSTEM,
                    localCsys=merged_part.datums[csys_global_orientation.id],
                    axis=AXIS_1,
                    angle=angle + 90.0,  # Flange orientation adjustment
                    additionalRotationType=ROTATION_ANGLE,
                    stackDirection=STACK_1
                )

def _create_individual_sections(compositeModel, plyAngle, dL):
    for k, angle in enumerate(plyAngle):
        shell_section_name = 'Section_Shell_%d' % k
        compositeModel.HomogeneousSolidSection(
            name=shell_section_name,
            material='CFRP',
            thickness=None
        )
        
        flange_section_name = 'Section_Flange_%d' % k
        compositeModel.HomogeneousSolidSection(
            name=flange_section_name,
            material='CFRP',
            thickness=None
        )
    
    compositeModel.HomogeneousSolidSection(
        name='General_CFRP_Section',
        material='CFRP',
        thickness=None
    )

define_materials_sections_merged(compositeModel, plyAngle, dL, E1, E2, E3, Nu12, Nu13, Nu23, G12, G13, G23)

def identify_inner_surface(curvedCompositePart, rk, height, thetaEnd):
    """
    Identify inner surface faces at radius rk[0] using systematic coordinate-based approach
    Uses nested loops with theta_mid and yOffset_mid for flange, theta_mid and meridionalOffset_mid for shell
    Also creates a surface set for GUI verification
    """
    inner_faces = []
    
    # Inner flange surface faces - use theta_mid and yOffset_mid
    flange_test_points = []
    for theta in thetaEdges_mid:  # Use circumferential midpoints
        for y in yOffsets_mid:  # Use axial midpoints
            x = rk[0] * np.cos(np.pi/2+theta)
            z = rk[0] * np.sin(np.pi/2+theta)
            flange_test_points.append((x, y, z))
    
    # Inner shell surface faces - use theta_mid and meridionalOffset_mid with pol2cart3D transforms
    shell_test_points = []
    for theta in thetaEdges_mid:  # Use circumferential midpoints
        for m in meridional_Offsets_mid:  # Use meridional midpoints
            x = ct.pol2cart3D_x_intersect(rk[0], bk[0], m, np.pi/2+theta)
            y = ct.pol2cart3D_y_intersect(rk[0], bk[0], m, np.pi/2+theta)
            z = ct.pol2cart3D_z_intersect(rk[0], bk[0], m, np.pi/2+theta)
            shell_test_points.append((x, y, z))
    
    
    # Combine all test points
    all_inner_test_points = flange_test_points + shell_test_points
    
    for i, point in enumerate(all_inner_test_points):
        try:
            found_faces = curvedCompositePart.faces.findAt((point,))
            if found_faces:
                for face in found_faces:
                    if face not in inner_faces:
                        # Validate that this is actually an inner surface face at rk[0]
                        face_point = face.pointOn[0]
                        x, y, z = face_point
                        face_radius = np.sqrt(x**2 + z**2)
                        
                        # Check if it's truly at inner radius with tolerance
                        is_inner_radius = abs(face_radius - rk[0]) < 0.1 * dL
                        # Exclude symmetry planes and end faces
                        is_not_bottom = y > -0.9 * height
                        is_not_z_symmetry = abs(z) > 0.01 * rk[0]
                        is_not_theta_symmetry = abs(abs(np.arctan2(z, x)) - thetaEnd) > 0.05
                        
                        if is_inner_radius and is_not_bottom and is_not_z_symmetry and is_not_theta_symmetry:
                            inner_faces.append(face)
                            surface_type = "flange" if i < len(flange_test_points) else "shell"
        except Exception:
            continue
    
    
    # Create a surface set for GUI verification
    if inner_faces:
        try:
            # Delete existing surface if it exists
            surface_set_name = 'Set_Inner_Faces_Found'
            try:
                if surface_set_name in curvedCompositePart.surfaces.keys():
                    del curvedCompositePart.surfaces[surface_set_name]
            except Exception:
                pass
            
            # Create new surface set for verification
            curvedCompositePart.Surface(name=surface_set_name, side1Faces=tuple(inner_faces))
            
        except Exception as e:
            pass
    return inner_faces

def identify_boundary_surfaces(compositeModel, curvedCompositePart, rk, bk, phiEdges, height, thetaStart, thetaEnd, rm):
    """
    Identify boundary surfaces such as symmetry planes and end faces.
    """
    surfaces = {
        'symmetry_planes': [],
        'end_faces': [],  # Bottom faces (y = -height)
        'inner_faces': []  # Keep for compatibility, but will be empty
    }
    
    all_faces = curvedCompositePart.faces
    
    # Calculate middle theta angle between thetaStart and thetaEnd
    thetamid = (thetaStart + thetaEnd) / 2.0
    
    # Use specific coordinate points like edge identification for bottom faces
    # Test points at rm, y=yOffsets[-1] (which is -height), theta=thetamid
    bottom_face_test_points = []
    for r in rm:  # Use radial midpoints like edge identification
        for theta in thetaEdges_mid:
            x = r * math.cos(math.pi / 2 + theta)
            z = r * math.sin(math.pi / 2 + theta)
            y = -height  # Bottom face y coordinate
            bottom_face_test_points.append((x, y, z))
            print("Bottom face test point: x=%.3f, y=%.3f, z=%.3f (r=%.3f, theta=%.3f)" % (x, y, z, r, theta))
    
    # Find bottom faces using specific coordinate points
    for test_point in bottom_face_test_points:
        try:
            found_faces = curvedCompositePart.faces.findAt((test_point,))
            if found_faces:
                for face in found_faces:
                    if face not in surfaces['end_faces']:
                        surfaces['end_faces'].append(face)
                        print("Identified bottom face at test point: x=%.3f, y=%.3f, z=%.3f" % test_point)
        except Exception as e:
            print("Failed to find face at test point %s: %s" % (test_point, str(e)))
            continue
            
    # Use the robust symmetry plane identification
    surfaces['symmetry_planes'] = get_symmetry_plane_faces(curvedCompositePart, rk, height, thetaStart, thetaEnd)
    
    return surfaces

# Most robust approach: Direct coordinate-based identification
def get_symmetry_plane_faces(curvedCompositePart, rk, height, thetaStart, thetaEnd):
    """
    Find symmetry plane faces using precise coordinate matching.
    Symmetry planes are now defined at thetaStart and thetaEnd.
    """
    symmetry_faces = []
    
    # theta=thetaStart symmetry plane coordinates
    thetaStart_test_points = []
    for r in rm:  # Use all radial midpoints (positive and negative)
        for y in yOffsets_mid:  # Use all axial midpoints
            x = r * math.cos(math.pi / 2 + thetaStart)
            z = r * math.sin(math.pi / 2 + thetaStart)
            thetaStart_test_points.append((x, y, z))
    for r, b in zip(rm, bm):
        for m in meridional_Offsets_mid:
            x = ct.pol2cart3D_x_intersect(r, b, m, np.pi/2+thetaStart)
            y = ct.pol2cart3D_y_intersect(r, b, m, np.pi/2+thetaStart)
            z = ct.pol2cart3D_z_intersect(r, b, m, np.pi/2+thetaStart)
            thetaStart_test_points.append((x, y, z))
    
    # theta=thetaEnd symmetry plane coordinates
    thetaEnd_test_points = []
    for r in rm:  # Use all radial midpoints (positive and negative)
        for y in yOffsets_mid:  # Use all axial midpoints
            x = r * math.cos(math.pi / 2 + thetaEnd)
            z = r * math.sin(math.pi / 2 + thetaEnd)
            thetaEnd_test_points.append((x, y, z))
    for r, b in zip(rm, bm):
        for m in meridional_Offsets_mid:
            x = ct.pol2cart3D_x_intersect(r, b, m, np.pi/2+thetaEnd)
            y = ct.pol2cart3D_y_intersect(r, b, m, np.pi/2+thetaEnd)
            z = ct.pol2cart3D_z_intersect(r, b, m, np.pi/2+thetaEnd)
            thetaEnd_test_points.append((x, y, z))
    
    # Combine all test points
    all_test_points = thetaStart_test_points + thetaEnd_test_points
    
    # Find faces at these specific coordinates
    for i, point in enumerate(all_test_points):
        found_faces = curvedCompositePart.faces.findAt((point,))
        if found_faces:
            for face in found_faces:
                if face not in symmetry_faces:
                    # Validate that this is actually a symmetry plane face
                    face_point = face.pointOn[0]
                    x, y, z = face_point
                    
                    # Check if it's on thetaStart or thetaEnd plane (with tolerance)
                    # The angle calculation needs to account for the π/2 offset in the coordinate system
                    angle = math.atan2(z, x)
                    adjusted_angle = angle - math.pi/2  # Remove the π/2 offset
                    
                    # Normalize angle to [-π, π] range
                    while adjusted_angle > math.pi:
                        adjusted_angle -= 2*math.pi
                    while adjusted_angle < -math.pi:
                        adjusted_angle += 2*math.pi
                        
                    is_thetaStart_plane = abs(adjusted_angle - thetaStart) < 1e-2
                    is_thetaEnd_plane = abs(adjusted_angle - thetaEnd) < 1e-2
                    
                    if is_thetaStart_plane or is_thetaEnd_plane:
                        symmetry_faces.append(face)
    
    return symmetry_faces

def create_face_region_safely(face_list, instance):
    """
    Safely create a region from a list of faces, handling Abaqus GeomSequence requirements
    """
    if not face_list:
        print("create_face_region_safely: No faces provided")
        return None
    
    print("create_face_region_safely: Processing %d faces" % len(face_list))
    
    # Method 1: Use face indices to create a proper sequence
    face_indices = []
    for i, face in enumerate(face_list):
        for j, inst_face in enumerate(instance.faces):
            if inst_face.index == face.index:
                face_indices.append(j)
                print("  Face %d: Part index %d -> Instance index %d" % (i, face.index, j))
                break
        else:
            print("  Face %d: Part index %d -> NOT FOUND in instance!" % (i, face.index))
    
    print("Found %d face indices in instance: %s" % (len(face_indices), face_indices))
    
    if face_indices:
        # Create a sequence of faces from the instance
        selected_faces = []
        for idx in face_indices:
            selected_faces.append(instance.faces[idx:idx+1])
        
        # Combine all face sequences
        if len(selected_faces) == 1:
            combined_faces = selected_faces[0]
        else:
            combined_faces = selected_faces[0]
            for face_seq in selected_faces[1:]:
                combined_faces = combined_faces + face_seq
        
        print("Successfully created combined face sequence with %d faces" % len(combined_faces))
        return regionToolset.Region(faces=combined_faces)
    
    print("No valid face indices found!")
    return None

def create_boundary_conditions(compositeModel, curvedCompositeInstance, surfaces, curvedCompositeCylCoordSys):
    """
    STEP 1: Apply boundary conditions to symmetry planes and bottom faces.
    These are the FACE symmetry conditions applied first.
    Symmetry planes are now at thetaStart and thetaEnd.
    """
    if 'Static_Step' not in compositeModel.steps.keys():
        compositeModel.StaticStep(name='Static_Step', previous='Initial', 
                                  description='Static analysis step')
    # STEP 1: Organize faces into separate groups for boundary conditions
    bottom_faces_no_symmetry = []
    # Get all bottom faces first
    all_bottom_faces = []
    if 'end_faces' in surfaces and surfaces['end_faces']:
        print("Processing %d end_faces for bottom face identification..." % len(surfaces['end_faces']))
        for i, face in enumerate(surfaces['end_faces']):
            face_point = face.pointOn[0]
            print("  End face %d: y=%.3f (height=%.3f, -0.9*height=%.3f)" % (i, face_point[1], height, -0.9*height))
            if face_point[1] < -0.9 * height:  # Bottom faces
                all_bottom_faces.append(face)
                print("    -> Added to all_bottom_faces")
            else:
                print("    -> Rejected (y not < -0.9*height)")
    print("Total bottom faces found: %d" % len(all_bottom_faces))
    
    # Separate symmetry faces into thetaStart and thetaEnd groups
    symmetry_faces_thetaStart = []
    symmetry_faces_thetaEnd = []
    
    if 'symmetry_planes' in surfaces and surfaces['symmetry_planes']:
        for sym_face in surfaces['symmetry_planes']:
            face_point = sym_face.pointOn[0]
            x, y, z = face_point
            
            # Determine which symmetry plane this face belongs to
            angle = math.atan2(z, x)
            
            # Check if it's on thetaStart or thetaEnd plane (with tolerance)
            angle_thetaStart = math.pi/2 + thetaStart
            angle_thetaEnd = math.pi/2 + thetaEnd
            
            is_thetaStart_plane = abs(angle - angle_thetaStart) < 1e-1
            is_thetaEnd_plane = abs(angle - angle_thetaEnd) < 1e-1
            
            if is_thetaStart_plane:
                symmetry_faces_thetaStart.append(sym_face)
            elif is_thetaEnd_plane:
                symmetry_faces_thetaEnd.append(sym_face)
    
    # Process ALL bottom faces
    for bottom_face in all_bottom_faces:
        bottom_faces_no_symmetry.append(bottom_face)
        face_point = bottom_face.pointOn[0]
        print("Added bottom face to BC list: x=%.3f, y=%.3f, z=%.3f" % face_point)
    
    print("Total faces for bottom BC: %d" % len(bottom_faces_no_symmetry))
    
    # Use all symmetry faces without exclusion
    symmetry_faces_thetaStart_clean = symmetry_faces_thetaStart
    symmetry_faces_thetaEnd_clean = symmetry_faces_thetaEnd
    
    # Create local cylindrical coordinate system for symmetry boundary conditions
    local_csys_cylindrical = curvedCompositePart.DatumCsysByThreePoints(
        coordSysType=CYLINDRICAL,
        name='LocalCSYS_Cylindrical_Symmetry_BC',
        origin=(0.0, 0.0, 0.0),
        point1=(1.0, 0.0, 0.0),
        point2=(0.0, 0.0, 1.0)
    )
    # STEP 2: Apply SYMMETRY boundary conditions to thetaStart faces
    if symmetry_faces_thetaStart_clean:
        symmetry_region_thetaStart = create_face_region_safely(symmetry_faces_thetaStart_clean, curvedCompositeInstance)
        compositeModel.DisplacementBC(
            name='BC_Symmetry_ThetaStart',
            createStepName='Static_Step',
            region=symmetry_region_thetaStart,
            u2=0.0,   # Circumferential displacement constrained in cylindrical coordinate system
            ur1=0.0,  # Rotation about radial axis
            ur3=0.0,  # Rotation about axial axis
            amplitude=UNSET,
            distributionType=UNIFORM,
            localCsys=curvedCompositeInstance.datums[local_csys_cylindrical.id]
        )
    # STEP 3: Apply SYMMETRY boundary conditions to thetaEnd faces
    if symmetry_faces_thetaEnd_clean:
        symmetry_region_thetaEnd = create_face_region_safely(symmetry_faces_thetaEnd_clean, curvedCompositeInstance)
        compositeModel.DisplacementBC(
            name='BC_Symmetry_ThetaEnd',
            createStepName='Static_Step',
            region=symmetry_region_thetaEnd,
            u2=0.0,   # Circumferential displacement constrained in cylindrical coordinate system
            ur1=0.0,  # Rotation about radial axis
            ur3=0.0,  # Rotation about axial axis
            amplitude=UNSET,
            distributionType=UNIFORM,
            localCsys=curvedCompositeInstance.datums[local_csys_cylindrical.id]
        )
    # STEP 4: Apply BOTTOM SUPPORT to ALL bottom faces
    if bottom_faces_no_symmetry:
        print("Creating bottom face region with %d faces..." % len(bottom_faces_no_symmetry))
        bottom_region = create_face_region_safely(bottom_faces_no_symmetry, curvedCompositeInstance)
        if bottom_region:
            print("Successfully created bottom face region")
        else:
            print("Failed to create bottom face region!")
        compositeModel.DisplacementBC(
            name='BC_Bottom_Support_Faces',
            createStepName='Static_Step',
            region=bottom_region,
            u2=0.0,   # Vertical displacement
            ur1=0.0,  # Rotation about radial axis
            ur3=0.0,  # Rotation about circumferential axis
            amplitude=UNSET,
            distributionType=UNIFORM
        )
        print("Applied bottom face boundary condition")
    return True

def create_edge_region_safely(edge_list, instance):
    """
    Safely create a region from a list of edges, handling Abaqus GeomSequence requirements
    """
    if not edge_list:
        return None
    
    # Method 1: Use edge indices to create a proper sequence
    edge_indices = []
    for edge in edge_list:
        for i, inst_edge in enumerate(instance.edges):
            if inst_edge.index == edge.index:
                edge_indices.append(i)
                break
    
    if edge_indices:
        # Create a sequence of edges from the instance
        selected_edges = []
        for idx in edge_indices:
            selected_edges.append(instance.edges[idx:idx+1])
        
        # Combine all edge sequences
        if len(selected_edges) == 1:
            combined_edges = selected_edges[0]
        else:
            combined_edges = selected_edges[0]
            for edge_seq in selected_edges[1:]:
                combined_edges = combined_edges + edge_seq
        
        return regionToolset.Region(edges=combined_edges)
    
    return None

def create_edge_boundary_conditions(compositeModel, curvedCompositeInstance, curvedCompositePart, 
                                   surfaces, local_csys, rk, height, thetaStart, thetaEnd):
    """
    STEP 2: Apply boundary conditions to shared edges between different face types.
    These are the EDGE symmetry conditions applied after face symmetries.
    1. Edges shared between thetaStart and thetaEnd symmetry planes  
    2. Edges shared between thetaStart symmetry plane and bottom faces
    3. Edges shared between thetaEnd symmetry plane and bottom faces
    """
    if 'Static_Step' not in compositeModel.steps.keys():
        compositeModel.StaticStep(name='Static_Step', previous='Initial', 
                                  description='Static analysis step')
    
    # 1. Edges shared between thetaStart and thetaEnd symmetry planes (at inner/outer radii)
    symmetry_radial_edges = []
    for r,b in zip(rm, bm):  # Inner and outer radii
#        for y in meridional_Offsets[-1]:  # Axial midpoints
# Edge at the intersection of both symmetry planes (should be at x=0 if symmetric about z-axis)
# For quarter model, this would be at the axis of revolution
        theta=0.0
        m=np.pi/2
        x = ct.pol2cart3D_x_intersect(r, b, m, np.pi/2+theta)
        y = ct.pol2cart3D_y_intersect(r, b, m, np.pi/2+theta)
        z = ct.pol2cart3D_z_intersect(r, b, m, np.pi/2+theta)
        test_point = (x, y, z)
        found_edges = curvedCompositePart.edges.findAt((test_point,))
        if found_edges:
            symmetry_radial_edges.extend(found_edges)
    
    # 2. Edges shared between thetaStart symmetry plane and bottom face
    thetaStart_bottom_edges = []
    for r in rm:  # Use radial midpoints
        x = r * math.cos(math.pi / 2 + thetaStart)
        z = r * math.sin(math.pi / 2 + thetaStart)
        y = -height  # Bottom edge
        test_point = (x, y, z)
        found_edges = curvedCompositePart.edges.findAt((test_point,))
        if found_edges:
            thetaStart_bottom_edges.extend(found_edges)
    
    # 3. Edges shared between thetaEnd symmetry plane and bottom face
    thetaEnd_bottom_edges = []
    for r in rm:  # Use radial midpoints
        x = r * math.cos(math.pi / 2 + thetaEnd)
        z = r * math.sin(math.pi / 2 + thetaEnd)
        y = -height  # Bottom edge
        test_point = (x, y, z)
        found_edges = curvedCompositePart.edges.findAt((test_point,))
        if found_edges:
            thetaEnd_bottom_edges.extend(found_edges)
        
    # Apply edge-based boundary conditions
    
    # BC for symmetry-radial shared edges (thetaStart-thetaEnd intersection)
    if symmetry_radial_edges:
        edge_region = create_edge_region_safely(symmetry_radial_edges, curvedCompositeInstance)
        compositeModel.DisplacementBC(
            name='BC_Edge_Symmetry_Radial',
            createStepName='Static_Step',
            region=edge_region,
            u1=0.0,   # X-displacement constrained
            u3=0.0,   # Z-displacement constrained
            ur1=0.0,  # Rotation about X-axis constrained
            ur2=0.0,  # Rotation about Y-axis constrained
            ur3=0.0,  # Rotation about Z-axis constrained
            amplitude=UNSET,
            distributionType=UNIFORM,
        )
        
    # BC for thetaStart-bottom shared edges
    if thetaStart_bottom_edges:
        edge_region = create_edge_region_safely(thetaStart_bottom_edges, curvedCompositeInstance)
        compositeModel.DisplacementBC(
            name='BC_Edge_ThetaStart_Bottom',
            createStepName='Static_Step',
            region=edge_region,
            u2=0.0,   # Y-displacement constrained
            u3=0.0,   # Z-displacement constrained
            ur1=0.0,  # Rotation about X-axis constrained
            ur2=0.0,  # Rotation about Y-axis constrained
            ur3=0.0,  # Rotation about Z-axis constrained
            amplitude=UNSET,
            distributionType=UNIFORM,
            localCsys=curvedCompositeInstance.datums[local_csys.id]
        )
    
    # BC for thetaEnd-bottom shared edges
    if thetaEnd_bottom_edges:
        edge_region = create_edge_region_safely(thetaEnd_bottom_edges, curvedCompositeInstance)
        compositeModel.DisplacementBC(
            name='BC_Edge_ThetaEnd_Bottom',
            createStepName='Static_Step',
            region=edge_region,
            u2=0.0,   # Y-displacement constrained
            u3=0.0,   # Z-displacement constrained
            ur1=0.0,  # Rotation about X-axis constrained
            ur2=0.0,  # Rotation about Y-axis constrained
            ur3=0.0,  # Rotation about Z-axis constrained
            amplitude=UNSET,
            distributionType=UNIFORM,
            localCsys=curvedCompositeInstance.datums[local_csys.id]
        )
    
    return True

def apply_internal_pressure_load(compositeModel, curvedCompositeInstance, curvedCompositePart, surfaces, InnerPressure):
    """
    Apply internal pressure load to the inner surfaces
    Based on the methodology from BC.py file
    """
    if InnerPressure == 0.0:
        print("No internal pressure specified (InnerPressure = 0.0)")
        return True
    
    try:
        if 'Static_Step' not in compositeModel.steps.keys():
            compositeModel.StaticStep(name='Static_Step', previous='Initial',
                                    description='Static analysis step')
        
        print("Applying internal pressure load: %.3f" % InnerPressure)
        
        # Use the simplified main.py approach for inner pressure faces
        # Build coordinate string for single findAt call on instance faces
        curvedCompositeStringFaces = ''
        
        # Use thetaEdges_mid and yOffsets_mid for comprehensive face coverage
        # This matches the main.py approach but adapted for our 3D geometry
        for theta in thetaEdges_mid:  # Use circumferential midpoints
            for y in yOffsets_mid:  # Use axial midpoints  
                x = rk[0] * np.cos(np.pi/2+theta)
                z = rk[0] * np.sin(np.pi/2+theta)
                curvedCompositeStringFaces += str(((x, y, z),),) + ','
        
        # Add shell coordinates using meridional offsets for complete coverage
        for theta in thetaEdges_mid:  # Use circumferential midpoints
            for m in meridional_Offsets_mid:  # Use meridional midpoints
                x = ct.pol2cart3D_x_intersect(rk[0], bk[0], m, np.pi/2+theta)
                y = ct.pol2cart3D_y_intersect(rk[0], bk[0], m, np.pi/2+theta)
                z = ct.pol2cart3D_z_intersect(rk[0], bk[0], m, np.pi/2+theta)
                curvedCompositeStringFaces += str(((x, y, z),),) + ','
        
        print("Built coordinate string for inner pressure faces")
        print("Testing %d flange points (thetaEdges_mid × yOffsets_mid)" % (len(thetaEdges_mid) * len(yOffsets_mid)))
        print("Testing %d shell points (thetaEdges_mid × meridional_Offsets_mid)" % (len(thetaEdges_mid) * len(meridional_Offsets_mid)))
        
        # Execute single findAt call on instance faces (following main.py approach)
        curvedCompositeStringFacesExec = 'curvedCompositeInnerPressureFaces = curvedCompositeInstance.faces.findAt(' + curvedCompositeStringFaces + ')'
        
        # Create local namespace for exec
        local_vars = {'curvedCompositeInstance': curvedCompositeInstance}
        exec(curvedCompositeStringFacesExec, globals(), local_vars)
        curvedCompositeInnerPressureFaces = local_vars['curvedCompositeInnerPressureFaces']
        
        print("Found %d inner pressure faces using main.py method" % len(curvedCompositeInnerPressureFaces))
        
        # Create surface set for GUI verification (using the found instance faces)
        if curvedCompositeInnerPressureFaces:
            try:
                # Create surface set on the part for verification
                surface_set_name = 'Set_Inner_Faces_Found'
                # Convert instance faces to part faces for surface creation
                part_faces_for_surface = []
                for inst_face in curvedCompositeInnerPressureFaces:
                    try:
                        face_point = inst_face.pointOn[0]
                        part_faces = curvedCompositePart.faces.findAt((face_point,))
                        if part_faces:
                            part_faces_for_surface.extend(part_faces)
                    except Exception:
                        continue
                
                # Delete existing surface if it exists
                try:
                    if surface_set_name in curvedCompositePart.surfaces.keys():
                        del curvedCompositePart.surfaces[surface_set_name]
                        print("Deleted existing surface set: %s" % surface_set_name)
                except Exception:
                    pass
                
                if part_faces_for_surface:
                    curvedCompositePart.Surface(name=surface_set_name, side1Faces=tuple(part_faces_for_surface))
                    print("Created surface set '%s' with %d faces for GUI verification" % (surface_set_name, len(part_faces_for_surface)))
                    print("You can now check this surface set in Abaqus GUI to verify the correct faces are selected")
                
            except Exception as e:
                print("Warning: Could not create surface set for GUI verification: %s" % str(e))
        
        # Apply pressure directly to instance faces (main.py approach - no surface creation needed)
        if curvedCompositeInnerPressureFaces:
            try:
                compositeModel.Pressure(
                    name='Load_Inner_Pressure',
                    createStepName='Static_Step',
                    region=regionToolset.Region(side1Faces=curvedCompositeInnerPressureFaces),
                    magnitude=InnerPressure,
                    distributionType=UNIFORM
                )
                print("✓ Successfully applied inner pressure load: %.3f using main.py method" % InnerPressure)
                return True
                
            except Exception as e:
                print("✗ Failed to apply pressure load: %s" % str(e))
                return False
        else:
            print("✗ No inner pressure faces found")
            return False
                
    except Exception as e:
        print("✗ Error in apply_internal_pressure_load: %s" % str(e))
        return False

def setup_complete_boundary_conditions_and_loads(compositeModel, curvedCompositePart, curvedCompositeInstance, 
                                                curvedCompositeCylCoordSys, rk, bk, phiEdges, height, 
                                                thetaStart, thetaEnd, InnerPressure):
    # First identify all boundary surfaces
    surfaces = identify_boundary_surfaces(compositeModel, curvedCompositePart, rk, bk, 
                                        phiEdges, height, thetaStart, thetaEnd, rm)
    
    # STEP 1: Apply face symmetry boundary conditions first
    face_bc_success = create_boundary_conditions(compositeModel, curvedCompositeInstance, 
                                                surfaces, curvedCompositeCylCoordSys)
    
    # STEP 2: Apply edge symmetry boundary conditions second
    # Create local coordinate system for edge boundary conditions
    local_csys_edge = curvedCompositePart.DatumCsysByThreePoints(
        coordSysType=CYLINDRICAL,
        name='LocalCSYS_Edge_BC',
        origin=(0.0, 0.0, 0.0),
        point1=(1.0, 0.0, 0.0),
        point2=(0.0, 0.0, 1.0)
    )
    
    # Create boundary conditions for shared edges
    # edge_bc_success = create_edge_boundary_conditions(compositeModel, curvedCompositeInstance, 
    #                                                 curvedCompositePart, surfaces, local_csys_edge,
    #                                                 rk, height, thetaStart, thetaEnd)
    
    # STEP 3: Apply internal pressure load
    pressure_success = apply_internal_pressure_load(compositeModel, curvedCompositeInstance, 
                                                   curvedCompositePart, surfaces, InnerPressure)
    
    # Combined success of boundary conditions and pressure load
    overall_success = face_bc_success and pressure_success
    
    print("\nBoundary conditions and loads summary:")
    print("  Face BC success: %s" % face_bc_success)
    print("  Pressure load success: %s" % pressure_success)
    print("  Overall success: %s" % overall_success)
    
    return overall_success

# Apply boundary conditions
bc_success = setup_complete_boundary_conditions_and_loads(
    compositeModel=compositeModel,
    curvedCompositePart=curvedCompositePart,
    curvedCompositeInstance=curvedCompositeInstance,
    curvedCompositeCylCoordSys=curvedCompositeCylCoordSys,
    rk=rk, bk=bk,
    phiEdges=phiEdges,
    height=height,
    thetaStart=thetaStart,
    thetaEnd=thetaEnd,
    InnerPressure=InnerPressure,
)























def setsCylindricalSpheroidAnalysis(compositeAssembly, curvedCompositeInstance, rk, bk, 
                                   height, thetaEnd, N, iInterfaceEval):
    """
    Create evaluation sets for cylindrical spheroid geometry
    Adapts setsCurvedCompositeCirc from 2D_L_CFRP_FEM_circular_path.py for 3D geometry
    Creates FEInner sets for layer analysis and FEInterface sets for interface analysis
    """
    try:
        print("\nCreating evaluation sets for post-processing...")
        
        # Import coordinate transformation functions (if available)
        import coordinateTransformation_ellipse as ct
        print("  ✓ Successfully imported coordinateTransformation_ellipse")        
        # Define evaluation parameters  # Mid-angle for evaluation
        
        print("  Debug: angleEval = %s (type: %s)" % (angleEval, type(angleEval)))
        print("  Debug: rm = %s (length: %d)" % (rm, len(rm)))
        print("  Debug: yOffsets = %s (length: %d)" % (yOffsets, len(yOffsets)))
        
        # Verify required variables exist
        try:
            print("  Debug: Checking required variables...")
            print("    compositeAssembly exists: %s" % (compositeAssembly is not None))
            print("    curvedCompositeInstance exists: %s" % (curvedCompositeInstance is not None))
            print("    np module exists: %s" % ('np' in globals()))
            
            # Test iteration over angleEval to verify it works
            print("  Debug: Testing angleEval iteration...")
            test_count = 0
            for theta in [angleEval]:
                test_count += 1
                print("    Test iteration %d: theta = %s" % (test_count, theta))
            print("  Debug: angleEval iteration test completed with %d iterations" % test_count)
            
        except Exception as e:
            print("  Debug: Error checking variables: %s" % str(e))
            return False
        
        # STEP 1: Create FEInner sets for each layer (for layer-wise analysis)
        print("  Creating FEInner sets for layer analysis...")
        
        for ii, r in enumerate(rm):         
            print("    Debug: Processing layer %d with r = %s" % (ii + 1, r))
            # Build coordinate string using same logic as meshing
            curvedCompositeStringEdges = ''          
            # Use fixed angleEval and all yOffsets for edge coordinates
            print("    Debug: About to iterate over [angleEval] = %s" % [angleEval])
            for theta in [angleEval]:  # Use the specific evaluation angle (make it iterable)
                print("      Debug: Processing theta = %s" % theta)
                try:
                    x = ct.cyl2cart_x(r, theta)
                    y = yOffsets[-1]  # Use the specific layer's yOffset
                    z = ct.cyl2cart_z(r, theta)
                    print("        Debug: Calculated x=%s, z=%s" % (x, z))
                    curvedCompositeStringEdges += str(((x, y, z),),) + ','
                except Exception as e:
                    print("        Debug: Error in coordinate calculation: %s" % str(e))
                    raise e
            
            print("    Debug: Finished coordinate generation, string length: %d" % len(curvedCompositeStringEdges))
            print("    Debug: First 100 chars of coordinate string: %s..." % curvedCompositeStringEdges[:100])
            
            try:
                # Execute dynamic edge finding (same as meshing logic)
                curvedCompositeStringEdgesExec = ('compositeSetEvalEdge = ' +
                                                'curvedCompositeInstance.edges.findAt(' + 
                                                curvedCompositeStringEdges + ')')
                
                print("    Debug: About to execute edge finding...")
                # Create local namespace for exec
                local_vars = {'curvedCompositeInstance': curvedCompositeInstance}
                exec(curvedCompositeStringEdgesExec, globals(), local_vars)
                compositeSetEvalEdge = local_vars['compositeSetEvalEdge']
                
                print("    Debug: Found %d edges" % len(compositeSetEvalEdge))
                
                # Create set for this layer
                compositeAssembly.Set(edges=compositeSetEvalEdge, name='FEInner' + str(ii + 1))
                print("    ✓ Created FEInner%d set with %d edges" % (ii + 1, len(compositeSetEvalEdge)))
                
            except Exception as e:
                print("    ⚠ Could not create FEInner%d set: %s" % (ii + 1, str(e)))
        
        # STEP 2: Create FEInterface sets for interface analysis (between layers)
        print("  Creating FEInterface sets for interface analysis...")
        
        for kk in range(len(iInterfaceEval)):
            interface_idx = iInterfaceEval[kk]
            r_interface = rk[interface_idx]
            b_interface = bk[interface_idx]
            
            # Build coordinate string using same logic as meshing (but for interface radius)
            curvedCompositeStringEdges = ''
            
            # Use thetaEdges_mid and yOffsets_mid for flange coordinates
            for theta in [angleEval]:  # Use circumferential midpoints
                for y in yOffsets_mid:  # Use axial midpoints  
                    x = r_interface * np.cos(theta)
                    z = r_interface * np.sin(theta)
                    curvedCompositeStringEdges += str(((x, y, z),),) + ','
            
            # Use thetaEdges_mid and meridional_Offsets_mid for shell coordinates
            for theta in [angleEval]:  # Use circumferential midpoints
                for m in meridional_Offsets_mid:  # Use meridional midpoints
                    x = ct.pol2cart3D_x_intersect(r_interface, b_interface, m, theta)
                    y = ct.pol2cart3D_y_intersect(r_interface, b_interface, m, theta)
                    z = ct.pol2cart3D_z_intersect(r_interface, b_interface, m, theta)
                    curvedCompositeStringEdges += str(((x, y, z),),) + ','
            
            try:
                # Execute dynamic edge finding (same as meshing logic)
                curvedCompositeStringEdgesExec = ('compositeSetEdgeInterface = ' +
                                                'curvedCompositeInstance.edges.findAt(' + 
                                                curvedCompositeStringEdges + ')')
                
                # Create local namespace for exec
                local_vars = {'curvedCompositeInstance': curvedCompositeInstance}
                exec(curvedCompositeStringEdgesExec, globals(), local_vars)
                compositeSetEdgeInterface = local_vars['compositeSetEdgeInterface']
                
                # Create interface set
                compositeAssembly.Set(edges=compositeSetEdgeInterface, 
                                    name='FEInterface' + str(interface_idx))
                print("    ✓ Created FEInterface%d set with %d edges" % (interface_idx, len(compositeSetEdgeInterface)))
                
            except Exception as e:
                print("    ⚠ Could not create FEInterface%d set: %s" % (interface_idx, str(e)))
        print("  ✓ Set creation completed")
        return True
        
    except Exception as e:
        print("  ✗ Error creating evaluation sets: %s" % str(e))
        return False

sets_success = setsCylindricalSpheroidAnalysis(
    compositeAssembly=compositeAssembly,
    curvedCompositeInstance=curvedCompositeInstance,
    rk=rk,
    bk=bk,
    height=height,
    thetaEnd=thetaEnd,
    N=N,
    iInterfaceEval=iInterfaceEval,
)

print("Evaluation sets creation: %s" % sets_success)

def create_and_submit_job(compositeModel, modelName, analysis_newPath):
    """
    Create and submit Abaqus job with the same parameters as main.py
    """
    try:
        # Job creation parameters matching main.py
        jobName = modelName + '_Job'
        
        print("Creating job: %s" % jobName)
        
        # Create job with same parameters as main.py
        mdb.Job(name=jobName, 
                model=modelName,
                description='Run FE-analysis', 
                parallelizationMethodExplicit=DOMAIN, 
                numDomains=12, 
                numCpus=6, 
                memory=85,
                echoPrint=OFF, 
                modelPrint=OFF, 
                contactPrint=OFF, 
                historyPrint=OFF)
        
        print("Job created successfully")
        
        # Save the CAE file before submitting job
        print("Saving CAE file...")
        mdb.saveAs(pathName=analysis_newPath + '\\' + modelName + '_CAE_Model')
        print("CAE file saved: %s" % (analysis_newPath + '\\' + modelName + '_CAE_Model.cae'))
        
        # Submit the job
        print("Submitting job: %s" % jobName)
        mdb.jobs[jobName].submit(consistencyChecking=OFF)
        
        print("Job submitted successfully")
        print("Waiting for job completion...")
        
        # Wait for job completion
        mdb.jobs[jobName].waitForCompletion()
        
        print("Job completed successfully: %s" % jobName)
        return True
        
    except Exception as e:
        print("Error in create_and_submit_job: %s" % str(e))
        return False

# Create and submit the job
job_success = create_and_submit_job(compositeModel, modelName, analysis_newPath)
print("Job submission success: %s" % job_success)


def cylindrical_spherical_elemDefLocCenGrav_TS_TG_ZS_ZG(thetaEval, zEval, stressLastFrameElementNodal, compositeSetElements):
    """
    Assign elements to TS_ZS, TS_ZG, TG_ZS, TG_ZG by centroid theta and z coordinates.
    CYLINDRICAL PRESSURE VESSEL METHODOLOGY - designed for pressure vessel FEInner sets
    TS/TG: Theta Smaller/Greater relative to thetaEval
    ZS/ZG: Z Smaller/Greater relative to zEval
    Returns: elements (dict), elementLabels (dict)
    """
    elements = {'TS_ZS': [], 'TS_ZG': [], 'TG_ZS': [], 'TG_ZG': []}
    elementLabels = {'TS_ZS': [], 'TS_ZG': [], 'TG_ZS': [], 'TG_ZG': []}
    elementsExist = []
    
    print("    Classifying elements by TS/TG and ZS/ZG (cylindrical pressure vessel methodology)...")
    
    # OPTIMIZATION: Convert elements to hash map for O(1) lookup instead of O(n) iteration
    try:
        elementsList = list(compositeSetElements.elements[0])  # Access the first (and only) sequence
        elementsMap = {elem.label: elem for elem in elementsList}  # Create hash map
        print("    DEBUG: Successfully converted %d elements to hash map" % len(elementsMap))
    except Exception as e:
        print("    DEBUG: Error converting elements to hash map: %s" % str(e))
        return elements, elementLabels
    
    # OPTIMIZATION: Process stress values only once, use hash lookup for elements
    processed_element_labels = set()
    for stress_value in stressLastFrameElementNodal:
        elem_label = stress_value.elementLabel
        if elem_label in processed_element_labels:
            continue  # Skip already processed elements
        processed_element_labels.add(elem_label)
        
        # O(1) hash lookup instead of O(n) iteration
        if elem_label in elementsMap:
            element = elementsMap[elem_label]
            try:
                # Ensure connectivity is properly converted to list
                connectivity_raw = element.connectivity
                nodeLabelConnect = np.array(list(connectivity_raw))
            except Exception as e:
                print("    DEBUG: Error in connectivity conversion: %s" % str(e))
                continue
            
            try:
                nodeLabelConnectInst = stress_value.instance
            except Exception as e:
                print("    DEBUG: Error getting instance: %s" % str(e))
                continue
            
            try:
                # Get node coordinates
                xCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[0] for kk in nodeLabelConnect])
                yCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[1] for kk in nodeLabelConnect])
                zCoord = np.array([nodeLabelConnectInst.getNodeFromLabel(kk).coordinates[2] for kk in nodeLabelConnect])
            except Exception as e:
                print("    DEBUG: Error getting coordinates: %s" % str(e))
                continue
            
            try:
                # Calculate cylindrical coordinates (theta from X-Z plane, Z is Y-coordinate)
                thetaCoord = np.arctan2(zCoord, xCoord)  # Theta angle from X-axis
                zCylCoord = yCoord  # Z coordinate is Y in global system
                
                thetaCentroid = np.mean(thetaCoord)
                zCentroid = np.mean(zCylCoord)
            except Exception as e:
                print("    DEBUG: Error calculating centroids: %s" % str(e))
                continue
            
            try:
                # TS/TG split by theta angle, ZS/ZG split by z coordinate
                if thetaCentroid <= thetaEval:
                    if zCentroid <= zEval:
                        elements['TS_ZS'].append(stress_value)
                        elementLabels['TS_ZS'].append(stress_value.elementLabel)
                    else:
                        elements['TS_ZG'].append(stress_value)
                        elementLabels['TS_ZG'].append(stress_value.elementLabel)
                else:
                    if zCentroid <= zEval:
                        elements['TG_ZS'].append(stress_value)
                        elementLabels['TG_ZS'].append(stress_value.elementLabel)
                    else:
                        elements['TG_ZG'].append(stress_value)
                        elementLabels['TG_ZG'].append(stress_value.elementLabel)
            except Exception as e:
                print("    DEBUG: Error classifying element: %s" % str(e))
                continue
    
    try:
        print("    DEBUG: About to print classification results")
        # Print classification results with explicit string conversion to avoid generator issues
        ts_zs_count = len(elements['TS_ZS'])
        ts_zg_count = len(elements['TS_ZG'])
        tg_zs_count = len(elements['TG_ZS'])
        tg_zg_count = len(elements['TG_ZG'])
        total_elements = ts_zs_count + ts_zg_count + tg_zs_count + tg_zg_count
        
        print("    ✓ Classified " + str(total_elements) + " elements:")
        print("      TS_ZS: " + str(ts_zs_count) + " elements")
        print("      TS_ZG: " + str(ts_zg_count) + " elements")
        print("      TG_ZS: " + str(tg_zs_count) + " elements")
        print("      TG_ZG: " + str(tg_zg_count) + " elements")
        print("    DEBUG: Successfully printed classification results")

        # DETAILED ELEMENT AND NODE ANALYSIS: Find nodes with (x≈0, y≈0) in each group
        print("      ===== DETAILED ELEMENT AND NODE ANALYSIS (OPTIMIZED) =====")
        axis_tol = 1e-6
        for group_name in ['TS_ZS', 'TS_ZG', 'TG_ZS', 'TG_ZG']:
            group_elems = elements[group_name]
            if not group_elems:
                continue
            print("      DEBUG: Group {} has {} elements".format(group_name, len(group_elems)))
            # OPTIMIZATION: Only analyze first element per group to avoid excessive output
            for elem in group_elems[:1]:  # Process only first element
                if not hasattr(elem, 'elementLabel'):
                    continue
                elemLabel = elem.elementLabel
                print("      DEBUG:   Sample Element {} from group {}".format(elemLabel, group_name))
        print("      ===== END DETAILED ANALYSIS =====")
    except Exception as e:
        print("    DEBUG: Error in printing classification results: " + str(e))
    
    return elements, elementLabels

def postProcessFEInner_TS_TG_ZS_ZG(odbPath, stepName, thetaEval, zEval, N):
    """
    Post-process FEInner sets with TS/TG and ZS/ZG classification
    CYLINDRICAL PRESSURE VESSEL METHODOLOGY - designed for cylindrical pressure vessel FEInner sets
    """
    try:
        print("\nStarting TS/TG and ZS/ZG post-processing (cylindrical pressure vessel methodology)...")
        
        # Open ODB
        odbObject = session.openOdb(name=odbPath)
        
        # Stress components for output
        stress_components = ['S11', 'S22', 'S33', 'S12', 'S13', 'S23']
        
        # Process each layer
        for layerIndex in range(N):
            layer_name = "FEInner{}".format(layerIndex + 1)
            print("  Processing layer %d (%s)..." % (layerIndex + 1, layer_name))
            
            try:
                # Get element set
                elementSet = odbObject.rootAssembly.elementSets[layer_name.upper()]
                lastFrame = odbObject.steps[stepName].frames[-1]
                
                # ELEMENT_NODAL stress (list of FieldValue objects)
                stressField = lastFrame.fieldOutputs['S'].getSubset(region=elementSet, position=ELEMENT_NODAL)
                
                # Convert generator to list to avoid "expecting a recognized type" error
                stressValues = list(stressField.values)
                
                print("    DEBUG: About to call cylindrical_spherical_elemDefLocCenGrav_TS_TG_ZS_ZG...")
                
                # Classify elements by geometry regions (cylindrical pressure vessel methodology)
                elements_dict, elementLabels_dict = cylindrical_spherical_elemDefLocCenGrav_TS_TG_ZS_ZG(
                    thetaEval, zEval, stressValues, elementSet
                )
                
                print("    DEBUG: Classification completed successfully!")
                print("    DEBUG: elements_dict type:", type(elements_dict))
                print("    DEBUG: elements_dict keys:", elements_dict.keys())
                for key in elements_dict.keys():
                    print("    DEBUG: elements_dict[%s] type:" % key, type(elements_dict[key]))
                    print("    DEBUG: elements_dict[%s] length:" % key, len(elements_dict[key]))
                    if len(elements_dict[key]) > 0:
                        print("    DEBUG: elements_dict[%s][0] type:" % key, type(elements_dict[key][0]))
                
                # Process each region (TS/TG + ZS/ZG combinations)
                for group in ['TS_ZS', 'TS_ZG', 'TG_ZS', 'TG_ZG']:
                    if elements_dict[group]:
                        print("    Processing group %s (%d elements)..." % (group, len(elements_dict[group])))
                        
                        # Collect stress data for this group
                        group_results = []
                        
                        print("    DEBUG: About to iterate through elements_dict[%s]..." % group)
                        for stress_value in elements_dict[group]:
                            try:
                                node_label = stress_value.nodeLabel
                                element_label = stress_value.elementLabel
                                stress_data = stress_value.data
                                
                                # Store stress data
                                result_entry = {
                                    'nodeLabel': node_label,
                                    'elementLabel': element_label,
                                    'group': group,
                                    'layer': layerIndex + 1,
                                    'stresses': {
                                        'S11': stress_data[0],
                                        'S22': stress_data[1], 
                                        'S33': stress_data[2],
                                        'S12': stress_data[3],
                                        'S13': stress_data[4],
                                        'S23': stress_data[5]
                                    }
                                }
                                group_results.append(result_entry)
                                
                            except Exception as e:
                                print("      Warning: Error processing stress data: %s" % str(e))
                                continue
                        
                        # Write results to file
                        if group_results:
                            filename = layer_name + '_' + group + '.txt'
                            
                            with open(filename, 'w') as f:
                                f.write("# Layer %d Group %s Stress Results (3D_plate_v1 methodology)\n" % (layerIndex + 1, group))
                                f.write("# NodeLabel\tElementLabel\tS11\tS22\tS33\tS12\tS13\tS23\n")
                                
                                for result in group_results:
                                    f.write("%d\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % (
                                        result['nodeLabel'],
                                        result['elementLabel'],
                                        result['stresses']['S11'],
                                        result['stresses']['S22'],
                                        result['stresses']['S33'],
                                        result['stresses']['S12'],
                                        result['stresses']['S13'],
                                        result['stresses']['S23']
                                    ))
                        
                            print("      ✓ Results written to: %s" % filename)
                        else:
                            print("      ⚠ No results found for group %s" % group)
                    else:
                        print("    ⚠ No elements found in group %s" % group)
                
            except Exception as e:
                print("    ✗ Error processing layer %d: %s" % (layerIndex + 1, str(e)))
                continue
        
        # Close ODB
        odbObject.close()
        
        print("  ✓ TS/TG and ZS/ZG post-processing completed (3D_plate_v1 methodology)")
        return True
        
    except Exception as e:
        print("  ✗ Error in TS/TG post-processing: %s" % str(e))
        return False


def setup_dual_postprocessing_coordinate_systems(odbObject, rk, height):
    """
    Create both cylindrical and spherical coordinate systems for post-processing
    Based on 2D_L_CFRP_FEM_circular_path.py dual coordinate approach
    """
    try:
        # Cylindrical coordinate system for cylindrical shell region
        postProc_Cyl_CS = odbObject.rootAssembly.DatumCsysByThreePoints(
            coordSysType=CYLINDRICAL, 
            name='PostProcessing_Cylindrical_CS',
            point1=(rk[-1], 0, 0),
            point2=(0, 0, rk[-1]),
            origin=(0, 0, 0)
        )
        print("✓ Created cylindrical coordinate system for cylindrical shell region")
        
        # Spherical coordinate system for hemispherical dome region  
        postProc_Sph_CS = odbObject.rootAssembly.DatumCsysByThreePoints(
            coordSysType=SPHERICAL,
            name='PostProcessing_Spherical_CS', 
            point1=(rk[-1], 0, 0),
            point2=(0, 0,rk[-1]),
            origin=(0, 0, 0)
        )
        print("✓ Created spherical coordinate system for hemispherical dome region")
        
        return postProc_Cyl_CS, postProc_Sph_CS
        
    except Exception as e:
        print("✗ Error creating post-processing coordinate systems: %s" % str(e))
        return None, None

def classify_elements_within_set(elementSet, odbObject, stepName):
    """
    Classify elements within a specific element set into cylindrical vs spherical regions
    Based on element centroid Y-coordinate: Y<0 = cylindrical, Y>=0 = hemispherical
    """
    try:
        cylindrical_elements = []
        spherical_elements = []
        
        # Get stress field to access instances (same way as 3D_plate_v1)
        lastFrame = odbObject.steps[stepName].frames[-1]
        stressField = lastFrame.fieldOutputs['S'].getSubset(region=elementSet, position=ELEMENT_NODAL)
        
        print("    Processing %d elements for classification..." % len(elementSet.elements[0]))
        
        for element in elementSet.elements[0]:
            # Get element node coordinates using instance method (same as 3D_plate_v1)
            nodeLabelConnect = list(element.connectivity)
            
            # Find the instance for this element by checking stress field values
            instance = None
            for stress_value in stressField.values:
                if stress_value.elementLabel == element.label:
                    instance = stress_value.instance
                    break
            
            if instance is not None:
                # Calculate element centroid from node coordinates (same method as 3D_plate_v1)
                node_coords = []
                for nodeLabel in nodeLabelConnect:
                    coords = instance.getNodeFromLabel(nodeLabel).coordinates
                    node_coords.append(coords)
                
                if node_coords and len(node_coords) > 0:
                    # Calculate centroid Y-coordinate
                    centroid_y = sum([coord[1] for coord in node_coords]) / len(node_coords)
                    
                    # Classification based on y-coordinate
                    # Elements in cylindrical shell: y <= 0.0 (straight walls + boundary)  
                    # Elements in hemispherical dome: y > 0.0 (curved dome only)
                    # Note: Y=0 boundary nodes should be processed as cylindrical
                    if centroid_y <= 0.0:
                        cylindrical_elements.append(element.label)
                    else:
                        spherical_elements.append(element.label)
        
        return cylindrical_elements, spherical_elements
        
    except Exception as e:
        print("    ✗ Error classifying elements in set: %s" % str(e))
        return [], []

def classify_elements_by_geometry_region(compositeSetElements, odbObject, stepName):
    """
    Classify elements into cylindrical shell vs hemispherical dome regions
    Based on element centroid Y-coordinate: Y<0 = cylindrical, Y>=0 = hemispherical
    """
    try:
        cylindrical_elements = []
        spherical_elements = []
        
        print("  Classifying elements by geometry region...")
        
        # Get stress field to access instances (same approach as 3D_plate_v1)
        lastFrame = odbObject.steps[stepName].frames[-1]
        stressField = lastFrame.fieldOutputs['S'].getSubset(region=compositeSetElements, position=ELEMENT_NODAL)
        
        for element in compositeSetElements.elements[0]:
            # Get element node coordinates using instance method (same as 3D_plate_v1)
            nodeLabelConnect = list(element.connectivity)
            
            # Find the instance for this element by checking stress field values
            instance = None
            for stress_value in stressField.values:
                if stress_value.elementLabel == element.label:
                    instance = stress_value.instance
                    break
            
            if instance is not None:
                # Calculate element centroid from node coordinates (same method as 3D_plate_v1)
                node_coords = []
                for nodeLabel in nodeLabelConnect:
                    coords = instance.getNodeFromLabel(nodeLabel).coordinates
                    node_coords.append(coords)
                
                if node_coords and len(node_coords) > 0:
                    # Calculate centroid
                    centroid_x = sum([coord[0] for coord in node_coords]) / len(node_coords)
                    centroid_y = sum([coord[1] for coord in node_coords]) / len(node_coords)  
                    centroid_z = sum([coord[2] for coord in node_coords]) / len(node_coords)
                
                # Classification based on y-coordinate
                # Elements in cylindrical shell: y <= 0 (straight walls + boundary)
                # Elements in hemispherical dome: y > 0 (curved dome only)
                # Note: Y=0 boundary nodes should be processed as cylindrical
                if centroid_y <= 0.0:
                    cylindrical_elements.append(element.label)
                else:
                    spherical_elements.append(element.label)
        
        print("  ✓ Classified %d cylindrical shell elements" % len(cylindrical_elements))
        print("  ✓ Classified %d hemispherical dome elements" % len(spherical_elements))
        
        return cylindrical_elements, spherical_elements
        
    except Exception as e:
        print("  ✗ Error classifying elements: %s" % str(e))
        return [], []

def process_group_cartesian_zline(odbObject, stepName, elementsInGroup, setName, group, stressField, transform_name, postProc_Cyl_CS=None, postProc_Sph_CS=None):
    """Process group elements and average stresses from both element groups at shared nodes"""
    try:
        print("      DEBUG: Starting processing for group %s" % group) 
        tol_cart = 1e-3
        node_data = {}
        
        # Get transformed stress field based on region 
        if transform_name == 'Cylindrical':
            if postProc_Cyl_CS is None:
                raise ValueError("Cylindrical coordinate system not provided")
            stressField = stressField.getTransformedField(datumCsys=postProc_Cyl_CS)
        elif transform_name == 'Spherical':
            if postProc_Sph_CS is None:
                raise ValueError("Spherical coordinate system not provided")  
            stressField = stressField.getTransformedField(datumCsys=postProc_Sph_CS)

        # Get complementary group name for averaging
        complementary_group = {
            'TS_ZS': 'TG_ZS',
            'TG_ZS': 'TS_ZS', 
            'TS_ZG': 'TG_ZG',
            'TG_ZG': 'TS_ZG'
        }[group]

        # OPTIMIZATION 1: Pre-filter by building node coordinate cache for axis nodes
        print("      Building axis node coordinate cache...")
        axis_node_coords = {}  # nodeLabel -> (x, y, z)
        axis_nodes_in_elements = set()  # Only nodes that belong to our elements
        
        # First pass: Build coordinate cache for nodes in our elements only
        for s in stressField.values:
            if hasattr(s, 'nodeLabel') and s.elementLabel in elementsInGroup:
                if s.nodeLabel not in axis_node_coords:
                    coords = s.instance.getNodeFromLabel(s.nodeLabel).coordinates
                    x, y, z = coords[0], coords[1], coords[2]
                    axis_node_coords[s.nodeLabel] = (x, y, z)
                    # Pre-filter: only cache nodes on axis
                    if abs(x) < tol_cart and abs(y) < tol_cart:
                        axis_nodes_in_elements.add(s.nodeLabel)
        
        print("      Found %d axis nodes in %d total element nodes" % (len(axis_nodes_in_elements), len(axis_node_coords)))
        
        # OPTIMIZATION 2: Process only pre-filtered axis nodes
        for s in stressField.values:
            if hasattr(s, 'nodeLabel') and s.nodeLabel in axis_nodes_in_elements:
                x, y, z = axis_node_coords[s.nodeLabel]  # Use cached coordinates
                pos_key = (round(x, 6), round(y, 6), round(z, 6))
                if pos_key not in node_data:
                    node_data[pos_key] = {'TS': [], 'TG': []}
                
                # Store stress based on element group (elements already filtered)
                if 'TS' in group:
                    node_data[pos_key]['TS'].append((x, y, z, s.data, s.elementLabel))
                else:
                    node_data[pos_key]['TG'].append((x, y, z, s.data, s.elementLabel))

        # Now process data, averaging between TS/TG groups at same position
        nodeData_sorted = []
        for pos_key, group_data in node_data.items():
            if group_data['TS'] or group_data['TG']:  # If we have data from either group
                mesh_x = group_data['TS'][0][0] if group_data['TS'] else group_data['TG'][0][0]
                mesh_y = group_data['TS'][0][1] if group_data['TS'] else group_data['TG'][0][1] 
                mesh_z = group_data['TS'][0][2] if group_data['TS'] else group_data['TG'][0][2]

                # Average stresses from both groups
                all_stresses = []
                if group_data['TS']:
                    all_stresses.extend([v[3] for v in group_data['TS']])
                if group_data['TG']:
                    all_stresses.extend([v[3] for v in group_data['TG']])
                
                if all_stresses:
                    avgStress = np.mean(all_stresses, axis=0)
                    # Use first element label as reference
                    elemLabel = group_data['TS'][0][4] if group_data['TS'] else group_data['TG'][0][4]
                    nodeData_sorted.append((elemLabel, pos_key, mesh_x, mesh_y, mesh_z, avgStress))

        # Sort by z coordinate
        nodeData_sorted.sort(key=lambda row: row[4])

        if nodeData_sorted:
            filename = setName + '_' + group + '.txt'
            with open(filename, 'w') as file:
                # Update header to show correct stress components                
                for row in nodeData_sorted:
                    elemLabel, pos_key, mesh_x, mesh_y, mesh_z, stress = row
                    file.write("{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{:.14e}\t{}\n".format(
                        stress[0], stress[1], stress[2],
                        stress[3], stress[4], stress[5],
                        mesh_x, mesh_y, mesh_z, elemLabel
                    ))
            print("      ✓ Results written: %s (%d positions)" % (filename, len(nodeData_sorted)))
        else:
            print("      DEBUG: No data to write for group %s" % group)
    except Exception as e:
        print("      ✗ Error in processing: %s" % str(e))

def process_feinner_group_with_dual_coordinates(odbObject, stepName, postProc_Cyl_CS, postProc_Sph_CS,
                                                cylindrical_elements, spherical_elements, setName, group,
                                                elements_dict, stressField, thetaEval, zEval):
    """
    Process a specific TS/TG ZS/ZG group with node-based averaging and dual coordinate systems
    Uses cylindrical pressure vessel methodology designed for FEInner sets
    """
    try:
        # Get elements in this group - OPTIMIZATION: Convert to set for O(1) lookup
        elemLabelsInGroup = set([fv.elementLabel for fv in elements_dict[group]])
        if not elemLabelsInGroup:
            print("    No elements found in group %s" % group)
            return True
        
        print("    Processing group %s with %d elements..." % (group, len(elemLabelsInGroup)))

        lastFrame = odbObject.steps[stepName].frames[-1]
        # Always use the original ELEMENT_NODAL field for node selection and output
        stressField_nodal = lastFrame.fieldOutputs['S'].getSubset(position=ELEMENT_NODAL)

        cyl_elements_in_group = [elem for elem in elemLabelsInGroup if elem in cylindrical_elements]
        sph_elements_in_group = [elem for elem in elemLabelsInGroup if elem in spherical_elements]

        # Cylindrical
        if cyl_elements_in_group:
            process_group_cartesian_zline(
                odbObject, stepName, cyl_elements_in_group, setName, group, stressField_nodal, 
                'Cylindrical', postProc_Cyl_CS, postProc_Sph_CS
            )
        # Spherical
        if sph_elements_in_group:
            process_group_cartesian_zline(
                odbObject, stepName, sph_elements_in_group, setName, group, stressField_nodal, 
                'Spherical', postProc_Cyl_CS, postProc_Sph_CS
            )
        return True
    except Exception as e:
        print("    ✗ Error processing group %s: %s" % (group, str(e)))
        return False

def extract_stress_spherical_region(odbObject, stepName, postProc_Sph_CS,
                                   spherical_elements, layer_index, extraction_points, output_filename):
    """
    Extract stress results in spherical coordinates for dome region  
    Based on coordinate transformation methodology for spherical stress extraction
    """
    try:
        print("  Extracting hemispherical dome stresses...")
        
        lastFrame = odbObject.steps[stepName].frames[-1]
        
        # Get stress field transformed to spherical coordinates
        stressField_Sph = lastFrame.fieldOutputs['S'].getTransformedField(datumCsys=postProc_Sph_CS)
        coordField = lastFrame.fieldOutputs['COORD']
        
        results = []
        stress_components_sph = ['S_RR', 'S_TT', 'S_PP', 'S_RT', 'S_RP', 'S_TP']  # Spherical stress components
        
        for point_name, (x_target, y_target, z_target, tolerance) in extraction_points.items():
            print("    Extracting spherical stresses at %s: (%.3f, %.3f, %.3f)" % (point_name, x_target, y_target, z_target))
            
            point_results = []
            
            # Find stress values near the target point in spherical elements
            for stress_value in stressField_Sph.values:
                if stress_value.elementLabel in spherical_elements:
                    # Get node coordinates
                    node_coord = None
                    for coord_value in coordField.values:
                        if coord_value.nodeLabel == stress_value.nodeLabel:
                            node_coord = coord_value.data
                            break
                    
                    if node_coord:
                        x, y, z = node_coord
                        distance = ((x - x_target)**2 + (y - y_target)**2 + (z - z_target)**2)**0.5
                        
                        if distance <= tolerance:
                            # Convert spherical stress tensor to named components
                            stress_dict = {
                                'S_RR': stress_value.data[0],    # Radial stress
                                'S_TT': stress_value.data[1],    # Meridional stress (theta direction)
                                'S_PP': stress_value.data[2],    # Circumferential stress (phi direction)  
                                'S_RT': stress_value.data[3],    # Radial-meridional shear
                                'S_RP': stress_value.data[4],    # Radial-circumferential shear
                                'S_TP': stress_value.data[5]     # Meridional-circumferential shear
                            }
                            point_results.append((x, y, z, stress_dict))
            
            if point_results:
                # Average results if multiple points found
                avg_stresses = {comp: sum([res[3][comp] for res in point_results]) / len(point_results) 
                              for comp in stress_components_sph}
                results.append((point_name, 'SPHERICAL', avg_stresses))
                print("    ✓ Found %d stress points for %s (spherical)" % (len(point_results), point_name))
            else:
                print("    ⚠ No stress data found near %s (spherical)" % point_name)
        
        # Write spherical results to file
        if results:
            with open(output_filename.replace('.txt', '_spherical.txt'), 'w') as f:
                f.write("# Hemispherical Dome Region Stress Results\n") 
                f.write("# Point_Name\tRegion_Type\tS_RR\tS_TT\tS_PP\tS_RT\tS_RP\tS_TP\n")
                for point_name, region_type, stresses in results:
                    f.write("%s\t%s\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % 
                           (point_name, region_type, stresses['S_RR'], stresses['S_TT'], stresses['S_PP'],
                            stresses['S_RT'], stresses['S_RP'], stresses['S_TP']))
            print("  ✓ Spherical stress results written to: %s" % output_filename.replace('.txt', '_spherical.txt'))
        
        return True
        
    except Exception as e:
        print("  ✗ Error extracting spherical stresses: %s" % str(e))
        return False

def extract_stress_results_by_layer_and_location(odbObject, stepName, layer_index, extraction_points, output_filename):
    """
    Extract stress results at specific locations for a given layer
    Combines approaches from main.py (layer-based) and 3D_plate_v1.py (location-based)
    Updated for dual coordinate system approach
    """
    try:
        # Get the last frame (final results)
        lastFrame = odbObject.steps[stepName].frames[-1]
        
        # Get stress field for all elements
        stressField = lastFrame.fieldOutputs['S']
        coordField = lastFrame.fieldOutputs['COORD']
        
        # For compatibility - simplified extraction without coordinate transformation
        results = []
        stress_components = ['S11', 'S22', 'S33', 'S12', 'S13', 'S23']
        
        for point_name, (x_target, y_target, z_target, tolerance) in extraction_points.items():
            print("  Extracting results near point %s: (%.3f, %.3f, %.3f)" % (point_name, x_target, y_target, z_target))
            
            point_results = []
            
            # Find nodes/elements near the target point
            for i, stress_value in enumerate(stressField.values):
                try:
                    # Get coordinates for this stress point
                    coord_value = coordField.values[i]
                    x, y, z = coord_value.data
                    
                    # Check if point is within tolerance
                    if (abs(x - x_target) < tolerance and 
                        abs(y - y_target) < tolerance and 
                        abs(z - z_target) < tolerance):
                        
                        # Extract stress components
                        stress_data = stress_value.data
                        point_results.append({
                            'nodeLabel': stress_value.nodeLabel,
                            'elementLabel': stress_value.elementLabel,
                            'coordinates': (x, y, z),
                            'stress': stress_data,
                            'location': point_name
                        })
                        
                except Exception:
                    continue
        
        # Write results to file
        if results:
            with open(output_filename, 'w') as f:
                f.write("# Post-processing results for layer %d\n" % layer_index)
                f.write("# Location\tNodeLabel\tElementLabel\tX\tY\tZ\tS11\tS22\tS33\tS12\tS13\tS23\n")
                
                for result in results:
                    f.write("%s\t%d\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % (
                        result['location'],
                        result['nodeLabel'],
                        result['elementLabel'],
                        result['coordinates'][0],
                        result['coordinates'][1], 
                        result['coordinates'][2],
                        result['stress'][0],  # S11
                        result['stress'][1],  # S22
                        result['stress'][2],  # S33
                        result['stress'][3],  # S12
                        result['stress'][4],  # S13
                        result['stress'][5]   # S23
                    ))
            
            print("  ✓ Results written to %s" % output_filename)
            return True
            
        else:
            print("  ✗ No results found for layer %d" % layer_index)
            return False
            
    except Exception as e:
        print("  ✗ Error extracting stress results: %s" % str(e))
        return False

def extract_displacement_results(odbObject, stepName, extraction_points, output_filename):
    """
    Extract displacement results at critical locations
    """
    try:
        lastFrame = odbObject.steps[stepName].frames[-1]
        dispField = lastFrame.fieldOutputs['U']
        coordField = lastFrame.fieldOutputs['COORD']
        
        results = []
        
        for point_name, (x_target, y_target, z_target, tolerance) in extraction_points.items():
            print("  Extracting displacements near point %s" % point_name)
            
            for i, disp_value in enumerate(dispField.values):
                try:
                    coord_value = coordField.values[i]
                    x, y, z = coord_value.data
                    
                    if (abs(x - x_target) < tolerance and 
                        abs(y - y_target) < tolerance and 
                        abs(z - z_target) < tolerance):
                        
                        disp_data = disp_value.data
                        results.append({
                            'nodeLabel': disp_value.nodeLabel,
                            'coordinates': (x, y, z),
                            'displacement': disp_data,
                            'location': point_name
                        })
                        
                except Exception:
                    continue
        
        if results:
            with open(output_filename, 'w') as f:
                f.write("# Displacement results\n")
                f.write("# Location\tNodeLabel\tX\tY\tZ\tU1\tU2\tU3\n")
                
                for result in results:
                    f.write("%s\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % (
                        result['location'],
                        result['nodeLabel'],
                        result['coordinates'][0],
                        result['coordinates'][1],
                        result['coordinates'][2],
                        result['displacement'][0],  # U1
                        result['displacement'][1],  # U2
                        result['displacement'][2]   # U3
                    ))
            
            print("  ✓ Displacement results written to %s" % output_filename)
            return True
        else:
            print("  ✗ No displacement results found")
            return False
            
    except Exception as e:
        print("  ✗ Error extracting displacement results: %s" % str(e))
        return False

def postprocess_cylindrical_spheroid(modelName, analysis_newPath, rk, bk, height, thetaEnd, N):
    """
    Comprehensive post-processing for cylindrical spheroid geometry with dual coordinate systems
    Adapts 2D_L_CFRP_FEM_circular_path.py methodology for cylinder+hemisphere geometry
    """
    try:
        import time
        startTime = time.time()

        # Define interface evaluation indices
        iInterfaceEval = list(range(1, N))  # Interface indices between layers

        print("\n" + "="*70)
        print("STARTING DUAL COORDINATE SYSTEM POST-PROCESSING")
        print("="*70)

        # Open ODB file
        odbPath = analysis_newPath + '\\' + modelName + '_Job.odb'
        print("Opening ODB file: %s" % odbPath)

        try:
            odbObject = session.openOdb(name=odbPath)
            print("✓ ODB file opened successfully")

            # Get the actual step name from ODB
            stepName_actual = list(odbObject.steps.keys())[0]  # Use actual step name from ODB
            print("Using step: %s" % stepName_actual)
            stepName = stepName_actual  # Use the actual step name from ODB

        except Exception as e:
            print("✗ Error opening ODB file: %s" % str(e))
            return False

        # Setup dual coordinate systems for post-processing (following 2D methodology)
        print("\nSetting up dual coordinate systems...")
        postProc_Cyl_CS, postProc_Sph_CS = setup_dual_postprocessing_coordinate_systems(odbObject, rk, height)

        if postProc_Cyl_CS is None or postProc_Sph_CS is None:
            print("✗ Failed to setup coordinate systems")
            odbObject.close()
            return False

        # Classify elements by geometry region (cylindrical shell vs hemispherical dome)
        print("\nClassifying elements by geometry region...")

        # We'll classify elements within each FEInner set, not globally
        # This avoids the 'ALL ELEMENTS' issue and focuses on the sets we actually use
        cylindrical_elements = []
        spherical_elements = []

        print("Available element sets: %s" % list(odbObject.rootAssembly.elementSets.keys()))

        print("✓ Element classification will be done per FEInner set during processing")

        print("\nExtracting stress results using FEInner sets with dual coordinate systems...")

        # Define geometry-adaptive evaluation parameters for cylindrical pressure vessel
        # Use fixed θ and z evaluation points (NOT radial evaluation)
        thetaEval = angleEval  # Mid-circumferential angle for TS/TG classification
        zEval = 0.0  # Axial evaluation point for ZS/ZG classification (at mid-height)

        print("Cylindrical pressure vessel evaluation points: theta_eval = %.3f rad, z_eval = %.3f (TS/TG + ZS/ZG methodology)" % (thetaEval, zEval))

        # Process each layer using the FEInner sets - adapted from 3D_plate_v1 methodology
        for layerIndex in range(N):  # 0-indexed like 3D_plate_v1
            setName = "FEInner{}".format(layerIndex + 1)
            print("\nProcessing layer %d (%s)..." % (layerIndex + 1, setName))

            try:
                elementSet = odbObject.rootAssembly.elementSets[setName.upper()]
                lastFrame = odbObject.steps[stepName].frames[-1]

                # Classify elements within this set into cylindrical vs spherical regions
                print("  Classifying elements in %s..." % setName)
                cylindrical_elements, spherical_elements = classify_elements_within_set(
                    elementSet, odbObject, stepName
                )
                print("    Cylindrical: %d, Spherical: %d elements" % (len(cylindrical_elements), len(spherical_elements)))

                # Get ELEMENT_NODAL stress field
                stressField = lastFrame.fieldOutputs['S'].getSubset(region=elementSet, position=ELEMENT_NODAL)

                # Convert generator to list to avoid "expecting a recognized type" error
                stressValues = list(stressField.values)

                # Classify elements using TS/TG and ZS/ZG (cylindrical pressure vessel methodology)
                elements_dict, elementLabels_dict = cylindrical_spherical_elemDefLocCenGrav_TS_TG_ZS_ZG(
                    thetaEval, zEval, stressValues, elementSet
                )

                # DEBUG: Show detailed element and node analysis
                print("      DEBUG: ===== DETAILED ELEMENT AND NODE ANALYSIS =====")

                # Find instance for node analysis
                instance = None
                for stress_value in stressValues:
                    instance = stress_value.instance
                    break

                # Build a mapping from (elementLabel, nodeLabel) to stress data for fast lookup
                stress_map = {}
                for s in stressValues:
                    if hasattr(s, 'elementLabel') and hasattr(s, 'nodeLabel'):
                        stress_map[(s.elementLabel, s.nodeLabel)] = s.data

                if instance:
                    # Show elements in each group
                    for group in ['TS_ZS', 'TS_ZG', 'TG_ZS', 'TG_ZG']:
                        group_elements = elementLabels_dict[group]
                        print("      DEBUG: Group %s has elements: %s" % (group, group_elements))

                        # Show all 8 nodes for each element in this group
                        for elemLabel in group_elements[:2]:  # Show first 2 elements to avoid too much output
                            try:
                                element = instance.getElementFromLabel(elemLabel)
                                print("      DEBUG:   Element %d - 8 nodes:" % elemLabel)
                                for i, nodeLabel in enumerate(element.connectivity):
                                    try:
                                        node = instance.getNodeFromLabel(nodeLabel)
                                        x, y, z = node.coordinates[0], node.coordinates[1], node.coordinates[2]
                                        r = (x**2 + z**2)**0.5
                                        if r > 1e-12:
                                            theta = math.atan2(z, x)
                                            if theta < 0:
                                                theta += 2 * math.pi
                                        else:
                                            theta = math.pi/2
                                        theta_diff = abs(theta - thetaEval)
                                        z_diff = abs(y - zEval)
                                        on_eval_line = (theta_diff < 0.05 and z_diff < 0.05)
                                        gui_marker = " [GUI_NODE]" if (abs(x) < 0.001 and abs(y) < 0.001) else ""
                                        eval_marker = " ✓EVAL" if on_eval_line else " ✗eval"
                                        # Get stress data for this node in this element, if available
                                        stress_data = stress_map.get((elemLabel, nodeLabel), None)
                                        if stress_data is not None:
                                            stress_str = " | Stress: [%s]" % (
                                                ", ".join("%.3e" % v for v in stress_data)
                                            )
                                        else:
                                            stress_str = " | Stress: [None]"
                                        print("      DEBUG:     Node %d [%d]: (%.3f,%.3f,%.3f) -> θ=%.3f%s%s%s" % 
                                              (i+1, nodeLabel, x, y, z, theta, eval_marker, gui_marker, stress_str))
                                    except Exception as e:
                                        print("      DEBUG:     Node %d [%d]: ERROR - %s" % (i+1, nodeLabel, str(e)))
                            except Exception as e:
                                print("      DEBUG:   Element %d: ERROR - %s" % (elemLabel, str(e)))

                print("      DEBUG: ===== END DETAILED ANALYSIS =====")

                # Process each classification group (cylindrical pressure vessel methodology)
                for group in ['TS_ZS', 'TS_ZG', 'TG_ZS', 'TG_ZG']:
                    success = process_feinner_group_with_dual_coordinates(
                        odbObject, stepName, postProc_Cyl_CS, postProc_Sph_CS,
                        cylindrical_elements, spherical_elements, setName, group, 
                        elements_dict, stressField, thetaEval, zEval
                    )
                    print("  Group %s processing success: %s" % (group, success))
                    
            except Exception as e:
                print("  ✗ Error processing layer %d: %s" % (layerIndex + 1, str(e)))
        
        # Extract spherical stress results for dome region
        print("Processing completed using integrated dual coordinate methodology")
        
        # Process RG/RS and ZS/ZG classification using integrated approach
        print("\nUsing integrated RG/RS and ZS/ZG methodology with dual coordinates...")
        print("Processing completed - results written to individual group files")
        
        # Process FEInterface sets for interface analysis
        # print("\nProcessing FEInterface sets for interface analysis...")
        # try:
        #     for interface_idx in iInterfaceEval:
        #         interface_set_name = 'FEInterface%d' % interface_idx
        #         try:
        #             interface_set = odbObject.rootAssembly.elementSets[interface_set_name]
        #             print("  Found %s with %d elements" % (interface_set_name, len(interface_set.elements[0])))
        #             # Extract interface stress results here if needed
        #         except KeyError:
        #             print("  ⚠ Interface set '%s' not found" % interface_set_name)
        # except Exception as e:
        #     print("  ✗ Error processing FEInterface sets: %s" % str(e))
        
        # Close ODB
        odbObject.close()
        
        elapsed_time = time.time() - startTime
        print("\n" + "="*70)
        print("FE INNER SET-BASED POST-PROCESSING COMPLETED (CYLINDRICAL PRESSURE VESSEL METHODOLOGY)")
        print("="*70)
        print("Processing time: %.2f seconds" % elapsed_time)
        return True

    except Exception as e:
        print("✗ Error in dual coordinate post-processing: %s" % str(e))
        return False

#Execute post-processing if job was successful
if job_success:
    print("\nJob completed successfully. Starting post-processing...")
    postprocess_success = postprocess_cylindrical_spheroid(
        modelName=modelName,
        analysis_newPath=analysis_newPath, 
        rk=rk,
        bk=bk,
        height=height,
        thetaEnd=thetaEnd,
        N=N
    )
    print("Post-processing success: %s" % postprocess_success)
else:
    print("Job failed. Skipping post-processing.")
