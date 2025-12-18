# -*- coding: utf-8 -*-
"""
@author: anair
"""

import __main__
import os
import time
import datetime as dt
import operator

from abaqus import *
from abaqusConstants import *
#session.setValues(kernelMemoryLimit=16384)
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
import csv


# path_modules = 'U:\\Sachdeva\\MT_Nair\\ABAQUS\\MT\\Macros'
# # path_modules = 'D:\\psingh\\MT\\ABAQUS\\MT\\Macros'
# # path_modules = r"C:\Users\lenovo\Desktop\Aerospace\Thesis\ABAQUS\MT\Macros"
# if path_modules not in sys.path:
#     sys.path.append(path_modules)



session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


##############################   MODEL CREATION   #############################

# Create a new model
Mdb()
modelName = 'EllipseModel_2D'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
mdb.models[modelName].rootAssembly.clearGeometryCache()
mdb.models[modelName].rootAssembly.regenerate()
myModel = mdb.models[modelName]


#############################   PARAMETERS    #############################

r_inner = 200     # semi-axis in a
z_inner = 100     # semi-axis in c
thick   = 1

r_outer = r_inner + thick
z_outer = z_inner + thick

n = 1.0         # superellipse exponent

n_spline = 80  # number of spline points
N_theta = 6   # number of meridional regions
N_part = 4    # Number of partitions
plyAngle = [0, 45, 60, 90]

mesh_size = 0.5  # Mesh size

Press = 1 # Pressure load
compositeMaterialName = 'cfk'  # 'cfk', 'AL', 'gfk', 'cfknew'

##############################   Geometry   #############################

s = myModel.ConstrainedSketch(
    name='__profile__',
    sheetSize=300.0
)

s.sketchOptions.setValues(viewStyle=AXISYM)
s.setPrimaryObject(option=STANDALONE)

# Axis of symmetry
s.ConstructionLine(point1=(0.0, -200.0), point2=(0.0, 200.0))

def superellipse_rz(theta, a, c, n):
    r = a * (math.cos(theta) ** (n))
    z = c * (math.sin(theta) ** (n))
    return r, z


theta_vals = np.linspace(0.0, math.pi/2.0, n_spline)

inner_pts = []
outer_pts = []

for th in theta_vals:

    inner_pts.append(
        superellipse_rz(th, r_inner, z_inner, n)
    )
    outer_pts.append(
        superellipse_rz(th, r_outer, z_outer, n)
    )

# Inner profile
s.Spline(points=inner_pts)

# Outer profile
s.Spline(points=outer_pts[::-1])

# Close bottom (z = 0)
s.Line(point1=inner_pts[0], point2=outer_pts[0])

# Close top (at max z)
s.Line(point1=inner_pts[-1], point2=outer_pts[-1])

s.unsetPrimaryObject()

p = myModel.Part(
    name='SuperEllipsoid_2D',
    dimensionality=AXISYMMETRIC,
    type=DEFORMABLE_BODY,
    twist=ON
)
p.BaseShell(sketch=s)

#################### Face Partitioning ##################

theta_pick = math.radians(10.0)

r_pick, z_pick = superellipse_rz(
    theta_pick, r_inner + 0.5*thick, z_inner + 0.5*thick, n)

f = p.faces
pickedFaces = f.findAt(((r_pick, z_pick, 0.0),))

theta_vals = np.linspace(0.0, math.pi/2.0, 30)

for i in range(1, N_part):

    # --- current offset ---
    a_i = r_inner + i * thick / N_part
    c_i = z_inner + i * thick / N_part
    # --- create sketch for this partition ---
    s1 = myModel.ConstrainedSketch(
        name='partition_%d' % i,
        sheetSize=300.0
    )
    s1.sketchOptions.setValues(viewStyle=AXISYM)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    # axis of symmetry
    s1.ConstructionLine(point1=(0.0, -300.0), point2=(0.0, 300.0))
    # --- generate open spline ---
    part_pts = [
        superellipse_rz(th, a_i, c_i, n)
        for th in theta_vals
    ]
    s1.Spline(points=part_pts)
    # --- apply partition ---
    p.PartitionFaceBySketch(
        faces=pickedFaces,
        sketch=s1
    )
    s1.unsetPrimaryObject()

def surface_normal_rz(r, z, a, c, n, eps=1e-12):
    ra = max(r / a, eps)
    zc = max(z / c, eps)
    e = 2/n
    dr = (e / a) * (ra ** (e - 1.0))
    dz = (e / c) * (zc ** (e - 1.0))
    mag = math.sqrt(dr*dr + dz*dz)
    return (dr / mag, dz / mag)   # (n_r, n_z)

delta = 2.0 * thick

theta_vals = np.linspace(0.0, math.pi/2.0, N_theta + 1)

theta_cut_vals = [
    0.5 * (theta_vals[i] + theta_vals[i+1])
    for i in range(len(theta_vals) - 1)
]

for i, theta in enumerate(theta_cut_vals):
    # --- point on mid-thickness surface ---
    r0, z0 = superellipse_rz(
        theta,
        r_inner,
        z_inner,
        n
    )
    # --- surface normal at that point ---
    nr, nz = surface_normal_rz(
        r0, z0,
        r_inner, z_inner,
        n
    )
    # -- points along the normal ---
    p_out = (r0 + delta * nr, z0 + delta * nz, 0.0)
    p_mid = (r0, z0, 0.0)
    p_in  = (r0 - delta * nr, z0 - delta * nz, 0.0)
    # --- datum points ---
    d1 = p.DatumPointByCoordinate(coords=p_out)
    d2 = p.DatumPointByCoordinate(coords=p_in)
    # d3 = p.DatumPointByCoordinate(coords=p_mid)
    # --- shortest-path partition ---
    faces = p.faces.getByBoundingBox(
    xMin=-1e6, yMin=-1e6, zMin=-1e6,
    xMax=+1e6, yMax=+1e6, zMax=+1e6)
    p.PartitionFaceByShortestPath(
        faces=faces,
        point1=p.datums[d1.id],
        point2=p.datums[d2.id]
    )
    print("Created normal partition at theta = %.2f deg" %
          math.degrees(theta))

################# SET creation ###################

def unique(seq):
    out, seen = [], set()
    for obj in seq:
        key = obj.index
        if key not in seen:
            seen.add(key)
            out.append(obj)
    return out

e1 = p.edges
theta0_edges = []
for k in range(N_part):
    r_mid = r_inner + (k + 0.5) * thick / float(N_part)
    z_mid = 0.0
    theta0_edges.append(e1.findAt(((r_mid, z_mid, 0.0),)))

theta90_edges = []
for k in range(N_part):
    r_mid = 0.0
    z_mid = z_inner + (k + 0.5) * thick / float(N_part)
    theta90_edges.append(e1.findAt(((r_mid, z_mid, 0.0),)))

theta0_edges = unique(theta0_edges)
theta90_edges = unique(theta90_edges)
p.Set(name='set_bottom', edges=theta0_edges)
p.Set(name='set_top', edges=theta90_edges)

N_sample = 40
theta_vals = np.linspace(0.0, math.pi/2.0, N_sample+2)[1:-1]
eps_theta = 1e-6
f1 = p.faces

for k in range(N_part):
    r_mid = r_inner + (k + 0.5) * thick / N_part
    z_mid = z_inner + (k + 0.5) * thick / N_part
    face_layer = []
    for theta in theta_vals:
        if theta < eps_theta or theta > math.pi/2 - eps_theta:
            continue
        r0, z0 = superellipse_rz(theta, r_mid, z_mid, n)
        pt = (r0, z0, 0.0)
        try:
            seq_face = f1.findAt((pt,),)
            if seq_face:
                face_layer.append(seq_face)
        except:
            pass
    # --- create set for this layer ---
    if face_layer:
        set_name = 'Set-Layer-%d' % k
        p.Set(faces=tuple(face_layer), name=set_name)
    else:
        print("WARNING: No faces found for layer", k)

################# Surface creation ###################
N_sample = 40   # more than enough
theta_vals = np.linspace(0.0, math.pi/2.0, N_sample+2)[1:-1]

s2 = p.edges

search_tol = 1e-3
eps_theta = 1e-6            # avoid theta = 0, pi/2
edge_combined = []
for theta in theta_vals:
    # avoid boundaries
    if theta < eps_theta or theta > math.pi/2 - eps_theta:
        continue
    # point on inner surface
    r0, z0 = superellipse_rz(theta, r_outer, z_outer, n)
    pt = (r0, z0, 0.0)
    found = s2.getClosest(
        coordinates=(pt,),
        searchTolerance=search_tol
    )
    if found:
        face_pt = found[0][1]
        seq_face = s2.findAt((face_pt,),)
        if seq_face:
            edge_combined.append(seq_face)

# --- Combine faces into one surface ---
if edge_combined:
    p.Surface(side1Edges=tuple(edge_combined), name='flux_Load')
else:
    print("No faces found for surface creation.")

edge_combined = []
for theta in theta_vals:
    if theta < eps_theta or theta > math.pi/2 - eps_theta:
        continue
    # point on inner surface
    r0, z0 = superellipse_rz(theta, r_inner, z_inner, n)
    pt = (r0, z0, 0.0)
    found = s2.getClosest(
        coordinates=(pt,),
        searchTolerance=search_tol
    )
    if found:
        face_pt = found[0][1]
        seq_face = s2.findAt((face_pt,),)
        if seq_face:
            edge_combined.append(seq_face)

# --- Combine faces into one surface ---
if edge_combined:
    p.Surface(side1Edges=tuple(edge_combined), name='Surf_Load')
else:
    print("No faces found for surface creation.")

#################### Material properties ####################
mdb.models['EllipseModel_2D'].parts['SuperEllipsoid_2D'].setValues(
    geometryRefinement=FINE)
mdb.models[modelName].Material(name='AL')
mdb.models[modelName].materials['AL'].Elastic(
    temperatureDependency=ON, 
    table=((85700.0, 0.318, 0.0),
           (84500.0, 0.32, 100.0),
           (81200.0, 0.325, 200.0),
           (77400.0, 0.33, 300.0))
)
mdb.models[modelName].materials['AL'].Conductivity(
    temperatureDependency=ON, 
    table=((0.0306, 20.0),
           (0.056, 73.0),
           (0.077, 123.0),
           (0.107, 223.0),
           (0.123, 300.0))
)
mdb.models[modelName].materials['AL'].Expansion(
    table=((1.44e-05, 20.0),
           (1.72e-05, 73.0),
           (1.94e-05, 123.0),
           (2.13e-05, 223.0),
           (2.28e-05, 300.0)),
    zero=293.0,
    temperatureDependency=ON
)

mdb.models[modelName].HomogeneousSolidSection(
    name='Section-AL', 
    material='AL', 
    thickness=None
)

p = mdb.models[modelName].parts['SuperEllipsoid_2D']

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

def materialComposite():
    """
    Defines the material properties for composite materials using CORRECT rotation about r-axis.
    Only creates materials for angles actually used in plyAngle list (optimized).
    
    For axisymmetric pressure vessel with RADIAL PLIES:
    - Direction 1 = radial (r) - THROUGH-THICKNESS
    - Direction 2 = axial (z)  
    - Direction 3 = hoop/circumferential (θ)
    - Fibers lie in θ-z plane (tangent to the surface)
    - Fiber angle φ measured about r-axis (rotation in θ-z plane)
    
    Material coordinate system (NEW mapping for STACK_3):
    - Axis 1 (fiber) = axial/z at 0°, rotates in Y-Z plane
    - Axis 2 (transverse in-plane) = hoop/θ
    - Axis 3 (through-thickness/STACK) = radial/r (unchanged by rotation)
    """
    def transform_compliance_axisymmetric(phi_deg):
        """
        Transform compliance matrix for rotation about LOCAL axis 1 (radial/thickness).
        
        Local coordinate system (Material-Local-CS) per Abaqus requirement:
        - Local 1 = Global X (radial/thickness) - ROTATION AXIS
        - Local 2 = Global Y (axial/fiber at 0°)
        - Local 3 = Global Z (hoop/circumferential) - STACK_3 direction (REQUIRED by Abaqus)
        
        Material properties in local coordinates:
        - Material 1 = through-thickness/radial (E1=10000)
        - Material 2 = fiber direction at 0° (E2=147000)
        - Material 3 = transverse/hoop (E3=10000)
        
        Rotation: About local-1 (radial axis) by angle phi in 2-3 plane (axial-hoop)
        - phi=0°: fiber in local-2 (axial)
        - phi=90°: fiber in local-3 (hoop)
        
        Args:
            phi_deg: Fiber angle in degrees from axial (local-2) toward hoop (local-3)
        
        Returns:
            Sr_phi: Transformed 6x6 compliance matrix in local coordinate system
        """
        phi_rad = math.radians(phi_deg)
        c = math.cos(phi_rad)
        s = math.sin(phi_rad)
        
        # Base compliance matrix in material coordinates
        # Material axes: 1=radial/thickness, 2=fiber/axial, 3=transverse/hoop
        S0 = np.zeros((6, 6))
        S0[0, 0] = 1.0 / E1  # Radial/thickness (through-thickness direction)
        S0[1, 1] = 1.0 / E2  # Fiber direction at 0°
        S0[2, 2] = 1.0 / E3  # Transverse/hoop
        S0[0, 1] = S0[1, 0] = -Nu12 / E1  # Radial-Fiber coupling
        S0[0, 2] = S0[2, 0] = -Nu13 / E1  # Radial-Hoop coupling
        S0[1, 2] = S0[2, 1] = -Nu23 / E2  # Fiber-Hoop coupling
        S0[3, 3] = 1.0 / G23  # Shear in 2-3 plane (fiber-hoop)
        S0[4, 4] = 1.0 / G13  # Shear in 1-3 plane (radial-hoop)
        S0[5, 5] = 1.0 / G12  # Shear in 1-2 plane (radial-fiber)
        
        # Transformation matrix: rotation about axis 1 (radial) in 2-3 plane (axial-hoop)
        # Maps local stresses [σ_1, σ_2, σ_3, τ_23, τ_13, τ_12] to rotated material stresses
        T = np.array([
            [1,        0,        0,        0,        0,        0],         # σ'_1 (radial - unchanged)
            [0,        c*c,      s*s,      2*c*s,    0,        0],         # σ'_2 (fiber direction)
            [0,        s*s,      c*c,     -2*c*s,    0,        0],         # σ'_3 (transverse)
            [0,       -c*s,      c*s,   c*c-s*s,    0,        0],         # τ'_23 (fiber-transverse)
            [0,        0,        0,        0,        c,       -s],         # τ'_13 (radial-transverse)
            [0,        0,        0,        0,        s,        c]          # τ'_12 (radial-fiber)
        ])
        
        # Transform: Sr = T^T * S0 * T
        Sr_phi = np.dot(np.dot(T.T, S0), T)
        
        return Sr_phi
    
    # Create materials for each unique ply angle using transformation
    unique_angles = set(plyAngle)
    print("Creating {} pre-rotated materials in local coordinate system...".format(len(unique_angles)))
    
    for phi in sorted(unique_angles):
        # Get transformed compliance matrix for this angle
        Sr_phi = transform_compliance_axisymmetric(phi)
        
        # Extract engineering constants from transformed compliance
        E1_rot = 1.0 / Sr_phi[0, 0]  # Local-1 (radial/thickness direction)
        E2_rot = 1.0 / Sr_phi[1, 1]  # Local-2 (fiber direction after rotation)
        E3_rot = 1.0 / Sr_phi[2, 2]  # Local-3 (hoop direction - STACK_3)
        
        Nu12_rot = -Sr_phi[1, 0] / Sr_phi[0, 0]
        Nu13_rot = -Sr_phi[2, 0] / Sr_phi[0, 0]
        Nu23_rot = -Sr_phi[2, 1] / Sr_phi[1, 1]
        
        G12_rot = 1.0 / Sr_phi[5, 5]
        G13_rot = 1.0 / Sr_phi[4, 4]
        G23_rot = 1.0 / Sr_phi[3, 3]
        
        print("  Angle {}°: E_local1={:.1f}, E_local2={:.1f}, E_local3={:.1f}".format(
            phi, E1_rot, E2_rot, E3_rot))
        
        # Create material with pre-rotated properties in local coordinate system
        material_name = "{}_{:.0f}".format(compositeMaterialName, phi)
        compositeMaterial = myModel.Material(material_name)
        compositeMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=(
            (E1_rot, E2_rot, E3_rot,
             Nu12_rot, Nu13_rot, Nu23_rot,
             G12_rot, G13_rot, G23_rot),))
        
        # Transform thermal expansion for this angle
        # Rotation in 2-3 plane (axial-hoop), axis 1 (radial) unchanged
        c = math.cos(math.radians(phi))
        s = math.sin(math.radians(phi))
        alpha_1_rot = alpha33  # Radial unchanged (axis 1)
        alpha_2_rot = alpha11 * c*c + alpha22 * s*s  # Fiber direction (rotates in 2-3 plane)
        alpha_3_rot = alpha11 * s*s + alpha22 * c*c  # Hoop direction
        
        compositeMaterial.Expansion(type=ORTHOTROPIC, table=((alpha_1_rot, alpha_2_rot, alpha_3_rot),))
        
        # Create section for this material
        section_name = "Composite_Section_{:.0f}".format(phi)
        myModel.HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)

if compositeMaterialName == 'AL':

    faces = p.faces.getByBoundingBox(
        xMin=-1e6, yMin=-1e6, zMin=-1e6,
        xMax= 1e6, yMax= 1e6, zMax= 1e6
    )
    region = regionToolset.Region(faces=faces)
    p.SectionAssignment(
        region=region,
        sectionName='Section-AL',
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField='',
        thicknessAssignment=FROM_SECTION
    )
else:
    materialComposite()
    # ---- Composite: ply-wise assignment ----
    for k, phi in enumerate(plyAngle):
        set_name = 'Set-Layer-%d' % k
        section_name = 'Composite_Section_%d' % phi
        if set_name not in p.sets:
            print("WARNING: {} not found, skipping".format(set_name))
            continue
        if section_name not in myModel.sections:
            print("ERROR: {} does not exist".format(section_name))
            continue
        region = p.sets[set_name]
        p.SectionAssignment(
            region=region,
            sectionName=section_name,
            offset=0.0,
            offsetType=MIDDLE_SURFACE,
            offsetField='',
            thicknessAssignment=FROM_SECTION
        )


################# ASSEMBLY ######################
a = mdb.models[modelName].rootAssembly
a.DatumCsysByThreePoints(coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), 
    point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))
p = mdb.models[modelName].parts['SuperEllipsoid_2D']
a.Instance(name='SuperEllipsoid_2D-1', part=p, dependent=ON)

################ Step Creation ###################
# mdb.models[modelName].CoupledTempDisplacementStep(name='LoadingStep', 
#     previous='Initial', description='LoadingStep', response=STEADY_STATE, 
#     deltmx=None, cetol=None, creepIntegration=None, amplitude=RAMP)

mdb.models[modelName].StaticStep(name='LoadingStep', previous='Initial', 
    initialInc=1, minInc=1e-05, maxInc=1.0)

################# Loading ###################
a1 = mdb.models['EllipseModel_2D'].rootAssembly
region = a1.instances['SuperEllipsoid_2D-1'].surfaces['Surf_Load']
mdb.models['EllipseModel_2D'].Pressure(name='Press_load', 
    createStepName='LoadingStep', region=region, distributionType=UNIFORM, 
    field='', magnitude=Press, amplitude=UNSET)

region = a1.instances['SuperEllipsoid_2D-1'].sets['set_bottom']
mdb.models['EllipseModel_2D'].YsymmBC(name='Bottom', createStepName='Initial', 
    region=region, localCsys=None)

region = a1.instances['SuperEllipsoid_2D-1'].sets['set_top']
mdb.models['EllipseModel_2D'].XsymmBC(name='Top', createStepName='Initial', 
    region=region, localCsys=None)

################# MESHING ###################
p = mdb.models['EllipseModel_2D'].parts['SuperEllipsoid_2D']
p.seedPart(size=mesh_size, deviationFactor=0.01, minSizeFactor=0.1)
elemType1 = mesh.ElemType(elemCode=CGAX8R, elemLibrary=STANDARD)
faces1 = p.faces.getByBoundingBox(
    xMin=-1e6, yMin=-1e6, zMin=-1e6,
    xMax=+1e6, yMax=+1e6, zMax=+1e6
)
pickedRegions = (faces1,)
p.setElementType(regions=pickedRegions, elemTypes=(elemType1,))
p.generateMesh()

############ local coordinate system  #############
def F_family(r, z, a0, c0, e, delta):
    a = a0 + delta
    c = c0 + delta
    return (r/a)**e + (z/c)**e - 1.0

def solve_delta_for_point(r, z, a0, c0, thick, n, tol=1e-10, it=60):
    e=2/n
    lo, hi = 0.0, float(thick)
    f_lo = F_family(r, z, a0, c0, e, lo)
    f_hi = F_family(r, z, a0, c0, e, hi)
    if f_lo * f_hi > 0.0:
        return None
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        f_mid = F_family(r, z, a0, c0, e, mid)
        if abs(f_mid) < tol:
            return mid
        if f_lo * f_mid <= 0.0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)

def layer_index_from_delta(delta, thick, N_part):
    k = int((delta / thick) * N_part)
    if k < 0: k = 0
    if k >= N_part: k = N_part - 1
    return k

def layer_mid_delta(k, thick, N_part):
    return (k + 0.5) * thick / float(N_part)


# elem = p.elements[0]        # just one element
# elemLabel = elem.label

# node_ids = elem.connectivity
# coords = [p.nodes[i].coordinates for i in node_ids]

# node_ids = elem.connectivity
# coords = [p.nodes[i].coordinates for i in node_ids]

# r = [c[0] for c in coords]   # X = radial (r)
# z = [c[1] for c in coords]   # Y = axial  (z)

# node_ids = elem.connectivity
# coords = [p.nodes[i].coordinates for i in node_ids]

# r = [c[0] for c in coords]   # X = radial (r)
# z = [c[1] for c in coords]   # Y = axial  (z)

# r_c = sum(r) / len(r)
# z_c = sum(z) / len(z)

# delta = solve_delta_for_point(r_c, z_c, r_inner, z_inner, thick, 2/n)
# if delta is None:
#     print("Centroid not between inner/outer surfaces (no root in [0,thick]).")
# else:
#     a_ref = r_inner + delta
#     c_ref = z_inner + delta

# delta = solve_delta_for_point(r_c, z_c, r_inner, z_inner, thick, 2/n)
# k = layer_index_from_delta(delta, thick, N_part)
# delta_ref = layer_mid_delta(k, thick, N_part)

# a_ref = r_inner + delta_ref
# c_ref = z_inner + delta_ref


# # a_ref = r_inner + 0.5*thick
# # c_ref = z_inner + 0.5*thick
# n_r, n_z = surface_normal_rz(r_c, z_c, a_ref, c_ref, n)

# # tangent = 90° rotation in r–z plane
# t_r = -n_z
# t_z =  n_r

# dcsys = p.DatumCsysByThreePoints(
#     name='CSYS_try_%d' % elemLabel,
#     coordSysType=CARTESIAN,
#     origin=(r_c, z_c, 0.0),
#     point1=(r_c + t_r, z_c + t_z, 0.0),   # X = tangent
#     point2=(r_c + n_r, z_c + n_z, 0.0)    # Y = normal
# )

# elemSet = p.SetFromElementLabels(
#     name='ONE_try',
#     elementLabels=(elemLabel,)
# )

# p.MaterialOrientation(
#     region=elemSet,
#     orientationType=SYSTEM,
#     localCsys=p.datums[dcsys.id],
#     axis=AXIS_3,
#     angle=90.0,
#     additionalRotationType=ROTATION_ANGLE,
#     stackDirection=STACK_3
# )

for elem in p.elements:
    elemLabel = elem.label
    # --- element centroid ---
    node_ids = elem.connectivity
    coords = [p.nodes[i].coordinates for i in node_ids]
    r = [c[0] for c in coords]
    z = [c[1] for c in coords]
    r_c = sum(r) / len(r)
    z_c = sum(z) / len(z)
    # --- solve for thickness offset ---
    delta = solve_delta_for_point(
        r_c, z_c,
        r_inner, z_inner,
        thick, n
    )
    if delta is None:
        continue   # skip elements outside the wall
    # --- map to thickness layer ---
    k = layer_index_from_delta(delta, thick, N_part)
    delta_ref = layer_mid_delta(k, thick, N_part)
    a_ref = r_inner + delta_ref
    c_ref = z_inner + delta_ref
    # --- normal + tangent ---
    n_r, n_z = surface_normal_rz(r_c, z_c, a_ref, c_ref, n)
    t_r, t_z = -n_z, n_r
    # --- local CSYS ---
    dcsys = p.DatumCsysByThreePoints(
        name='CSYS_ELEM_%d' % elemLabel,
        coordSysType=CARTESIAN,
        origin=(r_c, z_c, 0.0),
        point1=(r_c + t_r, z_c + t_z, 0.0),   # X = tangent
        point2=(r_c + n_r, z_c + n_z, 0.0)    # Y = normal
    )
    # --- element set ---
    set_name = 'ELEM_SET_%d' % elemLabel
    elemSet = p.SetFromElementLabels(
        name=set_name,
        elementLabels=(elemLabel,)
    )
    # --- assign material orientation ---
    p.MaterialOrientation(
        region=elemSet,
        orientationType=SYSTEM,
        localCsys=p.datums[dcsys.id],
        axis=AXIS_3,
        angle=90.0,
        additionalRotationType=ROTATION_ANGLE,
        stackDirection=STACK_3
    )

session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    datumCoordSystems=OFF)

################# JOB Creation ###################
job = mdb.Job(name='Job-1', model='EllipseModel_2D', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='C:\\abaqus_tmp', resultsFormat=ODB, numThreadsPerMpiProcess=0, numCpus=6, 
        numDomains=6, numGPUs=0)
job.submit(consistencyChecking=OFF)
job.waitForCompletion()