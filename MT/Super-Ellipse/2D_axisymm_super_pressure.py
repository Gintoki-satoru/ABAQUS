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
from __future__ import print_function
from abaqusConstants import ELEMENT_NODAL
from odbAccess import openOdb
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

r_inner = 103.59     # semi-axis in a
z_inner = 932.33      # semi-axis in c

n = 0.65         # superellipse exponent

n_spline = 150  # number of spline points
N_theta = 10   # number of meridional regions

plyAngle = [60, -30, -60, 30]  # stacking sequence (degrees)
thick   = 0.16*plyAngle.__len__()  # total thickness
N_part = plyAngle.__len__()  # number of partitions through thickness

r_outer = r_inner + thick
z_outer = z_inner + thick

mesh_size = 0.5  # Mesh size

Press = 1 # Pressure load
compositeMaterialName = 'car_epx'  # 'cfk', 'AL', 'gfk', 'cfknew', 'car_epx', 'im7_epx'

# Strength parameters
Xt = 3179.2    # Longitudinal tensile strength
Yt = 55.701      # Transverse tensile strength
Xc = -1705.3    # Longitudinal compressive strength
Yc = -367.44     # Transverse compressive strength
S = 199.11       # Shear strength

# Puck inclination parameters
p12_t = 0.30
p23_t = 0.3
p12_c = 0.35
p23_c = 0.3

# Optional fiber correction for FF (set False if you don't want it)
USE_FIBER_CORR = True
nu12     = 0.44     # lamina major Poisson ratio
nu_f12   = 0.26    # fiber major Poisson ratio
E11      = 173700     # lamina E11
Ef11     = 294000.0   # fiber E11
m_sigma_f = 1.1       # mean magnification factor
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

theta_vals = np.linspace(0.0, (math.pi / 2.0), n_spline)

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

theta_vals = np.linspace(0.0, math.pi/2.0, 150)

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
# val = (88 / 90.0) * (math.pi / 2.0)
# theta_cut_vals.append(val)

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
    cfk_table = [
        # (T[K], E1, E2, E3, G12, G13, G23, nu21, nu23, nu13, alpha_fiber(2), alpha_trans(1=3))
        (293.0, 11380.0, 161000.0, 11380.0, 5200.0, 3900.0, 5200.0, 0.32, 0.32, 0.45, -9e-7, 2.88e-5),
    ]
elif compositeMaterialName == 'im7_epx':
    im7_epx_table = [
        # (T[K], E1, E2, E3, G12, G13, G23, nu21, nu23, nu13, alpha_fiber(2), alpha_trans(1=3))
        (293.0, 7370.0, 144990.0, 7370.0, 4920.0, 2729.0, 4920.0, 0.34, 0.34, 0.34, -9e-7, 2.88e-5),
    ]
elif compositeMaterialName == 'car_epx':
    # Temperature-dependent base properties (units: MPa, 1/°C)
    # Table order: (T, E1, E2, E3, G12, G13, G23, nu21, nu23, nu13, alpha_fiber(2), alpha_trans(1=3))
    car_epx_table = [
        # T[K],  E1     E2      E3     G12    G13    G23    nu21   nu23   nu13   a2(×1e-7)   a1(×1e-6)
        (293.0,  8.53e3, 171e3,  8.53e3, 5.63e3, 2.64e3, 5.63e3, 0.287, 0.287, 0.369, -1.06e-7, 25.8e-6),
        (193.0, 10.46e3, 172e3, 10.46e3, 8.14e3, 3.46e3, 8.14e3, 0.328, 0.328, 0.369, -0.45e-7, 19.1e-6),
        (93.0,  13.55e3, 173e3, 13.55e3, 9.31e3, 4.80e3, 9.31e3, 0.389, 0.389, 0.369, -0.51e-7, 11.9e-6),
    ]
    car_epx_k_table = [
    # (T[K], k22_fiber, k11_radial, k33_hoop)
    (292.0, 0.006317, 0.000632, 0.000632),
    (70.0,  0.000857, 0.000234, 0.000234),
    (23.0,  0.000355, 0.000151, 0.000151),
    ]
else:
	pass

if compositeMaterialName == 'car_epx':
    prop_table = car_epx_table
elif compositeMaterialName == 'im7_epx':
    prop_table = im7_epx_table
elif compositeMaterialName == 'cfk':
    prop_table = cfk_table
else:
    prop_table = None

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
    def transform_compliance_axisymmetric(phi_deg, E1, E2, E3, Nu12, Nu13, Nu23, G12, G13, G23):
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
        
        return np.dot(np.dot(T.T, S0), T)
    
    # Create materials for each unique ply angle using transformation
    unique_angles = set(plyAngle)
    print("Creating {} pre-rotated materials in local coordinate system...".format(len(unique_angles)))
    
    for phi in sorted(unique_angles):
        material_name = "{}_{:.0f}".format(compositeMaterialName, phi)
        compositeMaterial = myModel.Material(material_name)
        elastic_table = []
        expansion_table = []
        # --- loop over temperature rows (for car_epx) ---
        for (tempK, E1t, E2t, E3t, G12t, G13t, G23t, nu21t, nu23t, nu13t, alpha_fiber, alpha_trans) in prop_table:
            # reciprocity: nu12 from nu21
            nu12t = nu21t * (E1t / E2t)
            # Get transformed compliance for this angle + temperature
            Sr_phi = transform_compliance_axisymmetric(
                phi, E1t, E2t, E3t, nu12t, nu13t, nu23t, G12t, G13t, G23t
            )
            # Extract engineering constants from transformed compliance
            E1_rot = 1.0 / Sr_phi[0, 0]
            E2_rot = 1.0 / Sr_phi[1, 1]
            E3_rot = 1.0 / Sr_phi[2, 2]
            Nu12_rot = -Sr_phi[1, 0] / Sr_phi[0, 0]
            Nu13_rot = -Sr_phi[2, 0] / Sr_phi[0, 0]
            Nu23_rot = -Sr_phi[2, 1] / Sr_phi[1, 1]
            G12_rot = 1.0 / Sr_phi[5, 5]
            G13_rot = 1.0 / Sr_phi[4, 4]
            G23_rot = 1.0 / Sr_phi[3, 3]
            # Temperature-dependent elastic constants: temperature is LAST entry
            elastic_table.append(
                (E1_rot, E2_rot, E3_rot,
                Nu12_rot, Nu13_rot, Nu23_rot,
                G12_rot, G13_rot, G23_rot,
                tempK)
            )
            # ---- thermal expansion: map base values (local axes) then rotate ----
            # base: alpha11=alpha33=alpha_trans, alpha22=alpha_fiber
            alpha11_t = alpha_trans
            alpha22_t = alpha_fiber
            alpha33_t = alpha_trans
            c = math.cos(math.radians(phi))
            s = math.sin(math.radians(phi))
            alpha_1_rot = alpha11_t
            alpha_2_rot = alpha22_t * c*c + alpha33_t * s*s
            alpha_3_rot = alpha22_t * s*s + alpha33_t * c*c
            # Temperature-dependent expansion: temperature is LAST entry
            expansion_table.append((alpha_1_rot, alpha_2_rot, alpha_3_rot, tempK))
        # Create material with temperature-dependent tables
        if compositeMaterialName == 'car_epx':
            compositeMaterial.Elastic(
                type=ENGINEERING_CONSTANTS,
                temperatureDependency=ON,
                table=tuple(elastic_table)
            )
        else:
            # temperature-independent (single-row table assumed)
            compositeMaterial.Elastic(
                type=ENGINEERING_CONSTANTS,
                table=(elastic_table[0][:9],)  # strip temperature column
            )
        if compositeMaterialName == 'car_epx':
            compositeMaterial.Expansion(
                type=ORTHOTROPIC,
                temperatureDependency=ON,
                table=tuple(expansion_table)
            )
        else:
            compositeMaterial.Expansion(
                type=ORTHOTROPIC,
                table=(expansion_table[0][:3],)
            )
        if compositeMaterialName == 'car_epx':
            k_table = []
            for (tempK, k22, k11, k33) in car_epx_k_table:
                # Rotation about axis 1 in 2-3 plane
                c = math.cos(math.radians(phi))
                s = math.sin(math.radians(phi))
                k11_rot = k11
                k22_rot = k22 * c*c + k33 * s*s
                k33_rot = k22 * s*s + k33 * c*c
                # Abaqus ORTHOTROPIC order: (k11, k22, k33, T)
                k_table.append((k11_rot, k22_rot, k33_rot, tempK))
            compositeMaterial.Conductivity(
                type=ORTHOTROPIC,
                temperatureDependency=ON,
                table=tuple(k_table)
            )
        print("  Angle {}°: created temp-dependent material with {} rows".format(phi, len(elastic_table)))
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

r_cutoff = 7.0
z_cutoff = 4.0
theta_min = 0
theta_min = math.radians(theta_min)
exclude_elems = set()
elem_to_phi = {}
elem_to_k   = {}
layer_to_elems = {}

for elem in p.elements:
    elemLabel = elem.label
    # --- element centroid ---
    node_ids = elem.connectivity
    coords = [p.nodes[i].coordinates for i in node_ids]
    r = [c[0] for c in coords]
    z = [c[1] for c in coords]
    r_c = sum(r) / len(r)
    z_c = sum(z) / len(z)
    if r_c < r_cutoff:
        exclude_elems.add(elemLabel)
    if z_c < z_cutoff:
        exclude_elems.add(elemLabel)
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
    if k not in layer_to_elems:
        layer_to_elems[k] = []
    layer_to_elems[k].append(elemLabel)
    phi = plyAngle[k]
    elem_to_k[elemLabel] = k
    elem_to_phi[elemLabel] = phi
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
    angle_deg = 90.0 if elemLabel in exclude_elems else 90.0
    p.MaterialOrientation(
        region=elemSet,
        orientationType=SYSTEM,
        localCsys=p.datums[dcsys.id],
        axis=AXIS_3,
        angle=angle_deg,
        additionalRotationType=ROTATION_ANGLE,
        stackDirection=STACK_3
    )

p.SetFromElementLabels(
    name='EXCLUDE_TS_WU',
    elementLabels=list(exclude_elems)
)

session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    datumCoordSystems=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(
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

################# POST-PROCESSING ###################
odbPath = 'Job-1.odb'

H1  = (1.0 / Xt) - (1.0 / Xc)
H2  = (1.0 / Yt) - (1.0 / Yc)
H11 = 1.0 / (Xt * Xc)
H22 = 1.0 / (Yt * Yc)
H6  = 0.0
H66 = 1.0 / (S * S)
H12 = -0.5 * math.sqrt(H11 * H22)

def tsai_wu_index(s1, s2, t12):
    return (H1*s1 + H2*s2 + H6*t12 +
            H11*s1*s1 + H22*s2*s2 + H66*t12*t12 +
            2.0*H12*s1*s2)

def voigt6_to_tensor(s11, s22, s33, s12, s13, s23):
    return np.array([
        [s11, s12, s13],
        [s12, s22, s23],
        [s13, s23, s33]
    ], dtype=float)

def tensor_to_voigt6(T):
    return (T[0,0], T[1,1], T[2,2], T[0,1], T[0,2], T[1,2])

def R_about_axis1(phi_deg):
    phi = math.radians(phi_deg)
    c = math.cos(phi)
    s = math.sin(phi)
    # Basis (1,2,3) -> (1,2',3')
    return np.array([
        [1.0, 0.0, 0.0],
        [0.0,  c,   s ],
        [0.0, -s,   c ]
    ], dtype=float)

def rotate_tensor_about1(T_global, phi_deg):
    R = R_about_axis1(phi_deg)
    return R.dot(T_global).dot(R.T)

def to_puck_axes_from_rotated(T_rot):
    P = np.array([
        [0.0, 1.0, 0.0],  # puck-1 (fiber)     <- rot-2'
        [0.0, 0.0, 1.0],  # puck-2 (trans)     <- rot-3'
        [1.0, 0.0, 0.0]   # puck-3 (thickness) <- rot-1
    ], dtype=float)
    return P.dot(T_rot).dot(P.T)

def global_to_puck_voigt(s11,s22,s33,s12,s13,s23, phi_deg):
    Tg = voigt6_to_tensor(s11,s22,s33,s12,s13,s23)
    Trot = rotate_tensor_about1(Tg, phi_deg)    # (1,2',3')
    Tp = to_puck_axes_from_rotated(Trot)        # (fiber,trans,thk)
    return tensor_to_voigt6(Tp)


# Fracture plane search
theta_min_deg = -90
theta_max_deg = 90
theta_step_deg = 1

try:
    exclude_elems
except NameError:
    exclude_elems = set()

EPS = 1e-12

# -------------------------
# PUCK FUNCTIONS
# -------------------------
def puck_ff(s11, s22, s33):
    Xc_pos = -Xc
    if USE_FIBER_CORR:
        coeff = (nu12 - nu_f12 * m_sigma_f * (E11 / max(Ef11, EPS)))
        s11_eff = s11 - coeff * (s22 + s33)
    else:
        s11_eff = s11
    if s11_eff >= 0.0:
        f_t = s11_eff / max(Xt, EPS)
        f_c = 0.0
    else:
        f_c = abs(s11_eff) / max(Xc_pos, EPS)
        f_t = 0.0
    return max(f_t, f_c), f_t, f_c, s11_eff


def puck_iff(s22, s33, s12, s13, s23):
    Yc_pos = -Yc
    R23A = Yc_pos / max(2.0 * (1.0 + p23_c), EPS)
    fmax = -1.0
    theta_crit = 0.0
    branch_crit = "tension"
    n_steps = int(round((theta_max_deg - theta_min_deg) / float(theta_step_deg))) + 1
    for k in range(n_steps):
        theta_deg = theta_min_deg + k * theta_step_deg
        th = math.radians(theta_deg)
        c = math.cos(th)
        s = math.sin(th)
        sigma_n = s22 * c*c + s33 * s*s + 2.0 * s23 * s * c
        tau_nt  = (s33 - s22) * s * c + 2.0 * s23 * (c*c - s*s)
        tau_n1  = s13 * s + s12 * c
        denom = tau_nt*tau_nt + tau_n1*tau_n1
        if denom < EPS:
            cos2 = 0.0
            sin2 = 0.0
        else:
            cos2 = (tau_nt*tau_nt) / denom
            sin2 = (tau_n1*tau_n1) / denom
        k_t = (p23_t / max(R23A, EPS)) * cos2 + (p12_t / max(S, EPS)) * sin2
        k_c = (p23_c / max(R23A, EPS)) * cos2 + (p12_c / max(S, EPS)) * sin2
        if sigma_n >= 0.0:
            # matrix tension branch
            term_n = (1.0 / max(Yt, EPS) - k_t) * sigma_n
            f_here = math.sqrt(
                term_n*term_n +
                (tau_nt / max(R23A, EPS))**2 +
                (tau_n1 / max(S, EPS))**2
            ) + k_t * sigma_n
            branch_here = "tension"
        else:
            # matrix compression branch
            f_here = math.sqrt(
                (k_c * sigma_n)**2 +
                (tau_nt / max(R23A, EPS))**2 +
                (tau_n1 / max(S, EPS))**2
            ) - k_c * sigma_n
            branch_here = "compression"
        if f_here > fmax:
            fmax = f_here
            theta_crit = theta_deg
            branch_crit = branch_here
    return fmax, branch_crit, theta_crit

def failure_analysis(elem_to_phi, strain_limit=0.005):
    odb = openOdb(odbPath, readOnly=True)
    stepName = list(odb.steps.keys())[-1]  # safe for both
    step = odb.steps[stepName]
    frameIndex = -1
    frame = step.frames[frameIndex]
    if 'S' not in frame.fieldOutputs:
        raise RuntimeError("Stress output 'S' not found in ODB frame. Request stresses (S).")
    if 'E' not in frame.fieldOutputs:
        raise RuntimeError("Strain output 'E' not found in ODB frame. Request strains (E).")
    Sfield = frame.fieldOutputs['S'].getSubset(position=INTEGRATION_POINT)
    Efield = frame.fieldOutputs['E'].getSubset(position=INTEGRATION_POINT)
    for sv in Sfield.values:
        ip = getattr(sv, 'integrationPoint', None)
        if ip is None:
            continue
        sd = sv.data
        s12 = float(sd[3]) if len(sd) > 3 else 0.0
    E_by_key = {}
    for ev in Efield.values:
        if ev.elementLabel in exclude_elems:
            continue
        inst_name = ev.instance.name if ev.instance else ''
        ip = getattr(ev, 'integrationPoint', None)
        sp = getattr(ev, 'sectionPoint', None)
        if ip is None:
            continue
        ed = ev.data
        e11 = float(ed[0]) if len(ed) > 0 else 0.0
        e22 = float(ed[1]) if len(ed) > 1 else 0.0
        e33 = float(ed[2]) if len(ed) > 2 else 0.0
        e12 = float(ed[3]) if len(ed) > 3 else 0.0
        e13 = float(ed[4]) if len(ed) > 4 else 0.0
        e23 = float(ed[5]) if len(ed) > 5 else 0.0
        E_by_key[(inst_name, ev.elementLabel, ip, sp)] = (e11, e22, e33, e12, e13, e23)
    base = os.path.splitext(os.path.basename(odbPath))[0]
    out_csv = 'allCriteria_%s_IP_ply3D.csv' % base
    max_any = -1.0;  max_any_elem=None; max_any_ip=None; max_any_inst=''; max_any_crit=''
    max_tw  = -1.0;  max_tw_elem=None;  max_tw_ip=None;  max_tw_inst=''
    max_eps_tension = -1.0; max_eps_elem=None; max_eps_ip=None; max_eps_inst=''
    max_puck = -1.0; max_puck_elem=None; max_puck_ip=None; max_puck_inst=''
    missing_strain = 0
    # --- Python 2/3 compatible file open ---
    if sys.version_info[0] < 3:
        f = open(out_csv, 'wb')
    else:
        f = open(out_csv, 'w', newline='')
    with f:
        w = csv.writer(f)
        w.writerow([
            'odb','step','frameIndex','frameValue',
            'instance','elementLabel','integrationPoint','sectionPoint',
            'S11','S22','S33','S12','S13','S23',
            'E11','E22','E33','E12','E13','E23',
            'phi_deg',
            'S11_p','S22_p','S33_p','S12_p','S13_p','S23_p',
            'E11_p','E22_p','E33_p','E12_p','E13_p','E23_p',
            'TW','FAIL_TW(>=1)',
            'max_tension_strain(max(E11_p+,E22_p+))','strain_limit','FAIL_tension_strain(>limit)',
            'Puck_FF','Puck_FF_t','Puck_FF_c','S_f_eff',
            'Puck_IFF','Puck_IFF_branch','theta_crit_deg',
            'FAIL_FF(>=1)','FAIL_IFF(>=1)',
            'FAIL_any'
        ])
        for sv in Sfield.values:
            el = sv.elementLabel
            if el in exclude_elems:
                continue
            if el not in elem_to_phi:
                continue
            inst_name = sv.instance.name if sv.instance else ''
            ip = getattr(sv, 'integrationPoint', None)
            sp = getattr(sv, 'sectionPoint', None)
            if ip is None:
                continue
            sd = sv.data
            s11 = float(sd[0]) if len(sd) > 0 else 0.0
            s22 = float(sd[1]) if len(sd) > 1 else 0.0
            s33 = float(sd[2]) if len(sd) > 2 else 0.0
            s12 = float(sd[3]) if len(sd) > 3 else 0.0
            s13 = float(sd[4]) if len(sd) > 4 else 0.0
            s23 = float(sd[5]) if len(sd) > 5 else 0.0
            phi = float(elem_to_phi[el])
            s11p, s22p, s33p, s12p, s13p, s23p = global_to_puck_voigt(
                s11,s22,s33,s12,s13,s23, phi)
            key = (inst_name, el, ip, sp)
            if key in E_by_key:
                e11, e22, e33, e12, e13, e23 = E_by_key[key]
                e11p, e22p, e33p, e12p, e13p, e23p = global_to_puck_voigt(
                    e11,e22,e33,e12,e13,e23, phi)
            else:
                missing_strain += 1
                e11=e22=e33=e12=e13=e23=float('nan')
                e11p=e22p=e33p=e12p=e13p=e23p=float('nan')
            tw = tsai_wu_index(s11p, s22p, s12p)
            fail_tw = 1 if tw >= 1.0 else 0
            if tw > max_tw:
                max_tw = tw; max_tw_elem = el; max_tw_ip = ip; max_tw_inst = inst_name
            if not (math.isnan(e11p) or math.isnan(e22p)):
                e11t = e11p if e11p > 0.0 else 0.0
                e22t = e22p if e22p > 0.0 else 0.0
                max_tension = max(e11t, e22t)
                fail_eps = 1 if (e11p > strain_limit or e22p > strain_limit) else 0
                if max_tension > max_eps_tension:
                    max_eps_tension = max_tension
                    max_eps_elem = el; max_eps_ip = ip; max_eps_inst = inst_name
            else:
                max_tension = float('nan'); fail_eps = 0
            ff, ff_t, ff_c, s_f_eff = puck_ff(s11p, s22p, s33p)
            iff, iff_branch, theta_star = puck_iff(s22p, s33p, s12p, s13p, s23p)
            fail_ff = 1 if ff >= 1.0 else 0
            fail_iff = 1 if iff >= 1.0 else 0
            fail_any = 1 if (fail_tw or fail_eps or fail_ff or fail_iff) else 0
            worst = max(
                tw,
                max_tension if not math.isnan(max_tension) else -1.0,
                ff, iff
            )
            if worst > max_any:
                max_any = worst; max_any_elem = el; max_any_ip = ip; max_any_inst = inst_name
                if worst == tw:
                    max_any_crit = 'Tsai-Wu'
                elif worst == ff:
                    max_any_crit = 'Puck_FF'
                elif worst == iff:
                    max_any_crit = 'Puck_IFF'
                else:
                    max_any_crit = 'Strain_limit'
            worst_puck = max(ff, iff)
            if worst_puck > max_puck:
                max_puck = worst_puck
                max_puck_elem = el; max_puck_ip = ip; max_puck_inst = inst_name
            w.writerow([
                base, step.name, frameIndex, frame.frameValue,
                inst_name, el, ip, str(sp) if sp is not None else '',
                s11,s22,s33,s12,s13,s23,
                e11,e22,e33,e12,e13,e23,
                phi,
                s11p,s22p,s33p,s12p,s13p,s23p,
                e11p,e22p,e33p,e12p,e13p,e23p,
                tw, fail_tw,
                max_tension, strain_limit, fail_eps,
                ff, ff_t, ff_c, s_f_eff,
                iff, iff_branch, theta_star,
                fail_ff, fail_iff,
                fail_any
            ])
    odb.close()
    print("Wrote:", out_csv)
    if missing_strain:
        print("WARNING: %d integration-point stress entries had no matching integration-point strain entry." % missing_strain)
    print("Max Tsai–Wu (ply/IP): %.4f" % max_tw)
    print("  Instance      :", max_tw_inst)
    print("  Element label :", max_tw_elem)
    print("  Integration pt:", max_tw_ip)
    print("Max tension max(E11_p+,E22_p+) (ply/IP): %.6f (limit %.6f)" % (max_eps_tension, strain_limit))
    print("  Instance      :", max_eps_inst)
    print("  Element label :", max_eps_elem)
    print("  Integration pt:", max_eps_ip)
    print("Max(Puck FF/IFF) (ply/IP): %.4f" % max_puck)
    print("  Instance      :", max_puck_inst)
    print("  Element label :", max_puck_elem)
    print("  Integration pt:", max_puck_ip)
    print("Max(ANY criterion) (ply/IP): %.4f  [%s]" % (max_any, max_any_crit))
    print("  Instance      :", max_any_inst)
    print("  Element label :", max_any_elem)
    print("  Integration pt:", max_any_ip)

failure_analysis(elem_to_phi, strain_limit=0.005)
