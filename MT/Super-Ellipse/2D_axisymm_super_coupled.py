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


# path_modules = 'U:\\Sachdeva\\MT_Nair\\ABAQUS\\MT\\Super-Ellipse'
path_modules = r"C:\Users\lenovo\Desktop\Aerospace\Thesis\ABAQUS\MT\Super-Ellipse"
if path_modules not in sys.path:
    sys.path.append(path_modules)

import Heat_coeff

from Heat_coeff import (
    superellipsoid_area,
    superellipsoid_volume,
    shape_factor,
    equivalent_heat_coeff
)

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

r_inner = 62.74     # semi-axis in a
z_inner = 2823.48   # semi-axis in c

n = 0.5        # superellipse exponent

n_spline = 150  # number of spline points
N_theta = 10   # number of meridional regions

plyAngle = [90,45,-45,90]  # stacking sequence (degrees)
thick   = 0.16*plyAngle.__len__()  # total thickness
thick = 0.918
N_part = plyAngle.__len__()  # number of partitions through thickness

r_outer = r_inner + thick
z_outer = z_inner + thick

mesh_size = 0.5  # Mesh size

Press = 1 # Pressure load
compositeMaterialName = 'AL'  # 'cfk', 'AL', 'gfk', 'cfknew', 'car_epx', 'im7_epx'

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

t_ins = 16                      # insulation thickness                      
t_outer = 2                     # outer layer thickness
k_liner = 0.000151               # liner thermal conductivity
k_outer = 0.0306                # outer layer thermal conductivity
k_ins = 3.0300e-08              # insulation thermal conductivity
conv_coeff = 10 / 1e6         # convection coefficient W/mm^2K
##############################   Geometry   #############################

s = myModel.ConstrainedSketch(
    name='__profile__',
    sheetSize=300.0
)

s.sketchOptions.setValues(viewStyle=AXISYM)
s.setPrimaryObject(option=STANDALONE)

# Axis of symmetry
s.ConstructionLine(point1=(0, -200.0), point2=(0, 200.0))

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
    s1.ConstructionLine(point1=(0, -300.0), point2=(0, 300.0))
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
    p.Set(edges=tuple(edge_combined), name='Temp_Load')
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
elif compositeMaterialName == 'im7_epx':
    im7_epx_table = [
        # (T[K], E1, E2, E3, G12, G13, G23, nu21, nu23, nu13, alpha_fiber(2), alpha_trans(1=3))
        (293.0, 11380.0, 161000.0, 11380.0, 5200.0, 3900.0, 5200.0, 0.32, 0.32, 0.45, -9e-7, 2.88e-5),
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
else:
    prop_table = None

def materialComposite():
    def transform_compliance_axisymmetric(phi_deg, E1, E2, E3, Nu12, Nu13, Nu23, G12, G13, G23):
        """
        Rotate compliance about LOCAL axis-1 (radial/thickness) by phi in the 2-3 plane (axial-hoop).

        Local axes convention USED HERE:
          1 = radial/thickness
          2 = axial (fiber direction at phi=0)
          3 = hoop  (fiber direction at phi=90)

        Returns:
          Sr_phi : rotated 6x6 compliance matrix in the SAME 1-2-3 basis (Voigt: 11,22,33,23,13,12)
        """
        phi_rad = math.radians(phi_deg)
        c = math.cos(phi_rad)
        s = math.sin(phi_rad)
        # Base compliance in material axes (1=radial, 2=fiber@0, 3=hoop/transverse)
        S0 = np.zeros((6, 6))
        S0[0, 0] = 1.0 / E1
        S0[1, 1] = 1.0 / E2
        S0[2, 2] = 1.0 / E3
        S0[0, 1] = S0[1, 0] = -Nu12 / E1
        S0[0, 2] = S0[2, 0] = -Nu13 / E1
        S0[1, 2] = S0[2, 1] = -Nu23 / E2
        S0[3, 3] = 1.0 / G23  # 23
        S0[4, 4] = 1.0 / G13  # 13
        S0[5, 5] = 1.0 / G12  # 12
        # Stress transformation in Voigt for rotation about axis-1 in 2-3 plane
        # Voigt ordering assumed: [11,22,33,23,13,12]
        T = np.array([
            [1,        0,        0,        0,        0,        0],
            [0,        c*c,      s*s,      2*c*s,    0,        0],
            [0,        s*s,      c*c,     -2*c*s,    0,        0],
            [0,       -c*s,      c*s,   c*c-s*s,    0,        0],
            [0,        0,        0,        0,        c,       -s],
            [0,        0,        0,        0,        s,        c]
        ])
        # Rotate compliance: S' = T^T * S * T
        Sr_phi = np.dot(np.dot(T.T, S0), T)
        return Sr_phi
    def C_to_abaqus_anisotropic_table(C):
        """
        Map 6x6 stiffness matrix C (Voigt: [11,22,33,23,13,12]) to Abaqus
        *ELASTIC, TYPE=ANISOTROPIC data order:

        Line 1: D1111 D1122 D2222 D1133 D2233 D3333 D1112 D2212
        Line 2: D3312 D1212 D1113 D2213 D3313 D1213 D1313 D1123
        Line 3: D2223 D3323 D1223 D1323 D2323
        """
        # Voigt indices: 0->11, 1->22, 2->33, 3->23, 4->13, 5->12
        table = (
            C[0,0], C[0,1], C[1,1], C[0,2], C[1,2], C[2,2], C[0,5], C[1,5],
            C[2,5], C[5,5], C[0,4], C[1,4], C[2,4], C[5,4], C[4,4], C[0,3],
            C[1,3], C[2,3], C[5,3], C[4,3], C[3,3]
        )
        return table
    def assert_stiffness_ok(C, material_name, phi, tempK):
        if not np.allclose(C, C.T, atol=1e-6, rtol=1e-6):
            raise RuntimeError("C not symmetric for {} at phi={} T={}".format(material_name, phi, tempK))
        if np.any(np.diag(C) <= 0.0):
            raise RuntimeError("Non-positive stiffness diagonal for {} at phi={} T={}".format(material_name, phi, tempK))
    unique_angles = sorted(set(plyAngle))
    for phi in unique_angles:
        material_name = "{}_{:.0f}".format(compositeMaterialName, phi)
        compositeMaterial = myModel.Material(material_name)
        anisotropic_table_rows = []
        expansion_table = []
        cphi = math.cos(math.radians(phi))
        sphi = math.sin(math.radians(phi))
        for (tempK, E1t, E2t, E3t, G12t, G13t, G23t, nu21t, nu23t, nu13t, alpha_fiber, alpha_trans) in prop_table:
            # Reciprocity
            nu12t = nu21t * (E1t / E2t)
            Sr_phi = transform_compliance_axisymmetric(
                phi, E1t, E2t, E3t, nu12t, nu13t, nu23t, G12t, G13t, G23t
            )
            C_phi = np.linalg.inv(Sr_phi)
            assert_stiffness_ok(C_phi, material_name, phi, tempK)
            row = C_to_abaqus_anisotropic_table(C_phi)
            if compositeMaterialName == 'car_epx':
                anisotropic_table_rows.append(tuple(list(row) + [tempK]))
            else:
                anisotropic_table_rows.append(row)
            alpha11_t = alpha_trans
            alpha22_t = alpha_fiber
            alpha33_t = alpha_trans
            alpha_1_rot = alpha11_t
            alpha_2_rot = alpha22_t * cphi*cphi + alpha33_t * sphi*sphi
            alpha_3_rot = alpha22_t * sphi*sphi + alpha33_t * cphi*cphi
            if compositeMaterialName == 'car_epx':
                expansion_table.append((alpha_1_rot, alpha_2_rot, alpha_3_rot, tempK))
            else:
                expansion_table.append((alpha_1_rot, alpha_2_rot, alpha_3_rot))
        if compositeMaterialName == 'car_epx':
            compositeMaterial.Elastic(
                type=ANISOTROPIC,
                temperatureDependency=ON,
                table=tuple(anisotropic_table_rows)
            )
            compositeMaterial.Expansion(
                type=ORTHOTROPIC,
                temperatureDependency=ON,
                table=tuple(expansion_table)
            )
        else:
            compositeMaterial.Elastic(
                type=ANISOTROPIC,
                table=(anisotropic_table_rows[0],)
            )
            compositeMaterial.Expansion(
                type=ORTHOTROPIC,
                table=(expansion_table[0],)
            )
        if compositeMaterialName == 'car_epx':
            k_table = []
            for (tempK, k22, k11, k33) in car_epx_k_table:
                k11_rot = k11
                k22_rot = k22 * cphi*cphi + k33 * sphi*sphi
                k33_rot = k22 * sphi*sphi + k33 * cphi*cphi
                k_table.append((k11_rot, k22_rot, k33_rot, tempK))
            compositeMaterial.Conductivity(
                type=ORTHOTROPIC,
                temperatureDependency=ON,
                table=tuple(k_table)
            )
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
mdb.models[modelName].CoupledTempDisplacementStep(name='LoadingStep', 
    previous='Initial', description='LoadingStep', response=STEADY_STATE, 
    deltmx=None, cetol=None, creepIntegration=None, amplitude=RAMP)

mdb.models['EllipseModel_2D'].fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 
    'CDISP', 'NT', 'HFL', 'RFL', 'COORD'))

# mdb.models[modelName].StaticStep(name='LoadingStep', previous='Initial', 
#     initialInc=1, minInc=1e-05, maxInc=1.0)

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
# mdb.models['EllipseModel_2D'].XsymmBC(name='Top', createStepName='Initial', 
#     region=region, localCsys=None)

region = a.instances['SuperEllipsoid_2D-1'].sets['Temp_Load']
mdb.models['EllipseModel_2D'].TemperatureBC(name='Temp', 
    createStepName='LoadingStep', region=region, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', magnitude=20.0, amplitude=UNSET)

A_liner = superellipsoid_area(r_inner, r_inner, z_inner, n, 1)
A_ins = superellipsoid_area(r_inner + thick, r_inner + thick, z_inner + thick, n, 1)
A_outer = superellipsoid_area(r_inner + thick + t_ins, r_inner + thick + t_ins, z_inner + thick + t_ins, n, 1)#
A_outer_total = superellipsoid_area(r_inner + thick + t_ins + t_outer, r_inner + thick + t_ins + t_outer, z_inner + thick + t_ins + t_outer, n, 1)

# 2. Compute volume
V_inner = superellipsoid_volume(r_inner + thick, r_inner + thick, z_inner + thick, n, 1) - superellipsoid_volume(r_inner, r_inner, z_inner, n, 1)
V_ins = superellipsoid_volume(r_inner + thick + t_ins, r_inner + thick + t_ins, z_inner + thick + t_ins, n, 1) - superellipsoid_volume(r_inner + thick, r_inner + thick, z_inner + thick, n, 1)
V_outer = superellipsoid_volume(r_inner + thick + t_ins + t_outer, r_inner + thick + t_ins + t_outer, z_inner + thick + t_ins + t_outer, n, 1) - superellipsoid_volume(r_inner + thick + t_ins, r_inner + thick + t_ins, z_inner + thick + t_ins, n, 1)

# 3. Compute shape factors
S_liner = shape_factor(A_liner, V_inner, thick, r_inner, r_inner, z_inner)
S_ins   = shape_factor(A_ins,   V_ins, t_ins,   r_inner, r_inner, z_inner)
S_outer = shape_factor(A_outer, V_outer, t_outer, r_inner, r_inner, z_inner)

# 4. Compute heat transfer coefficient and temperatures
Q_total, T3, T2, T1, h_eq = equivalent_heat_coeff(
    T_air=300, T_LH2=20,
    A_outer=A_outer_total,
    hc=10,
    S_ins=S_ins, k_ins=k_ins,
    S_liner=S_liner, k_liner=k_liner,
    S_outer=S_outer, k_outer=k_outer,
    A_outer_liner=A_ins
)
region = a1.instances['SuperEllipsoid_2D-1'].surfaces['flux_Load']
mdb.models['EllipseModel_2D'].SurfaceHeatFlux(name='HeatFlux-1', 
    createStepName='LoadingStep', region=region, magnitude=Q_total/A_ins, 
    distributionType=UNIFORM)

################# MESHING ###################
p = mdb.models['EllipseModel_2D'].parts['SuperEllipsoid_2D']
p.seedPart(size=mesh_size, deviationFactor=0.01, minSizeFactor=0.1)
elemType1 = mesh.ElemType(elemCode=CGAX4RT, elemLibrary=STANDARD)
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
    theta = atan2(r_c, z_c)   # meridional angle
    if theta < theta_min:
        exclude_elems.add(elemLabel)
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

for k in sorted(layer_to_elems.keys()):
    set_name = 'PLY_%02d' % k
    p.SetFromElementLabels(
        name=set_name,
        elementLabels=layer_to_elems[k]
    )

session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    datumCoordSystems=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(
datumCoordSystems=OFF)


def build_node_label_map(part):
    return {nd.label: nd for nd in part.nodes}

def create_interface_node_set_by_delta_band(part, ply_k, N_part,
                                            r_inner, z_inner, thick, n,
                                            tol_delta=None, expand_factor=0.35):
    """
    Creates interface node set for boundary between ply k and k+1 by selecting
    ALL nodes whose delta is within a tolerance band around delta_int.
    Much more robust than per-element node picking.
    expand_factor: fraction of ply thickness used as default tolerance if tol_delta=None
    """
    if ply_k >= N_part - 1:
        print("PLY_%02d is last ply -> no interface set created." % ply_k)
        return None
    ply_set_name  = "PLY_%02d" % ply_k
    node_set_name = "PLY_%02d_IFACE" % ply_k
    if ply_set_name not in part.sets:
        raise RuntimeError("Set '%s' not found in part.sets" % ply_set_name)
    elems = part.sets[ply_set_name].elements
    if len(elems) == 0:
        print("Set '%s' is empty." % ply_set_name)
        return None
    # Interface delta
    delta_int = (ply_k + 1.0) * thick / float(N_part)
    # Tolerance
    tply = thick / float(N_part)
    if tol_delta is None:
        tol_delta = expand_factor * tply
    label_to_node = build_node_label_map(part)
    # Collect all nodes connected to elements in this ply
    candidate_nodes = set()
    for e in elems:
        for nlab in e.connectivity:
            candidate_nodes.add(nlab)
    # Compute deltas (cache)
    delta_cache = {}
    iface_nodes = set()
    n_none = 0
    for nlab in candidate_nodes:
        nd = label_to_node[nlab]
        r = nd.coordinates[0]
        z = nd.coordinates[1]
        d = delta_cache.get(nlab, None)
        if d is None:
            d = solve_delta_for_point(r, z, r_inner, z_inner, thick, n)
            delta_cache[nlab] = d
        if d is None:
            n_none += 1
            continue
        if abs(d - delta_int) <= tol_delta:
            iface_nodes.add(nlab)
    if not iface_nodes:
        raise RuntimeError(
            "No interface nodes found for %s. Increase tol_delta or check delta solver."
            % ply_set_name
        )
    part.SetFromNodeLabels(
        name=node_set_name,
        nodeLabels=tuple(sorted(iface_nodes))
    )
    return node_set_name


for k in range(N_part - 1):
    create_interface_node_set_by_delta_band(
        part=p,
        ply_k=k,
        N_part=N_part,
        r_inner=r_inner, z_inner=z_inner,
        thick=thick, n=n,
        tol_delta=0.1*mesh_size,
        expand_factor=0.1
    )

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
    stepName = list(odb.steps.keys())[-1]
    step = odb.steps[stepName]
    frameIndex = -1
    frame = step.frames[frameIndex]
    if 'S' not in frame.fieldOutputs:
        raise RuntimeError("Stress output 'S' not found in ODB frame. Request stresses (S).")
    if 'E' not in frame.fieldOutputs:
        raise RuntimeError("Strain output 'E' not found in ODB frame. Request strains (E).")
    Sfield = frame.fieldOutputs['S'].getSubset(position=INTEGRATION_POINT)
    Efield = frame.fieldOutputs['E'].getSubset(position=INTEGRATION_POINT)
    # Exclude elements based on S12
    # S12_LIMIT = 0.02
    # exclude_elems_s12 = set()
    for sv in Sfield.values:
        ip = getattr(sv, 'integrationPoint', None)
        if ip is None:
            continue
        sd = sv.data
        s12 = float(sd[3]) if len(sd) > 3 else 0.0
        # if abs(s12) > S12_LIMIT:
        #     exclude_elems_s12.add(sv.elementLabel)
    E_by_key = {}
    for ev in Efield.values:
        if ev.elementLabel in exclude_elems:
            continue
        # if ev.elementLabel in exclude_elems_s12:
        #     continue
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
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([
            'odb','step','frameIndex','frameValue',
            'instance','elementLabel','integrationPoint','sectionPoint',
            # global stresses/strains
            'S11','S22','S33','S12','S13','S23',
            'E11','E22','E33','E12','E13','E23',
            # ply
            'phi_deg',
            # stresses/strains in Puck lamina axes (1=fiber,2=trans,3=thk)
            'S11_p','S22_p','S33_p','S12_p','S13_p','S23_p',
            'E11_p','E22_p','E33_p','E12_p','E13_p','E23_p',
            # criteria
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
            # if el in exclude_elems_s12:
            #     continue
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
            # ---- full 3D transform to Puck lamina axes ----
            s11p, s22p, s33p, s12p, s13p, s23p = global_to_puck_voigt(s11,s22,s33,s12,s13,s23, phi)
            key = (inst_name, el, ip, sp)
            if key in E_by_key:
                e11, e22, e33, e12, e13, e23 = E_by_key[key]
                e11p, e22p, e33p, e12p, e13p, e23p = global_to_puck_voigt(e11,e22,e33,e12,e13,e23, phi)
            else:
                missing_strain += 1
                e11=e22=e33=e12=e13=e23=float('nan')
                e11p=e22p=e33p=e12p=e13p=e23p=float('nan')
            tw = tsai_wu_index(s11p, s22p, s12p)
            fail_tw = 1 if tw >= 1.0 else 0
            if tw > max_tw:
                max_tw = tw; max_tw_elem = el; max_tw_ip = ip; max_tw_inst = inst_name
            # ---- ply strain limit (tension in fiber/trans) ----
            if not (math.isnan(e11p) or math.isnan(e22p)):
                e11t = e11p if e11p > 0.0 else 0.0
                e22t = e22p if e22p > 0.0 else 0.0
                max_tension = max(e11t, e22t)
                fail_eps = 1 if (e11p > strain_limit or e22p > strain_limit) else 0
                if max_tension > max_eps_tension:
                    max_eps_tension = max_tension; max_eps_elem = el; max_eps_ip = ip; max_eps_inst = inst_name
            else:
                max_tension = float('nan'); fail_eps = 0
            ff, ff_t, ff_c, s_f_eff = puck_ff(s11p, s22p, s33p)
            iff, iff_branch, theta_star = puck_iff(s22p, s33p, s12p, s13p, s23p)
            fail_ff = 1 if ff >= 1.0 else 0
            fail_iff = 1 if iff >= 1.0 else 0
            fail_any = 1 if (fail_tw or fail_eps or fail_ff or fail_iff) else 0
            # Track worst across all criteria
            worst = max(tw,
                        max_tension if not math.isnan(max_tension) else -1.0,
                        ff, iff)
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
                max_puck = worst_puck; max_puck_elem = el; max_puck_ip = ip; max_puck_inst = inst_name
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


################################## Seperate - Tsai Wu, Strain and Pucks code ####################################
def tw_and_tension_strain_limit_integration_point(strain_limit=0.005):
    odb = openOdb(odbPath, readOnly=True)
    stepName = list(odb.steps.keys())[-1]
    step = odb.steps[stepName]
    frameIndex = -1
    frame = step.frames[frameIndex]
    if 'S' not in frame.fieldOutputs:
        raise RuntimeError("Stress output 'S' not found in ODB frame. Request stresses (S).")
    if 'E' not in frame.fieldOutputs:
        raise RuntimeError("Strain output 'E' not found in ODB frame. Request strains (E).")
    Sfield = frame.fieldOutputs['S'].getSubset(position=INTEGRATION_POINT)
    Efield = frame.fieldOutputs['E'].getSubset(position=INTEGRATION_POINT)
    E_by_key = {}
    for ev in Efield.values:
        if 'exclude_elems' in globals() and ev.elementLabel in exclude_elems:
            continue
        inst_name = ev.instance.name if ev.instance else ''
        ip = getattr(ev, 'integrationPoint', None)
        sp = getattr(ev, 'sectionPoint', None)
        if ip is None:
            continue
        e11 = float(ev.data[0])  # rr
        e22 = float(ev.data[1])  # zz
        e33 = float(ev.data[2])  # tt
        e12 = float(ev.data[3]) if len(ev.data) > 3 else 0.0  # rz
        E_by_key[(inst_name, ev.elementLabel, ip, sp)] = (e11, e22, e33, e12)
    base = os.path.splitext(os.path.basename(odbPath))[0]
    out_csv = 'tw_tensionStrainLimit_%s_axisym_IP.csv' % base
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([
            'odb','step','frameIndex','frameValue',
            'instance','elementLabel','integrationPoint','sectionPoint',
            'S_rr(S11)','S_zz(S22)','S_tt(S33)','T_rz(S12)',
            'E_rr(E11)','E_zz(E22)','E_tt(E33)','E_rz(E12)',
            'TW_zz_tt','fails_TW(>=1)',
            'max_tension(E22,E33)','strain_limit','fails_tension_strain(>limit)',
            'fails_any'
        ])
        max_tw = -1.0
        max_tw_elem = None
        max_tw_ip = None
        max_tw_inst = ''
        max_eps_tension = -1.0
        max_eps_elem = None
        max_eps_ip = None
        max_eps_inst = ''
        missing_strain = 0
        for sv in Sfield.values:
            if 'exclude_elems' in globals() and sv.elementLabel in exclude_elems:
                continue
            inst_name = sv.instance.name if sv.instance else ''
            ip = getattr(sv, 'integrationPoint', None)
            if ip is None:
                continue
            s11 = float(sv.data[0])  # rr
            s22 = float(sv.data[1])  # zz
            s33 = float(sv.data[2])  # tt
            s12 = float(sv.data[3]) if len(sv.data) > 3 else 0.0  # rz
            s23 = 0.0
            tw = tsai_wu_index(s22, s33, s23)
            fail_tw = 1 if tw >= 1.0 else 0
            if tw > max_tw:
                max_tw = tw
                max_tw_elem = sv.elementLabel
                max_tw_ip = ip
                max_tw_inst = inst_name
            key = (inst_name, sv.elementLabel, ip, sp)
            if key in E_by_key:
                e11, e22, e33, e12 = E_by_key[key]
            else:
                missing_strain += 1
                e11 = e22 = e33 = e12 = float('nan')
            if not (math.isnan(e22) or math.isnan(e33)):
                e22_t = e22 if e22 > 0.0 else 0.0
                e33_t = e33 if e33 > 0.0 else 0.0
                max_tension = max(e22_t, e33_t)
                fail_eps = 1 if (e22 > strain_limit or e33 > strain_limit) else 0
                if max_tension > max_eps_tension:
                    max_eps_tension = max_tension
                    max_eps_elem = sv.elementLabel
                    max_eps_ip = ip
                    max_eps_inst = inst_name
            else:
                max_tension = float('nan')
                fail_eps = 0
            fail_any = 1 if (fail_tw or fail_eps) else 0
            w.writerow([
                base, step.name, frameIndex, frame.frameValue,
                inst_name, sv.elementLabel, ip, sp,
                s11, s22, s33, s12,
                e11, e22, e33, e12,
                tw, fail_tw,
                max_tension, strain_limit, fail_eps,
                fail_any
            ])
    odb.close()
    print("Wrote:", out_csv)
    if missing_strain:
        print("WARNING: %d integration-point stress entries had no matching integration-point strain entry." % missing_strain)
    print("Max Tsai–Wu (IP): %.4f" % max_tw)
    print("  Instance      :", max_tw_inst)
    print("  Element label :", max_tw_elem)
    print("  Integration pt:", max_tw_ip)
    print("Max tension max(E22+,E33+) (IP): %.6f (limit %.6f)" % (max_eps_tension, strain_limit))
    print("  Instance      :", max_eps_inst)
    print("  Element label :", max_eps_elem)
    print("  Integration pt:", max_eps_ip)

# tw_and_tension_strain_limit_integration_point(strain_limit=0.005)


def tsai_wu_failure_lamina(elem_to_phi):
    odb = openOdb(odbPath, readOnly=True)
    stepName = list(odb.steps.keys())[-1]
    step = odb.steps[stepName]
    frameIndex = -1
    frame = step.frames[frameIndex]
    if 'S' not in frame.fieldOutputs:
        raise RuntimeError("Stress output 'S' not found in ODB frame.")
    Sfield = frame.fieldOutputs['S'].getSubset(position=INTEGRATION_POINT)
    base = os.path.splitext(os.path.basename(odbPath))[0]
    out_csv = 'tsaiwu_%s_lamina_elementNodal.csv' % base
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([
            'odb','step','frameIndex','frameValue',
            'instance','elementLabel',
            'phi_deg',
            'S22_mer','S33_hoop','S23_mer_hoop',
            'sig_fiber','sig_trans','tau_f_t',
            'TsaiWu','fails(>=1)'
        ])
        max_tw = -1.0
        max_elem = None
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
            phi = float(elem_to_phi[el])
            s11 = float(sd[0]) if len(sd) > 0 else 0.0
            s22 = float(sd[1]) if len(sd) > 1 else 0.0
            s33 = float(sd[2]) if len(sd) > 2 else 0.0
            s12 = float(sd[3]) if len(sd) > 3 else 0.0
            s13 = float(sd[4]) if len(sd) > 4 else 0.0
            s23 = float(sd[5]) if len(sd) > 5 else 0.0
            phi = float(elem_to_phi[el])
            # ---- full 3D transform to Puck lamina axes ----
            s11p, s22p, s33p, s12p, s13p, s23p = global_to_puck_voigt(s11,s22,s33,s12,s13,s23, phi)
            tw = tsai_wu_index(s11p, s22p, s12p)
            fail = 1 if tw >= 1.0 else 0
            if tw > max_tw:
                max_tw = tw
                max_elem = el
            w.writerow([
                base, step.name, frameIndex, frame.frameValue,
                inst_name, el,
                phi,
                s22, s33, s23,
                s11p, s22p, s12p,
                tw, fail
            ])
    odb.close()
    print("Wrote:", out_csv)
    print("Max Tsai–Wu (lamina): %.4f" % max_tw)
    print("  Element label :", max_elem)

# tsai_wu_failure_lamina(elem_to_phi)


def plane_strain_limit_check(limit=0.005):
    """
    Checks ELEMENT_NODAL strains and flags failure if either E22 or E33 exceeds limit
    at any element nodal value.
    """
    odb = openOdb(odbPath, readOnly=True)
    stepName = list(odb.steps.keys())[-1]
    step = odb.steps[stepName]
    frameIndex = -1
    frame = step.frames[frameIndex]
    if 'E' not in frame.fieldOutputs:
        raise RuntimeError("Strain output 'E' not found in ODB frame. Request strains (E) in your output.")
    Efield = frame.fieldOutputs['E'].getSubset(position=ELEMENT_NODAL)
    base = os.path.splitext(os.path.basename(odbPath))[0]
    out_csv = 'planestrainLimit_%s_axisym_elementNodal.csv' % base
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([
            'odb','step','frameIndex','frameValue',
            'instance','elementLabel','nodeLabel',
            'E_rr(E11)','E_zz(E22)','E_tt(E33)','E_rz(E12)',
            'max(E22,E33)','limit','fail(>limit)'
        ])
        global_max = -1.0
        max_elem = None
        max_node = None
        max_inst = ''
        for v in Efield.values:
            if 'exclude_elems' in globals() and v.elementLabel in exclude_elems:
                continue
            e11 = float(v.data[0])  # rr
            e22 = float(v.data[1])  # zz  (your E22)
            e33 = float(v.data[2])  # tt  (your E33)
            e12 = float(v.data[3]) if len(v.data) > 3 else 0.0  # rz
            max_inplane = max(e22, e33)
            fail = 1 if (e22 > limit or e33 > limit) else 0
            if max_inplane > global_max:
                global_max = max_inplane
                max_elem = v.elementLabel
                max_node = v.nodeLabel
                max_inst = v.instance.name if v.instance else ''
            inst_name = v.instance.name if v.instance else ''
            w.writerow([
                base, step.name, frameIndex, frame.frameValue,
                inst_name, v.elementLabel, v.nodeLabel,
                e11, e22, e33, e12,
                max_inplane, limit, fail
            ])
    odb.close()
    print("Wrote:", out_csv)
    print("Global max max(E22,E33): %.6f (limit %.6f)" % (global_max, limit))
    print("  Instance      :", max_inst)
    print("  Element label :", max_elem)
    print("  Node label    :", max_node)

# plane_strain_limit_check(limit=0.005)
# tsai_wu_failure()

def puck_failure_integration_points():
    odb = openOdb(odbPath, readOnly=True)
    stepName = list(odb.steps.keys())[-1]
    step = odb.steps[stepName]
    frameIndex = -1
    frame = step.frames[frameIndex]
    if 'S' not in frame.fieldOutputs:
        raise RuntimeError("Stress output 'S' not found in ODB frame. Request stresses (S).")
    Sfield = frame.fieldOutputs['S'].getSubset(position=INTEGRATION_POINT)
    base = os.path.splitext(os.path.basename(odbPath))[0]
    out_csv = 'puck_%s_integrationPoint.csv' % base
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([
            'odb','step','frameIndex','frameValue',
            'instance','elementLabel',
            'integrationPoint',
            'S11','S22','S33','S12','S13','S23',
            'FF','FF_t','FF_c','S11_eff',
            'IFF','IFF_branch','theta_crit_deg',
            'FAIL_FF(>=1)','FAIL_IFF(>=1)','FAIL_any'
        ])
        max_any = -1.0
        max_elem = None
        max_ip = None
        for v in Sfield.values:
            if v.elementLabel in exclude_elems:
                continue
            data = v.data
            s11 = float(data[0])
            s22 = float(data[1])
            s33 = float(data[2])
            s12 = float(data[3]) if len(data) > 3 else 0.0
            s13 = float(data[4]) if len(data) > 4 else 0.0
            s23 = float(data[5]) if len(data) > 5 else 0.0
            ff, ff_t, ff_c, s11_eff = puck_ff(s22, s11, s33)
            iff, iff_branch, theta_star = puck_iff(s33, s11, s23, s12, s13)
            fail_ff = 1 if ff >= 1.0 else 0
            fail_iff = 1 if iff >= 1.0 else 0
            fail_any = 1 if (fail_ff or fail_iff) else 0
            worst = max(ff, iff)
            if worst > max_any:
                max_any = worst
                max_elem = v.elementLabel
                max_ip = getattr(v, 'integrationPoint', None)
            inst_name = v.instance.name if v.instance else ''
            ip = getattr(v, 'integrationPoint', None)
            w.writerow([
                base, step.name, frameIndex, frame.frameValue,
                inst_name, v.elementLabel,
                ip,
                s11, s22, s33, s12, s13, s23,
                ff, ff_t, ff_c, s11_eff,
                iff, iff_branch, theta_star,
                fail_ff, fail_iff, fail_any
            ])
    odb.close()
    print("Wrote:", out_csv)
    print("Max(Puck FF/IFF): %.4f" % max_any)
    print("  Element label     :", max_elem)
    print("  Integration point :", max_ip)

puck_failure_integration_points()


# Abaqus Python 2.7 script
# Extract COORD + Stress (S) for a node set and write CSV
#
# Usage (inside Abaqus/CAE or abaqus python):
#   - Set odb_path, instance_name, nodeset_name, step_name, frame_index, out_csv
#
# Notes:
# - Stress is typically stored at INTEGRATION_POINT. We request NODAL subset which
#   averages to nodes when available.
# - COORD comes from the ODB instance node coordinates.

from odbAccess import openOdb
from abaqusConstants import NODAL
import csv

# ------------------ USER INPUTS ------------------
odb_path      = r"job-1.odb"
step_name     = "LoadingStep"
frame_index   = -1 
instance_name = "SUPERELLIPSOID_2D-1" 
nodeset_name  = "PLY_01_IFACE" 
out_csv       = r"nodeset_stress_coord.csv"

# Stress components to export:
# For 3D: 'S11','S22','S33','S12','S13','S23'
# For Axisym: typically 'S11','S22','S33','S12' (depends on output)
stress_labels = ["S11","S22","S33","S12","S13","S23"]
# --------------------------------------------------

def safe_get_component(label_list, tensor_data, label):
    """Return tensor component by label if present; else blank."""
    # Abaqus returns stress as a 6-tuple for 3D (S11,S22,S33,S12,S13,S23)
    # In some cases it may be 4-tuple (axisym/planar): (S11,S22,S33,S12) or similar
    idx_map = {
        "S11": 0, "S22": 1, "S33": 2,
        "S12": 3, "S13": 4, "S23": 5
    }
    i = idx_map.get(label, None)
    if i is None:
        return ""
    if i < len(tensor_data):
        return tensor_data[i]
    return ""

def main():
    odb = openOdb(path=odb_path, readOnly=True)
    # Get step/frame
    step = odb.steps[step_name]
    frame = step.frames[frame_index]
    # Instance + nodeset
    inst = odb.rootAssembly.instances[instance_name]
    # Node set can exist under rootAssembly.nodeSets (global) or instance.nodeSets (local)
    if nodeset_name in inst.nodeSets:
        nset = inst.nodeSets[nodeset_name]
    elif nodeset_name in odb.rootAssembly.nodeSets:
        nset = odb.rootAssembly.nodeSets[nodeset_name]
    else:
        odb.close()
        raise ValueError("Node set not found: %s (checked instance and rootAssembly)" % nodeset_name)
    # COORD: from ODB nodes
    # Build node label -> (x,y,z)
    coord_map = {}
    for n in nset.nodes:
        # n.coordinates is a tuple (x,y,z) or (x,y) depending on model
        coord_map[n.label] = n.coordinates
    # Stress field
    if "S" not in frame.fieldOutputs:
        odb.close()
        raise ValueError("Field output 'S' not found in frame. Ensure stress output is requested.")
    S = frame.fieldOutputs["S"]
    # Try to get nodal stress directly/averaged-to-nodes
    try:
        S_nodal = S.getSubset(region=nset, position=ELEMENT_NODAL)
        s_vals = S_nodal.values
    except:
        # Fallback: just subset by region (may still be IP values)
        # If this happens, user likely needs to request NODAL output or accept IP export.
        S_sub = S.getSubset(region=nset)
        s_vals = S_sub.values
    # Build node label -> averaged stress tuple
    # If multiple contributions per node exist (from different elems), average them
    sum_map = {}
    cnt_map = {}
    for v in s_vals:
        # v.nodeLabel exists for nodal-position values (and for some averaged outputs)
        # If not present, this is likely integration-point-only data without node association.
        if not hasattr(v, "nodeLabel"):
            continue
        nl = v.nodeLabel
        data = v.data  # tuple
        if nl not in sum_map:
            sum_map[nl] = [0.0]*len(data)
            cnt_map[nl] = 0
        for i in range(len(data)):
            sum_map[nl][i] += data[i]
        cnt_map[nl] += 1
    avg_map = {}
    for nl in sum_map:
        c = float(cnt_map[nl])
        avg_map[nl] = tuple([x/c for x in sum_map[nl]])
    # Write CSV
    f = open(out_csv, "wb")  # Python 2.7 csv
    w = csv.writer(f)
    # Header
    # COORD size may be 2 or 3; write generic X,Y,(Z)
    # Determine max coordinate dimension present
    max_dim = 0
    for k in coord_map:
        if len(coord_map[k]) > max_dim:
            max_dim = len(coord_map[k])
    coord_cols = ["X","Y"] + (["Z"] if max_dim >= 3 else [])
    header = ["nodeLabel"] + coord_cols + stress_labels
    w.writerow(header)
    # Rows: iterate nodes in set order
    for n in nset.nodes:
        nl = n.label
        xyz = coord_map.get(nl, ())
        x = xyz[0] if len(xyz) > 0 else ""
        y = xyz[1] if len(xyz) > 1 else ""
        z = xyz[2] if len(xyz) > 2 else ""
        s = avg_map.get(nl, None)
        row = [nl, x, y]
        if max_dim >= 3:
            row.append(z)
        if s is None:
            # no stress found for this node (can happen at boundaries, missing output, etc.)
            row += [""] * len(stress_labels)
        else:
            for lab in stress_labels:
                row.append(safe_get_component(stress_labels, s, lab))
        w.writerow(row)
    f.close()
    odb.close()
    print("Wrote:", out_csv)

if __name__ == "__main__":
    main()