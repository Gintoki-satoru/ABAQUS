# -*- coding: utf-8 -*-

from abaqus import *
from abaqusConstants import *
from odbAccess import *
import sketch
import os
import numpy as np

import section
import regionToolset
import displayGroupMdbToolset as dgm
import mesh
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

r_inner = 158.71     # semi-axis in a
z_inner = 476.14     # semi-axis in c
thick   = 2.009

r_outer = r_inner + thick
z_outer = z_inner + thick

n = 0.85

n_theta = 80  # number of spline points

N_part = 4  # Number of partitions

mesh_size = 0.5  # Mesh size

Press = 1 # Pressure load

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


theta_vals = np.linspace(0.0, math.pi/2.0, n_theta)

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

# Outer profile (reverse so sketch closes cleanly)
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
    # axis of symmetry (required)
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

# def surface_normal_rz(r, z, a, c, n, eps=1e-12):
#     p = 0.5 * float(n)   # p = n/2
#     ra = max(r / a, eps)
#     zc = max(z / c, eps)
#     nr = (p / a) * (ra ** (p - 1.0))
#     nz = (p / c) * (zc ** (p - 1.0))
#     mag = math.sqrt(nr*nr + nz*nz)
#     return (nr / mag, nz / mag)

def surface_normal_rz(r, z, a, c, n, eps=1e-12):
    # gradient components
    # clamp to avoid (0)^(negative) when n<1
    ra = max(r / a, eps)
    zc = max(z / c, eps)
    e = 2/n
    dr = (e / a) * (ra ** (e - 1.0))
    dz = (e / c) * (zc ** (e - 1.0))
    mag = math.sqrt(dr*dr + dz*dz)
    return (dr / mag, dz / mag)   # (n_r, n_z)

def surface_tangent_rz(r, z, a, c, e, eps=1e-12):
    n_r, n_z = surface_normal_rz(r, z, a, c, e, eps)
    # 90° rotation in r–z plane
    t_r = -n_z
    t_z =  n_r
    return (t_r, t_z)


N_theta = 6          # number of meridional regions you want
delta = 2.0 * thick # normal offset distance (safe choice)

theta_vals = np.linspace(0.0, math.pi/2.0, N_theta + 1)

# use interior midpoints only
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
    d3 = p.DatumPointByCoordinate(coords=p_mid)
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
    name='Section-1', 
    material='AL', 
    thickness=None
)

p = mdb.models[modelName].parts['SuperEllipsoid_2D']

faces = p.faces.getByBoundingBox(
    xMin=-1e6, yMin=-1e6, zMin=-1e6,
    xMax= 1e6, yMax= 1e6, zMax= 1e6
)

region = regionToolset.Region(faces=faces)

p.SectionAssignment(
    region=region,
    sectionName='Section-1',
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

################# SET creation ###################

def unique(seq):
    out, seen = [], set()
    for obj in seq:
        key = obj.index  # edge index is stable within the part
        if key not in seen:
            seen.add(key)
            out.append(obj)
    return out

e = p.edges
theta0_edges = []
for k in range(N_part):
    r_mid = r_inner + (k + 0.5) * thick / float(N_part)
    z_mid = 0.0  # tiny offset to avoid picking exactly on a vertex
    # In axisymmetric 2D geometry, use (r, z, 0.0)
    theta0_edges.append(e.findAt(((r_mid, z_mid, 0.0),)))

theta90_edges = []
for k in range(N_part):
    r_mid = 0.0  # tiny offset to avoid r=0 singular pick
    z_mid = z_inner + (k + 0.5) * thick / float(N_part)
    theta90_edges.append(e.findAt(((r_mid, z_mid, 0.0),)))

theta0_edges = unique(theta0_edges)
theta90_edges = unique(theta90_edges)
p.Set(name='set_bottom', edges=theta0_edges)
p.Set(name='set_top', edges=theta90_edges)

################# Surface creation ###################
N_sample = 40   # more than enough
theta_vals = np.linspace(0.0, math.pi/2.0, N_sample+2)[1:-1]

s2 = p.edges

search_tol = 1e-3          # start reasonably large
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
    # avoid boundaries
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

######## local coordinate system  ##########
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

# n_r, n_z = surface_normal_rz(r_c, z_c, r_inner, z_inner, n)

# # tangent = 90° rotation in r–z plane
# t_r = -n_z
# t_z =  n_r

# d1csys = p.DatumCsysByThreePoints(
#     name='CSYS_ELEM_%d' % elemLabel,
#     coordSysType=CARTESIAN,
#     origin=(r_c, z_c, 0.0),
#     point1=(r_c + t_r, z_c + t_z, 0.0),   # X = tangent
#     point2=(r_c + n_r, z_c + n_z, 0.0)    # Y = normal
# )

# elemSet = p.SetFromElementLabels(
#     name='ONE_ELEM_SET',
#     elementLabels=(elemLabel,)
# )

# p.MaterialOrientation(
#     region=elemSet,
#     orientationType=SYSTEM,
#     localCsys=p.datums[d1csys.id],
#     axis=AXIS_3,
#     angle=0.0,
#     additionalRotationType=ROTATION_ANGLE,
#     stackDirection=STACK_3
# )

for elem in p.elements:
    elemLabel = elem.label
    # --- element node coordinates ---
    node_ids = elem.connectivity
    coords = [p.nodes[i].coordinates for i in node_ids]
    r = [c[0] for c in coords]   # radial
    z = [c[1] for c in coords]   # axial
    # --- element centroid ---
    r_c = sum(r) / len(r)
    z_c = sum(z) / len(z)
    # --- surface normal at centroid ---
    n_r, n_z = surface_normal_rz(
        r_c, z_c,
        r_inner, z_inner,
        n
    )
    # --- tangent (90° rotation) ---
    t_r = -n_z
    t_z =  n_r
    # --- create local CSYS ---
    dcsys = p.DatumCsysByThreePoints(
        name='CSYS_ELEM_%d' % elemLabel,
        coordSysType=CARTESIAN,
        origin=(r_c, z_c, 0.0),
        point1=(r_c + t_r, z_c + t_z, 0.0),   # X = tangent
        point2=(r_c + n_r, z_c + n_z, 0.0)    # Y = normal
    )
    # --- element set ---
    elemSet = p.SetFromElementLabels(
        name='ELEM_SET_%d' % elemLabel,
        elementLabels=(elemLabel,)
    )
    # --- assign material orientation ---
    p.MaterialOrientation(
        region=elemSet,
        orientationType=SYSTEM,
        localCsys=p.datums[dcsys.id],
        axis=AXIS_3,
        angle=0.0,
        additionalRotationType=ROTATION_ANGLE,
        stackDirection=STACK_3
    )

# print("Assigned material orientation to", len(p.elements), "elements.")
