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


session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


################ Parameters #####################
a, b, c = 150.0, 100.0, 150.0   # semi-axes outer points
t = 2.5                         # thickness
n1, n2 = 2, 2                   # shape exponents
num_points = 20                 # resolution along curve
n_long = 4                      # number of longitudinal partitions
num_layers = 4                  # number of layers including inner & outer through thickness

a_out, b_out, c_out = a + t, b + t, c + t

Mdb()
modelName = 'SuperEllipse'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
model = mdb.models[modelName]
model.rootAssembly.clearGeometryCache()
model.rootAssembly.regenerate()
p = model.Part(name='SuperEllipsoid', dimensionality=THREE_D, type=DEFORMABLE_BODY)

################ Part Creation #####################

# ------------------------
# Geometry creation function
# ------------------------
def signed_power(base, exp):
    return math.copysign(abs(base)**exp, base)

def superellipsoid_point_3d(phi, theta, a, b, c, n1, n2):
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    x = a * signed_power(cos_phi, 2.0/n1) * signed_power(cos_theta, 2.0/n2)
    y = b * signed_power(cos_phi, 2.0/n1) * signed_power(sin_theta, 2.0/n2)
    z = c * signed_power(sin_phi, 2.0/n1)
    return x, y, z

def surface_normal(phi, theta, a, b, c, n1, n2, h=1e-6):
    """Approximate normal via finite differences."""
    # Base point
    p0 = superellipsoid_point_3d(phi, theta, a, b, c, n1, n2)
    # Slightly offset in phi and theta
    p_phi = superellipsoid_point_3d(phi + h, theta, a, b, c, n1, n2)
    p_theta = superellipsoid_point_3d(phi, theta + h, a, b, c, n1, n2)
    # Tangent vectors
    t_phi = [p_phi[i] - p0[i] for i in range(3)]
    t_theta = [p_theta[i] - p0[i] for i in range(3)]
    # Cross product t_phi × t_theta
    nx = t_phi[1]*t_theta[2] - t_phi[2]*t_theta[1]
    ny = t_phi[2]*t_theta[0] - t_phi[0]*t_theta[2]
    nz = t_phi[0]*t_theta[1] - t_phi[1]*t_theta[0]
    # Normalize
    norm = math.sqrt(nx**2 + ny**2 + nz**2)
    return (nx/norm, ny/norm, nz/norm)

def offset_point_along_normal(point, normal, t):
    return tuple(point[i] + t*normal[i] for i in range(3))

# ------------------------
# Layer offsets
# ------------------------
t_vals = [i * t/(num_layers-1) for i in range(num_layers)]  # 0 → t

# ------------------------
# Function to generate points for a case
# ------------------------
def generate_case_points(phi_vals, theta_vals):
    all_layers = [[] for _ in range(num_layers)]
    for phi in phi_vals:
        for j, theta in enumerate(theta_vals):
            p_base = superellipsoid_point_3d(phi, theta, a, b, c, n1, n2)
            n_vec = surface_normal(phi, theta, a, b, c, n1, n2)
            for i, t_i in enumerate(t_vals):
                all_layers[i].append(offset_point_along_normal(p_base, n_vec, t_i))
    return all_layers

# ------------------------
# Case 1: phi = 0, theta = 0 → 90 deg
# ------------------------
theta_vals_case1 = [math.radians(i*90/num_points) for i in range(num_points+1)]
phi_case1 = [0.0]
layers_case1 = generate_case_points(phi_case1, theta_vals_case1)

# ------------------------
# Case 2: phi = 0 → 90 deg, theta = 0
# ------------------------
phi_vals_case2 = [math.radians(i*90/num_points) for i in range(num_points+1)]
theta_case2 = [0.0]
layers_case2 = generate_case_points(phi_vals_case2, theta_case2)

# ------------------------
# Case 3: phi = 0 → 90 deg, theta = 90 deg
# ------------------------
theta_case3 = [math.radians(90)]
layers_case3 = generate_case_points(phi_vals_case2, theta_case3)

# ------------------------
# last point of case2 = last of case3
# ------------------------
for i in range(num_layers):
    layers_case2[i][-1] = layers_case3[i][-1]
    '''y_val = layers_case3[i][0][1]
    layers_case1[i][-1] = (0.0, y_val, 0.0)
    layers_case3[i][0] = (0.0, y_val, 0.0)
    layers_case1[i][0] = layers_case2[i][0]'''

# ------------------------
# Create datum points
# ------------------------
def create_datums(point_list):
    for pt in point_list:
        p.DatumPointByCoordinate(coords=pt)

# Add all layers for each case
for i in range(num_layers):
    create_datums(layers_case1[i])
    create_datums(layers_case2[i])
    create_datums(layers_case3[i])

all_wires = []

for layer_idx in range(num_layers):
    datum_layer = []
    # Case 1 points
    for pt in layers_case1[layer_idx]:
        dp = p.DatumPointByCoordinate(coords=pt)
        datum_layer.append(p.datums[dp.id])
    wire_case1 = p.WireSpline(points=datum_layer, mergeType=IMPRINT, meshable=ON, smoothClosedSpline=ON)
    all_wires.append(wire_case1)
    datum_layer = []
    # Case 2 points
    for pt in layers_case2[layer_idx]:
        dp = p.DatumPointByCoordinate(coords=pt)
        datum_layer.append(p.datums[dp.id])
    wire_case2 = p.WireSpline(points=datum_layer, mergeType=IMPRINT, meshable=ON, smoothClosedSpline=ON)
    all_wires.append(wire_case2)
    datum_layer = []
    # Case 3 points
    for pt in layers_case3[layer_idx]:
        dp = p.DatumPointByCoordinate(coords=pt)
        datum_layer.append(p.datums[dp.id])
    wire_case3 = p.WireSpline(points=datum_layer, mergeType=IMPRINT, meshable=ON, smoothClosedSpline=ON)
    all_wires.append(wire_case3)

e = p.edges
loft_sections = []
for layer_idx in range(num_layers):
    mid_idx = len(layers_case1[layer_idx]) // 2
    section = (
        e.findAt(coordinates=layers_case1[layer_idx][mid_idx]),
        e.findAt(coordinates=layers_case2[layer_idx][mid_idx]),
        e.findAt(coordinates=layers_case3[layer_idx][mid_idx])
    )
    loft_sections.append(section)

p.SolidLoft(loftsections=loft_sections, startCondition=NONE, endCondition=NONE)
p.regenerate()

############ Material properties ############
model.Material(name='Aluminium')
model.materials['Aluminium'].Elastic(table=((70000.0,0.34), ))

############ Section assignment ############
model.HomogeneousSolidSection(name='AL_section', 
    material='Aluminium', thickness=None)
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']

long = math.pi / 4      # longitude
lat = math.pi / 4       # latitude
q = a
w = b
e = c
x, y, z = superellipsoid_point_3d(lat, long, q, w, e, n1, n2)
picked_point = (x, y, z)
c1 = p.cells
cells = c1.findAt(((x,y,z), ))
region = regionToolset.Region(cells=cells)
p.SectionAssignment(region=region, sectionName='AL_section', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

############ Assembly ############
a1 = mdb.models['SuperEllipse'].rootAssembly
a1.DatumCsysByDefault(CARTESIAN)
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
a1.Instance(name='SuperEllipsoid-1', part=p, dependent=OFF)
session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)

############ Partitioning ############

#### Strat - 1 ####
'''p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
dp = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
datumPlane0 = p.datums[dp.id]

csys = p.DatumCsysByThreePoints(
    name='Datum csys-1',
    coordSysType=CARTESIAN,
    origin=(0.0, 0.0, 0.0),
    point1=(1.0, 0.0, 1.0),
    point2=(0.0, 1.0, 0.0)
)
csys_id = csys.id
d = p.datums

angle_increment = 90.0 / n_long
for i in range(1, n_long):
    angle = i * angle_increment
    dp_rot = p.DatumPlaneByRotation(
        plane=datumPlane0,
        axis=d[csys_id].axis2,
        angle=angle
    )
    dp_rot_obj = p.datums[dp_rot.id]
    current_cells = p.faces.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionFaceByDatumPlane(datumPlane=dp_rot_obj, faces=current_cells)

dp = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
datumPlane0 = p.datums[dp.id]

angle_increment = 90.0 / n_long
d = p.datums
for i in range(1, n_long):
    angle = i * angle_increment
    dp_rot = p.DatumPlaneByRotation(plane=datumPlane0, axis=d[csys_id].axis3, angle=angle)
    dp_rot_obj = p.datums[dp_rot.id]
    current_cells = p.faces.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionFaceByDatumPlane(datumPlane=dp_rot_obj, faces=current_cells)

#### Strat - 2 ####
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
dp_xy = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
datumPlane0 = p.datums[dp_xy.id]

n_partitions = 4
total_length = c
dz = total_length / n_partitions

for i in range(1, n_partitions):
    offset = i * dz
    dp_offset = p.DatumPlaneByOffset(plane=datumPlane0, flip=SIDE1, offset=offset)
    dp_offset_obj = p.datums[dp_offset.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_offset_obj, cells=current_cells)

dp_yz = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
datumPlane1 = p.datums[dp_yz.id]

for i in range(1, n_partitions):
    offset = i * dz
    dp_offset = p.DatumPlaneByOffset(plane=datumPlane1, flip=SIDE1, offset=offset)
    dp_offset_obj = p.datums[dp_offset.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_offset_obj, cells=current_cells)
    '''
    