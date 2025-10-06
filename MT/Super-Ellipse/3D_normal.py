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
a, b, c = 150.0, 100.0, 150.0   # semi-axes
t = 2.5                         # thickness
n1, n2 = 2, 2           # shape exponents
num_points = 30                 # resolution along curve
n_long = 4  # number of longitudinal partitions

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
# Case 1: phi = 0, theta = 0 → 90 deg
# ------------------------
theta_vals = [math.radians(i*90/num_points) for i in range(num_points+1)]
phi_case1 = 0.0

inner_points_case1 = []
outer_points_case1 = []
for theta in theta_vals:
    p_inner = superellipsoid_point_3d(phi_case1, theta, a, b, c, n1, n2)
    n_vec = surface_normal(phi_case1, theta, a, b, c, n1, n2)
    p_outer = offset_point_along_normal(p_inner, n_vec, t)
    inner_points_case1.append(p_inner)
    outer_points_case1.append(p_outer)

# ------------------------
# Case 2: phi = 0 → 90 deg, theta = 0
# ------------------------
phi_vals = [math.radians(i*90/num_points) for i in range(num_points+1)]
theta_case2 = 0.0

inner_points_case2 = []
outer_points_case2 = []
for phi in phi_vals:
    p_inner = superellipsoid_point_3d(phi, theta_case2, a, b, c, n1, n2)
    n_vec = surface_normal(phi, theta_case2, a, b, c, n1, n2)
    p_outer = offset_point_along_normal(p_inner, n_vec, t)
    inner_points_case2.append(p_inner)
    outer_points_case2.append(p_outer)

# ------------------------
# Case 3: phi = 0 → 90 deg, theta = 90 deg
# ------------------------
theta_case3 = math.radians(90)

inner_points_case3 = []
outer_points_case3 = []
for phi in phi_vals:
    p_inner = superellipsoid_point_3d(phi, theta_case3, a, b, c, n1, n2)
    n_vec = surface_normal(phi, theta_case3, a, b, c, n1, n2)
    p_outer = offset_point_along_normal(p_inner, n_vec, t)
    inner_points_case3.append(p_inner)
    outer_points_case3.append(p_outer)

# Case2[-1] = Case3[-1]
inner_points_case2[-1] = inner_points_case3[-1]
outer_points_case2[-1] = outer_points_case3[-1]

# ------------------------
# Create datum points
# ------------------------
def create_datums(point_list):
    for pt in point_list:
        p.DatumPointByCoordinate(coords=pt)

# Case 1
create_datums(inner_points_case1 + outer_points_case1)

# Case 2
create_datums(inner_points_case2 + outer_points_case2)

# Case 3
create_datums(inner_points_case3 + outer_points_case3)

# ------------------------
# Create a wire spline using these datum points
# ------------------------

datum_inner_case1 = []
for pt in inner_points_case1:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in outer_points_case1:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in inner_points_case2:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in outer_points_case2:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in inner_points_case3:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in outer_points_case3:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

dp_inner_case1 = inner_points_case1
dp_inner_case2 = inner_points_case2
dp_inner_case3 = inner_points_case3

dp_outer_case1 = outer_points_case1
dp_outer_case2 = outer_points_case2
dp_outer_case3 = outer_points_case3

mid_idx = len(dp_inner_case1) // 2

e1 = p.edges
p.SolidLoft(loftsections=((e1.findAt(coordinates=(dp_outer_case1[mid_idx])), e1.findAt(coordinates=(dp_outer_case2[mid_idx])), e1.findAt(coordinates=( dp_outer_case3[mid_idx]))), 
                          (e1.findAt(coordinates=(dp_inner_case1[mid_idx])), e1.findAt(coordinates=(dp_inner_case2[mid_idx])), e1.findAt( coordinates=(dp_inner_case3[mid_idx])))), 
                          startCondition=NONE, endCondition=NONE)

p.regenerate()

############ Material properties ############
model.Material(name='Aluminium')
model.materials['Aluminium'].Elastic(table=((70000.0,0.34), ))

############ Section assignment ############
model.HomogeneousSolidSection(name='AL_section', 
    material='Aluminium', thickness=None)
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']

long = math.pi / 4      # longitude
lat = math.pi / 4  # latitude
q = a + t/2
w = b + t/2
e = c + t/2
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
"""p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
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
    dp_rot = p.DatumPlaneByRotation(plane=datumPlane0, axis=d[262].axis3, angle=angle)
    dp_rot_obj = p.datums[dp_rot.id]
    current_cells = p.faces.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionFaceByDatumPlane(datumPlane=dp_rot_obj, faces=current_cells)"""

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
    
    