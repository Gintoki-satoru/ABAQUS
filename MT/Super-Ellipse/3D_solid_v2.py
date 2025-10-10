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
num_points = 20                 # resolution along curve
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

# ------------------------
# Case 1: phi = 0, theta = 0 → 90 deg
# ------------------------
theta_vals = [math.radians(i*90/num_points) for i in range(num_points+1)]
phi_case1 = 0.0

inner_points_case1 = []
outer_points_case1 = []
for theta in theta_vals:
    inner_points_case1.append(superellipsoid_point_3d(phi_case1, theta, a, b, c, n1, n2))
    outer_points_case1.append(superellipsoid_point_3d(phi_case1, theta, a_out, b_out, c_out, n1, n2))

# ------------------------
# Case 2: phi = 0 → 90 deg, theta = 0
# ------------------------
phi_vals = [math.radians(i*90/num_points) for i in range(num_points+1)]
theta_case2 = 0.0

inner_points_case2 = []
outer_points_case2 = []
for phi in phi_vals:
    inner_points_case2.append(superellipsoid_point_3d(phi, theta_case2, a, b, c, n1, n2))
    outer_points_case2.append(superellipsoid_point_3d(phi, theta_case2, a_out, b_out, c_out, n1, n2))

# ------------------------
# Case 3: phi = 0 → 90 deg, theta = 90 deg
# ------------------------
theta_case3 = math.radians(90)

inner_points_case3 = []
outer_points_case3 = []
for phi in phi_vals:
    inner_points_case3.append(superellipsoid_point_3d(phi, theta_case3, a, b, c, n1, n2))
    outer_points_case3.append(superellipsoid_point_3d(phi, theta_case3, a_out, b_out, c_out, n1, n2))

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
mdb.models['SuperEllipse'].parts['SuperEllipsoid'].setValues(
    geometryRefinement=FINE)

############ Partitioning ############
'''
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
csys = p.DatumCsysByThreePoints(
    name='Datum csys-1',
    coordSysType=CARTESIAN,
    origin=(0.0, 0.0, 0.0),
    point1=(1.0, 0.0, 1.0),
    point2=(0.0, 1.0, 0.0)
)
csys_id = csys.id

dp = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
datumPlane0 = p.datums[dp.id]
angle_increment = 90.0 / n_long
d = p.datums
for i in range(1, n_long):
    angle = i * angle_increment
    dp_rot = p.DatumPlaneByRotation(plane=datumPlane0, axis=d[csys_id].axis3, angle=angle)
    dp_rot_obj = p.datums[dp_rot.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_rot_obj, cells=current_cells)

dp = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
datumPlane0 = p.datums[dp.id]
angle_increment = 90.0 / n_long
d = p.datums
for i in range(1, n_long):
    angle = i * angle_increment
    dp_rot = p.DatumPlaneByRotation(plane=datumPlane0, axis=d[csys_id].axis2, angle=angle)
    dp_rot_obj = p.datums[dp_rot.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_rot_obj, cells=current_cells)
'''

#### Strat - 2 ####

def surface_normal(phi, theta, a, b, c, n1, n2, h=1e-6):
    """Approximate normal via finite differences."""
    p0 = superellipsoid_point_3d(phi, theta, a, b, c, n1, n2)
    p_phi = superellipsoid_point_3d(phi + h, theta, a, b, c, n1, n2)
    p_theta = superellipsoid_point_3d(phi, theta + h, a, b, c, n1, n2)
    t_phi = [p_phi[i] - p0[i] for i in range(3)]
    t_theta = [p_theta[i] - p0[i] for i in range(3)]
    nx = t_phi[1]*t_theta[2] - t_phi[2]*t_theta[1]
    ny = t_phi[2]*t_theta[0] - t_phi[0]*t_theta[2]
    nz = t_phi[0]*t_theta[1] - t_phi[1]*t_theta[0]
    norm = math.sqrt(nx**2 + ny**2 + nz**2)
    return (nx/norm, ny/norm, nz/norm)

def offset_point_along_normal(point, normal, t):
    return tuple(point[i] + t*normal[i] for i in range(3))

'''# Get the part and faces
f = p.faces
num_partitions = 4
theta_case = math.radians(90)
phi_vals = [i * math.pi/2 / num_partitions for i in range(1, num_partitions)]

for phi_example in phi_vals:
    pt_inner = superellipsoid_point_3d(phi_example, theta_case, a, b, c, n1, n2)
    pt_outer = superellipsoid_point_3d(phi_example, theta_case, a_out, b_out, c_out, n1, n2)
    n_vec = surface_normal(phi_example, theta_case, a, b, c, n1, n2)
    # Slight adjustment along normal (reduces Abaqus tolerance errors)
    pt_inner_adj = offset_point_along_normal(pt_inner, n_vec, -1e-4)
    pt_outer_adj = offset_point_along_normal(pt_outer, n_vec, +1e-4)
    datum_inner = p.DatumPointByCoordinate(coords=pt_inner_adj)
    datum_outer = p.DatumPointByCoordinate(coords=pt_outer_adj)
    # Attempt with increasing search tolerance
    closest_face_result = f.getClosest(coordinates=(pt_inner_adj,), searchTolerance=1e-3)
    if not closest_face_result:
        closest_face_result = f.getClosest(coordinates=(pt_inner_adj,), searchTolerance=1e-2)
    if closest_face_result:
        closest_face = closest_face_result[0][0]
        p.PartitionFaceByShortestPath(
            faces=closest_face,
            point1=p.datums[datum_inner.id],
            point2=p.datums[datum_outer.id]
        )
    else:
        print("No nearby face found at:", pt_inner_adj)

phi_case = 0.0
num_partitions = 4
theta_vals = [i * math.pi/2 / num_partitions for i in range(1, num_partitions + 1)]
for theta_example in theta_vals:   
    # 1️⃣ Compute points on the surface
    pt_inner = superellipsoid_point_3d(phi_case, theta_example, a, b, c, n1, n2)
    pt_outer = superellipsoid_point_3d(phi_case, theta_example, a_out, b_out, c_out, n1, n2)
    # 2️⃣ Compute surface normal
    n_vec = surface_normal(phi_case, theta_example, a, b, c, n1, n2)
    # 3️⃣ Adjust for numerical tolerance
    pt_inner_adj = offset_point_along_normal(pt_inner, n_vec, -1e-4)
    pt_outer_adj = offset_point_along_normal(pt_outer, n_vec, +1e-4)
    # 4️⃣ Create datum points
    datum_inner = p.DatumPointByCoordinate(coords=pt_inner_adj)
    datum_outer = p.DatumPointByCoordinate(coords=pt_outer_adj)
    # 5️⃣ Find closest face (with tolerance)
    closest_face_result = f.getClosest(coordinates=(pt_inner_adj,), searchTolerance=1e-3)
    if not closest_face_result:
        closest_face_result = f.getClosest(coordinates=(pt_inner_adj,), searchTolerance=1e-2)
    # 6️⃣ Partition along the line between inner & outer
    if closest_face_result:
        closest_face = closest_face_result[0][0]
        p.PartitionFaceByShortestPath(
            faces=closest_face,
            point1=p.datums[datum_inner.id],
            point2=p.datums[datum_outer.id]
        )
    else:
        print("No nearby face found at:", pt_inner_adj)
'''
####### Strat - 3 ######
# Get the part and faces
num_partitions = 4
theta_case = math.radians(90)
phi_vals = [i * math.pi/2 / num_partitions for i in range(1, num_partitions)]

for phi_example in phi_vals:
    # 1️⃣ Compute inner and outer points
    pt_inner = superellipsoid_point_3d(phi_example, theta_case, a, b, c, n1, n2)
    pt_outer = superellipsoid_point_3d(phi_example, theta_case, a_out, b_out, c_out, n1, n2)
    # 2️⃣ Compute normal and adjust points slightly
    n_vec = surface_normal(phi_example, theta_case, a, b, c, n1, n2)
    pt_inner_adj = offset_point_along_normal(pt_inner, n_vec, -1e-4)
    pt_outer_adj = offset_point_along_normal(pt_outer, n_vec, +1e-4)
    # 3️⃣ Create datum points
    datum_inner = p.DatumPointByCoordinate(coords=pt_inner_adj)
    datum_outer = p.DatumPointByCoordinate(coords=pt_outer_adj)
    datum_xy = p.DatumPointByCoordinate(coords=(pt_inner_adj[2], pt_inner_adj[1], pt_inner_adj[0]))
    # 4️⃣ Create datum plane by three points
    plane_datum = p.DatumPlaneByThreePoints(
        point1=p.datums[datum_inner.id],
        point2=p.datums[datum_outer.id],
        point3=p.datums[datum_xy.id]
    )
    # 5️⃣ Select entire part using bounding box
    cell_to_partition = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    # 6️⃣ Partition the selected cell using the datum plane
    p.PartitionCellByDatumPlane(
        cells=cell_to_partition,
        datumPlane=p.datums[plane_datum.id]
    )

