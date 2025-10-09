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
c = p.cells
cells = c.findAt(((x,y,z), ))
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


p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=152.5)
p.DatumPointByCoordinate(coords=(0.0, 152.5, 0.0))
f, e, d = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=d[4], sketchUpEdge=e.findAt(coordinates=(
    140.891557, 0.0, 58.359393)), sketchPlaneSide=SIDE1, origin=(16.38815, 
    152.5, 16.38815))
s = mdb.models['SuperEllipse'].ConstrainedSketch(name='__profile__', 
    sheetSize=914.04, gridSpacing=22.85, transform=t)
g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.CircleByCenterPerimeter(center=(-16.6340125191169, 16.1385423244747), 
    point1=(5.7125, 0.0))
s.ObliqueDimension(vertex1=v.findAt((-16.634013, 16.138542)), vertex2=v.findAt(
    (5.7125, 0.0)), textPoint=(18.6269860775, 55.8358914094314), value=50.0)
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
f = p.faces
pickedFaces = f.findAt(((65.559857, 120.914177, 65.559857), ), ((0.0, 
    152.229696, 2.494274), ), ((152.229696, 2.494274, 0.0), ), ((65.291169, 
    120.418627, 65.291169), ), ((0.0, 151.60495, 2.484051), ), ((151.60495, 
    2.484051, 0.0), ), ((65.022481, 119.923077, 65.022481), ), ((0.0, 
    150.980204, 2.473829), ), ((150.980204, 2.473829, 0.0), ), ((64.753793, 
    119.427527, 64.753793), ), ((64.485105, 118.931977, 64.485105), ), ((0.0, 
    150.355459, 2.463606), ), ((150.355459, 2.463606, 0.0), ))
f1, e1, d2 = p.faces, p.edges, p.datums
p.PartitionFaceBySketchThruAll(sketchPlane=d2[4], sketchUpEdge=e1.findAt(
    coordinates=(140.891557, 0.0, 58.359393)), faces=pickedFaces, 
    sketchPlaneSide=SIDE1, sketch=s)
s.unsetPrimaryObject()
del mdb.models['SuperEllipse'].sketches['__profile__']


