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

Mdb()
modelName = 'SuperEllipse'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
model = mdb.models[modelName]
model.rootAssembly.clearGeometryCache()
model.rootAssembly.regenerate()

################ Parameters #####################
a, b, c = 150.0, 150.0, 150.0   # inner semi-axes
t = 2.5                         # total thickness
n1, n2 = 2, 2                   # shape exponents
num_points = 30                 # points per curve
num_layers = 4                  # number of layers through thickness

a_out, b_out, c_out = a + t, b + t, c + t  # outer semi-axes

################ Helper Functions #####################
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

def create_datums(p, point_list):
    for pt in point_list:
        p.DatumPointByCoordinate(coords=pt)

def create_wire_splines(part, datum_points, smoothClosed=True):
    wire = part.WireSpline(points=datum_points,
                           mergeType=IMPRINT,
                           meshable=ON,
                           smoothClosedSpline=ON if smoothClosed else OFF)
    return wire

def create_layer_part(model, layer_index):
    frac1 = float(layer_index) / num_layers
    frac2 = float(layer_index + 1) / num_layers
    # Semi-axes for inner and outer surfaces of this layer
    a_i, b_i, c_i = a + frac1*(a_out - a), b + frac1*(b_out - b), c + frac1*(c_out - c)
    a_o, b_o, c_o = a + frac2*(a_out - a), b + frac2*(b_out - b), c + frac2*(c_out - c)
    # Angular values
    phi_vals = [math.radians(i*90/num_points) for i in range(num_points+1)]
    theta_vals = [math.radians(i*90/num_points) for i in range(num_points+1)]
    # --- Case 1 (phi=0) ---
    inner_case1 = [superellipsoid_point_3d(0.0, th, a_i, b_i, c_i, n1, n2) for th in theta_vals]
    outer_case1 = [superellipsoid_point_3d(0.0, th, a_o, b_o, c_o, n1, n2) for th in theta_vals]
    # --- Case 2 (theta=0) ---
    inner_case2 = [superellipsoid_point_3d(ph, 0.0, a_i, b_i, c_i, n1, n2) for ph in phi_vals]
    outer_case2 = [superellipsoid_point_3d(ph, 0.0, a_o, b_o, c_o, n1, n2) for ph in phi_vals]
    # --- Case 3 (theta=90°) ---
    theta_90 = math.radians(90)
    inner_case3 = [superellipsoid_point_3d(ph, theta_90, a_i, b_i, c_i, n1, n2) for ph in phi_vals]
    outer_case3 = [superellipsoid_point_3d(ph, theta_90, a_o, b_o, c_o, n1, n2) for ph in phi_vals]
    # Fix overlapping endpoints
    inner_case2[-1] = inner_case3[-1]
    outer_case2[-1] = outer_case3[-1]
    part_name = "Layer_" + str(layer_index + 1)
    p = model.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    # Create datum points for each case
    create_datums(p, inner_case1 + outer_case1)
    create_datums(p, inner_case2 + outer_case2)
    create_datums(p, inner_case3 + outer_case3)
    # Create wire splines for each case
    create_wire_splines(p, inner_case1)
    create_wire_splines(p, outer_case1)
    create_wire_splines(p, inner_case2)
    create_wire_splines(p, outer_case2)
    create_wire_splines(p, inner_case3)
    create_wire_splines(p, outer_case3)
    mid_idx = len(inner_case1) // 2
    e1 = p.edges
    p.SolidLoft(loftsections=((e1.findAt(coordinates=(inner_case1[mid_idx])), e1.findAt(coordinates=(inner_case2[mid_idx])), e1.findAt(coordinates=( inner_case3[mid_idx]))), 
                          (e1.findAt(coordinates=(outer_case1[mid_idx])), e1.findAt(coordinates=(outer_case2[mid_idx])), e1.findAt( coordinates=(outer_case3[mid_idx])))), 
                          startCondition=NONE, endCondition=NONE)
    return p

model = mdb.models['SuperEllipse']
for i in range(num_layers):
    create_layer_part(model, i)

############ Create SuperEllipsoid by merging layers ############

def assemble_and_merge_layers(num_layers):
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    instance_names = []
    for i in range(1, num_layers + 1):
        part_name = "Layer_" + str(i)
        inst_name = "Layer_" + str(i) + "-1"
        part = model.parts[part_name]
        assembly.Instance(name=inst_name, part=part, dependent=ON)
        instance_names.append(inst_name)
    instances_to_merge = tuple(assembly.instances[name] for name in instance_names)
    assembly.InstanceFromBooleanMerge(
        name='SuperEllipsoid',
        instances=instances_to_merge,
        keepIntersections=ON,
        originalInstances=DELETE,
        domain=GEOMETRY
    )

assemble_and_merge_layers(num_layers)

############ Material properties ############
model.Material(name='Aluminium')
model.materials['Aluminium'].Elastic(table=((70000.0,0.34), ))

############ Section assignment ############
model.HomogeneousSolidSection(name='AL_section', 
    material='Aluminium', thickness=None)
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']

current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
)
region = regionToolset.Region(cells=current_cells)
p.SectionAssignment(region=region, sectionName='AL_section', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

############ Partitioning ############

#### Strat - 1 ####
'''
n_long = 10  # number of longitudinal partitions
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
    p.PartitionCellByDatumPlane(datumPlane=dp_rot_obj, cells=current_cells)'''

#### Strat - 2 ####
'''p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
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
    p.PartitionCellByDatumPlane(datumPlane=dp_offset_obj, cells=current_cells)'''

#### Strat - 3 ####
