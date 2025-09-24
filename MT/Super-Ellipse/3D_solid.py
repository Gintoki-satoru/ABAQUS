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

# ------------------------
# Parameters
# ------------------------
a = 150.0
b = 100.0
n = 2.0
thick = 2.5
num_points = 100

################### Create model ###################
Mdb()
modelName = 'SuperEllipse'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
model = mdb.models[modelName]
model.rootAssembly.clearGeometryCache()
model.rootAssembly.regenerate()

################### Create sketch ###################
s = model.ConstrainedSketch(name='__profile__', sheetSize=2*max(a+thick,b+thick))

axis_line = s.Line(point1=(0.0, 0.0), point2=(0.0, b+thick))
s.setAsConstruction(objectList=(axis_line,))

# --- Inner profile ---
inner_points = []
for i in range(num_points + 1):
    y = b * i / num_points
    x = a * (1 - (y / b) ** n) ** (1.0 / n)
    inner_points.append((x, y))

s.Spline(points=inner_points)

# --- Outer profile ---
outer_points = []
for i in range(num_points + 1):
    y = (b + thick) * i / num_points
    x = (a + thick) * (1 - (y / (b + thick)) ** n) ** (1.0 / n)
    outer_points.append((x, y))

outer_points.reverse()
s.Spline(points=outer_points)

s.Line(point1=inner_points[-1], point2=outer_points[0])
s.Line(point1=outer_points[-1], point2=inner_points[0])

################### Create solid part ###################
part_name = 'Dome'
myPart = model.Part(name=part_name, dimensionality=THREE_D,
                    type=DEFORMABLE_BODY)

myPart.BaseSolidRevolve(sketch=s, angle=90.0, flipRevolveDirection=ON)
s.unsetPrimaryObject()

############ Function for selection in parametric form ############

def superellipse_point(u):
    """Parametric 2D superellipse curve point."""
    cos_u = math.cos(u)
    sin_u = math.sin(u)
    x = a * math.copysign(abs(cos_u)**(2.0/n), cos_u)
    y = b * math.copysign(abs(sin_u)**(2.0/n), sin_u)
    return x, y

def revolved_surface_point(phi, theta):
    """Returns 3D point on revolved surface."""
    x_curve, y_curve = superellipse_point(phi)
    x = x_curve * math.cos(theta)
    y = y_curve
    z = x_curve * math.sin(theta)
    return (x, y, z)

############ Material properties ############
model.Material(name='Aluminium')
model.materials['Aluminium'].Elastic(table=((70000.0,0.34), ))

############ Section assignment ############
mdb.models['SuperEllipse'].HomogeneousSolidSection(name='AL_section', 
    material='Aluminium', thickness=None)

p = mdb.models['SuperEllipse'].parts['Dome']
long = math.pi / 4      # longitude
lat = -math.pi / 4  # latitude
x, y, z = revolved_surface_point(long, lat)
picked_point = (x, y, z)
c = p.cells
cells = c.findAt(((x,y,z), ))
region = regionToolset.Region(cells=cells)
p = mdb.models['SuperEllipse'].parts['Dome']
p.SectionAssignment(region=region, sectionName='AL_section', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

########### Assembly ##################
a = mdb.models['SuperEllipse'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['SuperEllipse'].parts['Dome']
a.Instance(name='Dome-1', part=p, dependent=OFF)
session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)

# --------------------------
# Parameters
# --------------------------
n_partitions = 4  # number of longitudinal partitions
p = mdb.models['SuperEllipse'].parts['Dome']

# Optional: delete previous datum planes / partitions
p.deleteFeatures(('Datum plane-1', 'Datum plane-2', 'Partition cell-1'))

# --- Step 1: create datum plane at 0 offset ---
dp = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
datumPlane0 = p.datums[dp.id]

# --- Step 2: identify construction line (rotation axis) ---
# Usually first vertical construction line along y-axis in part sketch
axis_line = p.datums[1]  # adjust index if needed

# --- Step 3: loop to create rotated datum planes ---
angle_increment = 90.0 / n_partitions
rotated_planes = []

for i in range(1, n_partitions):
    angle = i * angle_increment
    dp_rot = p.DatumPlaneByRotation(plane=datumPlane0, axis=axis_line, angle=angle)
    rotated_planes.append(p.datums[dp_rot.id])

# --- Step 4: partition the cell using all rotated planes ---
cell = p.cells[0]  # assuming only one solid dome cell
for dp in rotated_planes:
    p.PartitionCellByDatumPlane(datumPlane=dp, cells=(cell,))

