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

################### Parameters #####################
a = 150.0
b = 100.0
c = a
n = 2
thick = 2.5
num_points = 100

################### Create model ###################
Mdb()
modelName = 'SuperEllipse'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
mdb.models[modelName].rootAssembly.clearGeometryCache()
mdb.models[modelName].rootAssembly.regenerate()
model = mdb.models[modelName]

################### Create sketch ###################
points = []
for i in range(num_points + 1):
    y = b * i / num_points
    x = a * (1 - (y / b) ** n) ** (1.0 / n)
    points.append((x, y))

s = model.ConstrainedSketch(name='__profile__', sheetSize=2*max(a,b))

axis_line = s.Line(point1=(0.0, 0.0), point2=(0.0, b))

s.setAsConstruction(objectList=(axis_line,))

s.Spline(points=points)

part_name = 'Dome'
myPart = model.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

myPart.BaseShellRevolve(sketch=s, angle=90.0, flipRevolveDirection=ON)
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
model.HomogeneousShellSection(name='Dome_shell', 
    preIntegrate=OFF, material='Aluminium', thicknessType=UNIFORM, 
    thickness=thick, thicknessField='', nodalThicknessField='', 
    idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
    thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
    integrationRule=SIMPSON, numIntPts=5)

p = mdb.models['SuperEllipse'].parts['Dome']
f = p.faces
long = math.pi / 4      # longitude
lat = -math.pi / 4  # latitude
x, y, z = revolved_surface_point(long, lat)
picked_point = (x, y, z)

faces = f.findAt(((x, y, z), ))
region = regionToolset.Region(faces=faces)
p = mdb.models['SuperEllipse'].parts['Dome']
p.SectionAssignment(region=region, sectionName='Dome_shell', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

