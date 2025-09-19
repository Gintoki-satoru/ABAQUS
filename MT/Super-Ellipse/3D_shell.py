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
n = 3
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

############ Function for coordinate transform ############
def ct_cart(a, b, c, r, s, u, v):
    # a,b,c superellipsoid radii
    # r,s shape exponents
    # u latitude angle
    # v longitude angle
    theta = np.radians(u)  # latitude
    phi = np.radians(v)
    eps1 = 2.0 / n
    eps2 = 2.0 / n
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    x = a * np.sign(cos_theta) * abs(cos_theta)**eps1 * np.sign(cos_phi) * abs(cos_phi)**eps2
    y = b * np.sign(cos_theta) * abs(cos_theta)**eps1 * np.sign(sin_phi) * abs(sin_phi)**eps2
    z = c * np.sign(sin_theta) * abs(sin_theta)**eps1
    return x, y, z

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
theta = -45
phi = 45
x, y, z = ct_cart(a, b, c, n, n, theta, phi)
faces = f.findAt(((x, y, z), ))
region = regionToolset.Region(faces=faces)
p = mdb.models['SuperEllipse'].parts['Dome']
p.SectionAssignment(region=region, sectionName='Dome_shell', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)