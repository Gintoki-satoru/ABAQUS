from abaqus import *
from abaqusConstants import *
import sketch
import os
import numpy as np

path_modules = 'N:\\Sachdeva\\MT_Nair\\ABAQUS\\MT\\Macros'

os.chdir(path_modules)

# Further packages:
import coordinateTransformation_ellipse as ct

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

m_a_inner = 150.0  # Semi-major axis of the inner ellipse
m_b_inner = 100.0  # Semi-minor axis of the inner ellipse

thick = 5.0  # Thickness

m_a_outer = m_a_inner + thick  # Semi-major axis of the outer ellipse
m_b_outer = m_b_inner + thick  # Semi-minor axis of the outer ellipse

# Create a new model
myModel = mdb.Model(name='EllipseModel_2D')

# Create a new sketch
mySketch = myModel.ConstrainedSketch(name='EllipseSketch', sheetSize=200.0)

center = (0.0, 0.0)
majorAxisPoint = (m_a_inner, 0.0)
minorAxisPoint = (0.0, m_b_inner)

mySketch.EllipseByCenterPerimeter(center=center, axisPoint1=majorAxisPoint, axisPoint2=minorAxisPoint)

outerMajorAxisPoint = (m_a_outer, 0.0)
outerMinorAxisPoint = (0.0, m_b_outer)

mySketch.EllipseByCenterPerimeter(center=center, axisPoint1=outerMajorAxisPoint, axisPoint2=outerMinorAxisPoint)

mySketch.Line(point1=(0.0, m_b_outer), point2=(0.0, -m_b_outer))
mySketch.Line(point1=(-m_a_outer, 0.0), point2=(m_a_outer, 0.0))

# Use the coordinate transformation functions to find points on the curve to trim
phi_values = np.radians([135, 225, 315])
a_inner, b_inner = m_a_inner, m_b_inner
a_outer, b_outer = m_a_outer, m_b_outer
g = mySketch.geometry

for phi in phi_values:
    x_inner = ct.pol2cart_x(a_inner, phi)
    y_inner = ct.pol2cart_y(b_inner, phi)
    x_outer = ct.pol2cart_x(a_outer, phi)
    y_outer = ct.pol2cart_y(b_outer, phi)
    
    try:
        mySketch.autoTrimCurve(curve1=g.findAt((x_inner, y_inner)), point1=(x_inner, y_inner))
        mySketch.autoTrimCurve(curve1=g.findAt((x_outer, y_outer)), point1=(x_outer, y_outer))
    except:
        pass

mySketch.autoTrimCurve(curve1=g.findAt((0, m_b_inner/2)), point1=(0, m_b_inner/2))
mySketch.autoTrimCurve(curve1=g.findAt((0, -m_b_inner/2)), point1=(0, -m_b_inner/2))
mySketch.autoTrimCurve(curve1=g.findAt((m_a_inner/2, 0)), point1=(m_a_inner/2, 0))


myPart = myModel.Part(name='EllipsePart', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
myPart.BaseShell(sketch=mySketch)

# Material properties

myModel.Material(name='Aluminium')
myModel.materials['Aluminium'].Elastic(table=((70000.0,0.35), ))


