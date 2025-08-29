# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

# Vertical Strut Creation
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-2.0, -2.0), point2=
    (2.02499999995343, -2.0))
mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
    addUndoState=False, entity=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2.02499999995343, 
    -2.0), point2=(2.02499999995344, 2.0))
mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
mdb.models['Model-1'].sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], entity2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2.02499999995344, 
    2.0), point2=(-2.0, -2.0))
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    0.209346607327461, -2.21029210090637), value=4.0, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[0], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[1])
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    2.46034598350525, -0.514440298080444), value=4.0, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[1], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2])
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    -0.532000243663788, -0.0635590553283691), value=5.65685424949239, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[0])
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], ))
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], ))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseWire(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

# Body strut creation
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=6.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-3.0, -3.0), point2=
    (3.0, -3.0))
mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
    addUndoState=False, entity=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(3.0, -3.0), point2=(
    3.0, 3.0))
mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
mdb.models['Model-1'].sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], entity2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(3.0, 3.0), point2=(
    -3.0, -3.0))
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    3.78300213813782, -0.497716873884201), value=4.0, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[1], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2])
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    0.790014147758484, -2.50974082946777), value=5.65685424949238, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[0], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[1])
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    -0.0277883037924767, -0.261998385190964), value=6.92820323027551, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[0])
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], ))
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], ))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='body_strut', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['body_strut'].BaseWire(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

# Material Properties
mdb.models['Model-1'].Material(name='Aluminium')
mdb.models['Model-1'].materials['Aluminium'].Elastic(table=((70000.0, 0.35), ))
mdb.models['Model-1'].Material(name='cfrp')
mdb.models['Model-1'].materials['cfrp'].Elastic(table=((150000.0, 9000.0, 
    9000.0, 0.34, 0.34, 0.4, 5000.0, 5000.0, 5000.0), ), type=
    ENGINEERING_CONSTANTS)

# Section assignment   
mdb.models['Model-1'].CircularProfile(name='Profile-1', r=0.2)
mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
    DURING_ANALYSIS, material='Aluminium', name='core_strut', poissonRatio=0.0, 
    profile='Profile-1', temperatureVar=LINEAR)
mdb.models['Model-1'].parts['body_strut'].Set(edges=
    mdb.models['Model-1'].parts['body_strut'].edges.getSequenceFromMask((
    '[#1 ]', ), ), name='Set-1')
mdb.models['Model-1'].parts['body_strut'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['body_strut'].sets['Set-1'], sectionName=
    'core_strut', thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].parts['Part-1'].Set(edges=
    mdb.models['Model-1'].parts['Part-1'].edges.getSequenceFromMask(('[#1 ]', 
    ), ), name='Set-1')
mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-1'].sets['Set-1'], sectionName=
    'core_strut', thicknessAssignment=FROM_SECTION)

# Unit cell creation
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-1', ), 
    vector=(1.975, 2.0, 0.0))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-2', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 
    10.0, 0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-1-1', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-1', ), 
    vector=(3.975, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-2', ), 
    vector=(1.975, 2.0, 0.0))
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(0.0, 0.0, 1.0), instanceList=('Part-1-1', 'Part-1-2'), 
    number1=1, number2=2, spacing1=4.025, spacing2=4.0)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-3', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-3', ), 
    vector=(1.975, 2.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=-90.0, axisDirection=(0.0, 
    10.0, 0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-1-3', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-4', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-4', ), 
    vector=(1.975, 2.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 10.0, 
    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-1-4', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-4', ), 
    vector=(0.0, 0.0, 4.0))
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(1.0, 0.0, 0.0), instanceList=('Part-1-3', 'Part-1-4'), 
    number1=2, number2=1, spacing1=4.0, spacing2=4.0)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-1', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('body_strut-1', ), 
    vector=(2.656854, 3.0, 0.0))
mdb.models['Model-1'].rootAssembly.deleteFeatures(('Part-1-1', 'Part-1-2', 
    'Part-1-1-lin-1-2', 'Part-1-2-lin-1-2', 'Part-1-3', 'Part-1-4', 
    'Part-1-3-lin-2-1', 'Part-1-4-lin-2-1', 'body_strut-1'))
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=
    mdb.models['Model-1'].parts['Part-1'].features['Wire-1'].sketch)
mdb.models['Model-1'].parts['Part-1'].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__edit__'], 
    upToFeature=mdb.models['Model-1'].parts['Part-1'].features['Wire-1'])
mdb.models['Model-1'].sketches['__edit__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__edit__'].geometry[4], ))
mdb.models['Model-1'].sketches['__edit__'].Line(point1=(-2.0, 2.0), point2=(
    -2.0, -2.0))
mdb.models['Model-1'].sketches['__edit__'].VerticalConstraint(addUndoState=
    False, entity=mdb.models['Model-1'].sketches['__edit__'].geometry[5])
mdb.models['Model-1'].sketches['__edit__'].Line(point1=(-2.0, -2.0), point2=(
    2.0, -2.0))
mdb.models['Model-1'].sketches['__edit__'].HorizontalConstraint(addUndoState=
    False, entity=mdb.models['Model-1'].sketches['__edit__'].geometry[6])
mdb.models['Model-1'].sketches['__edit__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    mdb.models['Model-1'].sketches['__edit__'].geometry[5], entity2=
    mdb.models['Model-1'].sketches['__edit__'].geometry[6])
mdb.models['Model-1'].sketches['__edit__'].Line(point1=(2.0, -2.0), point2=(
    -2.0, 2.0))
mdb.models['Model-1'].sketches['__edit__'].ObliqueDimension(textPoint=(
    -3.08787965774536, 0.0663106441497803), value=4.0, vertex1=
    mdb.models['Model-1'].sketches['__edit__'].vertices[3], vertex2=
    mdb.models['Model-1'].sketches['__edit__'].vertices[4])
mdb.models['Model-1'].sketches['__edit__'].ObliqueDimension(textPoint=(
    0.106067657470703, -2.30174207687378), value=4.0, vertex1=
    mdb.models['Model-1'].sketches['__edit__'].vertices[4], vertex2=
    mdb.models['Model-1'].sketches['__edit__'].vertices[5])
mdb.models['Model-1'].sketches['__edit__'].ObliqueDimension(textPoint=(
    1.20577907562256, -0.50660514831543), value=5.65685424949238, vertex1=
    mdb.models['Model-1'].sketches['__edit__'].vertices[5], vertex2=
    mdb.models['Model-1'].sketches['__edit__'].vertices[3])
mdb.models['Model-1'].sketches['__edit__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__edit__'].geometry[5], ))
mdb.models['Model-1'].sketches['__edit__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__edit__'].geometry[6], ))
mdb.models['Model-1'].parts['Part-1'].features['Wire-1'].setValues(sketch=
    mdb.models['Model-1'].sketches['__edit__'])
del mdb.models['Model-1'].sketches['__edit__']
mdb.models['Model-1'].parts['Part-1'].regenerate()
mdb.models['Model-1'].parts.changeKey(fromName='Part-1', toName='face_strut')
mdb.models['Model-1'].parts['face_strut'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    edges=mdb.models['Model-1'].parts['face_strut'].edges.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName='core_strut', thicknessAssignment=
    FROM_SECTION)
del mdb.models['Model-1'].parts['face_strut'].sectionAssignments[0]
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='face_strut-1', 
    part=mdb.models['Model-1'].parts['face_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_strut-1', ), 
    vector=(-2.0, 2.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_strut-1', ), 
    vector=(4.0, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='face_strut-2', 
    part=mdb.models['Model-1'].parts['face_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_strut-2', ), 
    vector=(-2.0, 2.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 
    10.0, 0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('face_strut-2', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='face_strut-3', 
    part=mdb.models['Model-1'].parts['face_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_strut-3', ), 
    vector=(-2.0, 2.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=-90.0, axisDirection=(0.0, 
    10.0, 0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('face_strut-3', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_strut-3', ), 
    vector=(0.0, 0.0, 4.0))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='face_strut-4', 
    part=mdb.models['Model-1'].parts['face_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_strut-4', ), 
    vector=(-2.0, 2.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 10.0, 
    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('face_strut-4', ))
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(0.0, 0.0, 1.0), instanceList=('face_strut-1', 
    'face_strut-2'), number1=1, number2=2, spacing1=4.0, spacing2=4.0)
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(0.0, 1.0, 0.0), instanceList=('face_strut-3', 
    'face_strut-4'), number1=2, number2=1, spacing1=4.0, spacing2=4.0)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-1', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('body_strut-1', ), 
    vector=(2.656854, 3.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=-45.0, axisDirection=(0.0, 
    10.0, 0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('body_strut-1', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-2', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('body_strut-2', ), 
    vector=(2.656854, 3.0, 4.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=45.0, axisDirection=(0.0, 4.0, 
    0.0), axisPoint=(0.0, 0.0, 4.0), instanceList=('body_strut-2', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-3', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('body_strut-3', ), 
    vector=(6.656854, 3.0, 4.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=135.0, axisDirection=(0.0, 4.0, 
    0.0), axisPoint=(4.0, 0.0, 4.0), instanceList=('body_strut-3', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-4', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('body_strut-4', ), 
    vector=(6.656854, 3.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=-135.0, axisDirection=(0.0, 
    4.0, 0.0), axisPoint=(4.0, 0.0, 0.0), instanceList=('body_strut-4', ))
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['face_strut-1'], 
    mdb.models['Model-1'].rootAssembly.instances['face_strut-2'], 
    mdb.models['Model-1'].rootAssembly.instances['face_strut-3'], 
    mdb.models['Model-1'].rootAssembly.instances['face_strut-4'], 
    mdb.models['Model-1'].rootAssembly.instances['face_strut-1-lin-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['face_strut-2-lin-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['face_strut-3-lin-2-1'], 
    mdb.models['Model-1'].rootAssembly.instances['face_strut-4-lin-2-1'], 
    mdb.models['Model-1'].rootAssembly.instances['body_strut-1'], 
    mdb.models['Model-1'].rootAssembly.instances['body_strut-2'], 
    mdb.models['Model-1'].rootAssembly.instances['body_strut-3'], 
    mdb.models['Model-1'].rootAssembly.instances['body_strut-4']), name=
    'fbcc_unitcell', originalInstances=DELETE)
mdb.models['Model-1'].parts['fbcc_unitcell'].assignBeamSectionOrientation(
    method=N1_COSINES, n1=(0.0, 0.0, -1.0), region=Region(
    edges=mdb.models['Model-1'].parts['fbcc_unitcell'].edges.getSequenceFromMask(
    mask=('[#ffffff ]', ), )))

# Core - 4mm creation
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(0.0, 0.0, 
    1.0), direction2=(1.0, 0.0, 0.0), instanceList=('fbcc_unitcell-1', ), 
    number1=6, number2=50, spacing1=4.0, spacing2=4.0)
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-3'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-4'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-5'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-6'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-7'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-8'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-9'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-10'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-11'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-12'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-13'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-14'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-15'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-16'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-17'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-18'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-19'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-20'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-21'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-22'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-23'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-24'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-25'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-26'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-27'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-28'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-29'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-30'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-31'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-32'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-33'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-34'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-35'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-36'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-37'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-38'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-39'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-40'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-41'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-42'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-43'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-44'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-45'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-46'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-47'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-48'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-49'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-1-50'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-1'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-2'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-3'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-4'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-5'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-6'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-7'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-8'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-9'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-10'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-11'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-12'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-13'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-14'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-15'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-16'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-17'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-18'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-19'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-20'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-21'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-22'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-23'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-24'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-25'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-26'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-27'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-28'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-29'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-30'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-31'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-32'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-33'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-34'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-35'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-36'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-37'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-38'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-39'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-40'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-41'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-42'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-43'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-44'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-45'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-46'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-47'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-48'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-49'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-2-50'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-1'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-2'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-3'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-4'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-5'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-6'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-7'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-8'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-9'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-10'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-11'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-12'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-13'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-14'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-15'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-16'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-17'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-18'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-19'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-20'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-21'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-22'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-23'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-24'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-25'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-26'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-27'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-28'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-29'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-30'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-31'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-32'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-33'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-34'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-35'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-36'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-37'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-38'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-39'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-40'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-41'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-42'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-43'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-44'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-45'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-46'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-47'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-48'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-49'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-3-50'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-1'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-2'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-3'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-4'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-5'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-6'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-7'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-8'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-9'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-10'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-11'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-12'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-13'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-14'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-15'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-16'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-17'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-18'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-19'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-20'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-21'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-22'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-23'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-24'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-25'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-26'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-27'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-28'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-29'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-30'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-31'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-32'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-33'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-34'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-35'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-36'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-37'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-38'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-39'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-40'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-41'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-42'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-43'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-44'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-45'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-46'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-47'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-48'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-49'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-4-50'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-1'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-2'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-3'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-4'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-5'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-6'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-7'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-8'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-9'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-10'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-11'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-12'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-13'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-14'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-15'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-16'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-17'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-18'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-19'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-20'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-21'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-22'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-23'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-24'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-25'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-26'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-27'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-28'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-29'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-30'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-31'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-32'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-33'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-34'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-35'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-36'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-37'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-38'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-39'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-40'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-41'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-42'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-43'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-44'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-45'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-46'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-47'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-48'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-49'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-5-50'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-1'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-2'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-3'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-4'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-5'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-6'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-7'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-8'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-9'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-10'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-11'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-12'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-13'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-14'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-15'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-16'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-17'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-18'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-19'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-20'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-21'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-22'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-23'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-24'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-25'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-26'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-27'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-28'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-29'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-30'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-31'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-32'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-33'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-34'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-35'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-36'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-37'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-38'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-39'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-40'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-41'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-42'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-43'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-44'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-45'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-46'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-47'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-48'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-49'], 
    mdb.models['Model-1'].rootAssembly.instances['fbcc_unitcell-1-lin-6-50']), 
    name='core_4mm', originalInstances=DELETE)
mdb.models['Model-1'].parts['core_4mm'].assignBeamSectionOrientation(method=
    N1_COSINES, n1=(0.0, 0.0, -1.0), region=Region(
    edges=mdb.models['Model-1'].parts['core_4mm'].edges.getSequenceFromMask(
    mask=('[#ffffffff:157 ]', ), )))

# Core 12mm creation
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(0.0, 1.0, 0.0), instanceList=('core_4mm-1', ), number1=1, 
    number2=3, spacing1=200.0, spacing2=4.0)
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['core_4mm-1'], 
    mdb.models['Model-1'].rootAssembly.instances['core_4mm-1-lin-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['core_4mm-1-lin-1-3']), name=
    'core_12mm', originalInstances=DELETE)
mdb.models['Model-1'].parts['core_12mm'].assignBeamSectionOrientation(method=
    N1_COSINES, n1=(0.0, 0.0, -1.0), region=Region(
    edges=mdb.models['Model-1'].parts['core_12mm'].edges.getSequenceFromMask(
    mask=('[#ffffffff:471 ]', ), )))

# Face sheet creation
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-20.0, -100.0), 
    point2=(28.75, 100.0))
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    59.8762817382812, 14.5294952392578), value=200.0, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[3])
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    15.25634765625, -110.374267578125), value=24.0, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[3], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[0])
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='face_sheet', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['face_sheet'].BaseShell(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

# Composite layup - Shell modelling
mdb.models['Model-1'].parts['face_sheet'].CompositeLayup(description='', 
    elementType=SHELL, name='CompositeLayup-1', offsetType=BOTTOM_SURFACE, 
    symmetric=False, thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].Section(
    integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
    temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].ReferenceOrientation(
    additionalRotationType=ROTATION_NONE, angle=0.0, axis=AXIS_3, fieldName='', 
    localCsys=None, orientationType=GLOBAL)
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].suppress(
    )
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].CompositePly(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, material='cfrp', numIntPoints=3, orientationType=
    SPECIFY_ORIENT, orientationValue=0.0, plyName='Ply-1', region=Region(
    faces=mdb.models['Model-1'].parts['face_sheet'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), suppressed=False, thickness=0.375, thicknessType=
    SPECIFY_THICKNESS)
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].CompositePly(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, material='cfrp', numIntPoints=3, orientationType=
    SPECIFY_ORIENT, orientationValue=90.0, plyName='Ply-2', region=Region(
    faces=mdb.models['Model-1'].parts['face_sheet'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), suppressed=False, thickness=0.375, thicknessType=
    SPECIFY_THICKNESS)
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].resume(
    )

# Sandwich Assembly
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='face_sheet-1', 
    part=mdb.models['Model-1'].parts['face_sheet'])
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 
    200.0, 0.0), axisPoint=(-20.0, -100.0, 0.0), instanceList=('face_sheet-1', 
    ))
mdb.models['Model-1'].rootAssembly.rotate(angle=-90.0, axisDirection=(0.0, 0.0, 
    -24.0), axisPoint=(-20.0, -100.0, 0.0), instanceList=('face_sheet-1', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_sheet-1', ), 
    vector=(220.0, 112.0, 24.0))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='face_sheet-2', 
    part=mdb.models['Model-1'].parts['face_sheet'])
mdb.models['Model-1'].rootAssembly.rotate(angle=-90.0, axisDirection=(0.0, 
    200.0, 0.0), axisPoint=(-20.0, -100.0, 0.0), instanceList=('face_sheet-2', 
    ))
mdb.models['Model-1'].rootAssembly.rotate(angle=-90.0, axisDirection=(0.0, 0.0, 
    -24.0), axisPoint=(-20.0, -100.0, 24.0), instanceList=('face_sheet-2', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_sheet-2', ), 
    vector=(220.0, 100.0, 0.0))

# Buckling step creation
mdb.models['Model-1'].BuckleStep(maxIterations=300, name='Step-1', numEigen=1, 
    previous='Initial', vectors=30)

# Reference point
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0.0, 6.0, 12.0))
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(200.0, 6.0, 12.0))

# Constraint
mdb.models['Model-1'].Tie(adjust=ON, main=Region(
    side2Faces=mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), name='Constraint-1', positionToleranceMethod=COMPUTED
    , secondary=Region(
    vertices=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].vertices.getSequenceFromMask(
    mask=('[#1050 #4050440 #10040028 #0 #4 #24800 #8000000', 
    ' #34021 #20000000 #a #0:4 #800a208 #80000410 #28200808', 
    ' #88 #0 #4240 #0:38 #6880a0 #c48c0004 #41410100', 
    ' #2 #21 #0:7 #8200000 #80 #0:2 #448005', 
    ' #140000 #112 #4480050 #1400000 #1120 #44800500 #14000000', 
    ' #11200 #48005000 #40000004 #112001 #80050000 #44 #1120014', 
    ' #500000 #448 #11200140 #5000000 #4480 #12001400 #50000001', 
    ' #44800 #20014000 #11 #448005 #140000 #112 #90028', 
    ' #44000000 #80405144 #8c214000 #84 #b00000 #90448401 #4429040', 
    ' #110a0209 #44280824 #10a02090 #42808241 #a020904 #28082411 #a0209044', 
    ' #80824110 #2090442 #824110a #20904428 #824110a0 #9044280 #24110a02', 
    ' #90442808 #4110a020 #4428082 #21820209 #48042290 #90482412 #4824120', 
    ' #48241209 #82412090 #24120904 #8 ]', ), )), thickness=ON, tieRotations=
    ON)
mdb.models['Model-1'].constraints.changeKey(fromName='Constraint-1', toName=
    'top')
mdb.models['Model-1'].Tie(adjust=ON, main=Region(
    side2Faces=mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), name='bottom', positionToleranceMethod=COMPUTED, 
    secondary=Region(
    vertices=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].vertices.getSequenceFromMask(
    mask=('[#0 #5 #0 #8a00a000 #22011140 #1a #a00000', 
    ' #2800000 #22008a0 #c4804180 #50720502 #55555555:2 #1155555 #400000', 
    ' #8200000 #0:2 #12801421 #55540000 #55555555:2 #881155 #1480', 
    ' #14104041 #4080000 #14104 #10404080 #8000014 #1410404 #40408000', 
    ' #1410 #41040408 #40800001 #141040 #4040800 #80000141 #14104040', 
    ' #4080000 #14104 #10404080 #8000014 #1410404 #40408000 #1410', 
    ' #41040408 #40800001 #141040 #4040800 #80000141 #14104040 #4080000', 
    ' #14104 #10404080 #8000014 #1410404 #8102000 #4000051 #0:3', 
    ' #8081000 #a285200 #1428a145 #1450a285 #2851428a #28a1450a #50a28514', 
    ' #51428a14 #61450a28 #2214 #20 ]', ), )), thickness=ON, tieRotations=ON)
mdb.models['Model-1'].Coupling(alpha=0.0, controlPoint=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[656], )), couplingType=
    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name='left', 
    surface=Region(
    edges=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].edges.getSequenceFromMask(
    mask=('[#0:26 #78000000 #6 #0:19 #19e0000 #0:50 #6033c000', 
    ' #66666 #bcc00 #0:34 #c0000000 #c0000001 #818 #0:118', 
    ' #40000000 #1800010 #c #833 #0:33 #100 #1824', 
    ' #0:97 #80000 #4800 #0:7 #f20c0000 #2 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].edges.getSequenceFromMask(
    mask=('[#4 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].edges.getSequenceFromMask(
    mask=('[#4 ]', ), ), 
    vertices=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].vertices.getSequenceFromMask(
    mask=('[#0:13 #37d00000 #0:8 #37d00 #0:34 #dbbadc0 #e0', 
    ' #0:11 #2000 #80 #0:32 #1c0 #0:2 #21800000', ' #0:20 #200050 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].vertices.getSequenceFromMask(
    mask=('[#c ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].vertices.getSequenceFromMask(
    mask=('[#c ]', ), )), u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
mdb.models['Model-1'].Coupling(alpha=0.0, controlPoint=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[657], )), couplingType=
    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name='right', 
    surface=Region(
    edges=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].edges.getSequenceFromMask(
    mask=('[#3000000 #6000c000 #180 #0:2 #80000000 #303661 #180', 
    ' #0 #658 #0 #6060000 #3280000 #0:5 #30', 
    ' #0:2 #8070040 #186181 #400 #0:3 #c000000 #2', 
    ' #6c0 #0:3 #4000 #0:9 #f000000 #0:5 #8000000', 
    ' #4ce #0:53 #40000000 #0:12 #70000000 #0:20 #4 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].edges.getSequenceFromMask(
    mask=('[#1 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].edges.getSequenceFromMask(
    mask=('[#1 ]', ), ), 
    vertices=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].vertices.getSequenceFromMask(
    mask=('[#c00001 #1c300740 #0:2 #18707c #300 #4000 #0', 
    ' #e0000000 #40000002 #22b1c2 #0:3 #60c0000 #0:4 #8000', 
    ' #0:2 #20000000 #0:39 #20 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].vertices.getSequenceFromMask(
    mask=('[#3 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].vertices.getSequenceFromMask(
    mask=('[#3 ]', ), )), u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

# Boundary Conditions
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
    region=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[656], )), u1=SET, u2=SET
    , u3=SET, ur1=SET, ur2=SET, ur3=SET)
mdb.models['Model-1'].boundaryConditions.changeKey(fromName='BC-1', toName=
    'left')
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='right', 
    region=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[657], )), u1=UNSET, u2=
    SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)

# Load
mdb.models['Model-1'].ConcentratedForce(cf1=-1.0, createStepName='Step-1', 
    distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
    Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[657], )))

# Mesh 
mdb.models['Model-1'].parts['face_sheet'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=0.5)
mdb.models['Model-1'].parts['face_sheet'].generateMesh()
mdb.models['Model-1'].parts['face_sheet'].setElementType(elemTypes=(ElemType(
    elemCode=S8R, elemLibrary=STANDARD), ElemType(elemCode=STRI65, 
    elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['face_sheet'].faces.getSequenceFromMask((
    '[#1 ]', ), ), ))
mdb.models['Model-1'].parts['core_12mm'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=1.0)
mdb.models['Model-1'].parts['core_12mm'].generateMesh()
mdb.models['Model-1'].rootAssembly.regenerate()

# Job creation
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='12_shell', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, numThreadsPerMpiProcess=1, queue=None, resultsFormat=
    ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, 
    waitMinutes=0)
# Save by anair on 2025_01_27-16.25.00; build 2023 2022_09_28-20.11.55 183150
