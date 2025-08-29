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
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
    0.0, 4.0))
mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='vert_strut', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['vert_strut'].BaseWire(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=6.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 4.0), point2=(
    0.0, 0.0))
mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
    2.95, 0.0))
mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
    addUndoState=False, entity=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3])
mdb.models['Model-1'].sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], entity2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2.95, 0.0), point2=(
    0.0, 4.0))
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    0.918419063091278, -0.456526756286621), value=5.65685424949238, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[1], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2])
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    3.96096634864807, 2.23051881790161), value=6.92820323027551, vertex1=
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
mdb.models['Model-1'].Material(name='Aluminium')
mdb.models['Model-1'].materials['Aluminium'].Elastic(table=((70000.0, 0.35), ))
mdb.models['Model-1'].Material(name='carbon_epoxy')
mdb.models['Model-1'].materials['carbon_epoxy'].Elastic(table=((125000.0, 
    5000.0, 5000.0, 0.25, 0.25, 0.4, 2500.0, 2500.0, 1800.0), ), type=
    ENGINEERING_CONSTANTS)
mdb.models['Model-1'].CircularProfile(name='circle', r=0.2)
mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
    DURING_ANALYSIS, material='Aluminium', name='core_strut', poissonRatio=0.0, 
    profile='circle', temperatureVar=LINEAR)
mdb.models['Model-1'].parts['body_strut'].Set(edges=
    mdb.models['Model-1'].parts['body_strut'].edges.getSequenceFromMask((
    '[#1 ]', ), ), name='Set-1')
mdb.models['Model-1'].parts['body_strut'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['body_strut'].sets['Set-1'], sectionName=
    'core_strut', thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].parts['vert_strut'].Set(edges=
    mdb.models['Model-1'].parts['vert_strut'].edges.getSequenceFromMask((
    '[#1 ]', ), ), name='Set-1')
mdb.models['Model-1'].parts['vert_strut'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['vert_strut'].sets['Set-1'], sectionName=
    'core_strut', thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].parts['vert_strut'].Set(edges=
    mdb.models['Model-1'].parts['vert_strut'].edges.getSequenceFromMask((
    '[#1 ]', ), ), name='Set-2')
mdb.models['Model-1'].parts['vert_strut'].assignBeamSectionOrientation(method=
    N1_COSINES, n1=(0.0, 0.0, -1.0), region=
    mdb.models['Model-1'].parts['vert_strut'].sets['Set-2'])
mdb.models['Model-1'].parts['body_strut'].assignBeamSectionOrientation(method=
    N1_COSINES, n1=(0.0, 0.0, -1.0), region=Region(
    edges=mdb.models['Model-1'].parts['body_strut'].edges.getSequenceFromMask(
    mask=('[#1 ]', ), )))
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='vert_strut-1', 
    part=mdb.models['Model-1'].parts['vert_strut'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='vert_strut-2', 
    part=mdb.models['Model-1'].parts['vert_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('vert_strut-2', ), 
    vector=(4.0, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(0.0, 0.0, 1.0), instanceList=('vert_strut-1', 
    'vert_strut-2'), number1=1, number2=2, spacing1=4.0, spacing2=4.0)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-1', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.rotate(angle=-45.0, axisDirection=(0.0, 4.0, 
    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('body_strut-1', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-2', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('body_strut-2', ), 
    vector=(0.0, 0.0, 4.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=45.0, axisDirection=(0.0, 4.0, 
    0.0), axisPoint=(0.0, 0.0, 4.0), instanceList=('body_strut-2', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-3', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('body_strut-3', ), 
    vector=(4.0, 0.0, 4.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=135.0, axisDirection=(0.0, 4.0, 
    0.0), axisPoint=(4.0, 0.0, 4.0), instanceList=('body_strut-3', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='body_strut-4', 
    part=mdb.models['Model-1'].parts['body_strut'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('body_strut-4', ), 
    vector=(4.0, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=-135.0, axisDirection=(0.0, 
    4.0, 0.0), axisPoint=(4.0, 0.0, 0.0), instanceList=('body_strut-4', ))
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['vert_strut-1'], 
    mdb.models['Model-1'].rootAssembly.instances['vert_strut-2'], 
    mdb.models['Model-1'].rootAssembly.instances['vert_strut-1-lin-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['vert_strut-2-lin-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['body_strut-1'], 
    mdb.models['Model-1'].rootAssembly.instances['body_strut-2'], 
    mdb.models['Model-1'].rootAssembly.instances['body_strut-3'], 
    mdb.models['Model-1'].rootAssembly.instances['body_strut-4']), name=
    'bccz_unit_cell', originalInstances=DELETE)
mdb.models['Model-1'].parts['bccz_unit_cell'].assignBeamSectionOrientation(
    method=N1_COSINES, n1=(0.0, 0.0, -1.0), region=Region(
    edges=mdb.models['Model-1'].parts['bccz_unit_cell'].edges.getSequenceFromMask(
    mask=('[#fff ]', ), )))
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(0.0, 0.0, 1.0), instanceList=('bccz_unit_cell-1', ), 
    number1=50, number2=6, spacing1=4.0, spacing2=4.0)
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-1-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-1-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-1-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-1-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-2-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-2-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-2-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-2-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-2-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-2-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-3-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-3-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-3-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-3-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-3-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-3-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-4-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-4-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-4-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-4-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-4-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-4-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-5-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-5-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-5-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-5-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-5-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-5-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-6-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-6-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-6-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-6-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-6-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-6-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-7-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-7-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-7-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-7-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-7-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-7-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-8-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-8-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-8-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-8-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-8-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-8-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-9-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-9-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-9-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-9-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-9-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-9-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-10-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-10-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-10-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-10-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-10-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-10-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-11-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-11-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-11-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-11-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-11-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-11-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-12-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-12-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-12-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-12-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-12-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-12-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-13-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-13-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-13-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-13-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-13-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-13-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-14-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-14-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-14-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-14-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-14-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-14-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-15-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-15-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-15-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-15-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-15-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-15-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-16-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-16-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-16-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-16-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-16-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-16-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-17-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-17-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-17-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-17-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-17-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-17-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-18-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-18-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-18-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-18-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-18-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-18-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-19-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-19-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-19-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-19-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-19-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-19-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-20-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-20-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-20-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-20-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-20-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-20-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-21-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-21-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-21-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-21-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-21-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-21-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-22-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-22-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-22-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-22-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-22-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-22-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-23-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-23-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-23-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-23-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-23-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-23-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-24-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-24-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-24-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-24-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-24-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-24-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-25-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-25-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-25-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-25-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-25-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-25-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-26-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-26-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-26-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-26-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-26-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-26-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-27-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-27-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-27-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-27-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-27-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-27-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-28-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-28-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-28-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-28-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-28-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-28-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-29-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-29-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-29-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-29-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-29-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-29-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-30-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-30-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-30-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-30-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-30-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-30-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-31-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-31-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-31-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-31-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-31-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-31-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-32-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-32-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-32-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-32-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-32-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-32-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-33-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-33-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-33-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-33-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-33-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-33-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-34-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-34-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-34-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-34-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-34-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-34-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-35-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-35-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-35-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-35-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-35-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-35-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-36-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-36-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-36-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-36-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-36-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-36-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-37-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-37-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-37-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-37-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-37-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-37-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-38-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-38-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-38-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-38-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-38-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-38-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-39-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-39-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-39-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-39-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-39-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-39-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-40-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-40-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-40-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-40-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-40-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-40-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-41-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-41-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-41-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-41-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-41-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-41-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-42-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-42-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-42-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-42-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-42-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-42-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-43-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-43-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-43-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-43-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-43-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-43-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-44-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-44-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-44-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-44-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-44-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-44-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-45-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-45-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-45-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-45-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-45-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-45-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-46-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-46-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-46-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-46-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-46-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-46-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-47-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-47-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-47-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-47-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-47-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-47-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-48-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-48-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-48-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-48-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-48-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-48-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-49-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-49-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-49-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-49-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-49-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-49-6'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-50-1'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-50-2'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-50-3'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-50-4'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-50-5'], 
    mdb.models['Model-1'].rootAssembly.instances['bccz_unit_cell-1-lin-50-6']), 
    name='core_4mm', originalInstances=DELETE)
mdb.models['Model-1'].parts['core_4mm'].assignBeamSectionOrientation(method=
    N1_COSINES, n1=(0.0, 0.0, -1.0), region=Region(
    edges=mdb.models['Model-1'].parts['core_4mm'].edges.getSequenceFromMask(
    mask=('[#ffffffff:86 #1f ]', ), )))
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
    mask=('[#ffffffff:258 #7fff ]', ), )))
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-100.0, 25.0), 
    point2=(100.0, -20.0))
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    -16.645133972168, -30.2778625488281), value=200.0, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[1], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2])
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    112.591865539551, -2.12805891036987), value=24.0, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[3])
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='face_sheet', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['face_sheet'].BaseShell(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
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
    , axis=AXIS_3, material='carbon_epoxy', numIntPoints=27, orientationType=
    SPECIFY_ORIENT, orientationValue=0.0, plyName='Ply-1', region=Region(
    faces=mdb.models['Model-1'].parts['face_sheet'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), suppressed=False, thickness=0.375, thicknessType=
    SPECIFY_THICKNESS)
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].CompositePly(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, material='carbon_epoxy', numIntPoints=27, orientationType=
    SPECIFY_ORIENT, orientationValue=90.0, plyName='Ply-2', region=Region(
    faces=mdb.models['Model-1'].parts['face_sheet'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), suppressed=False, thickness=0.375, thicknessType=
    SPECIFY_THICKNESS)
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].CompositePly(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, material='carbon_epoxy', numIntPoints=27, orientationType=
    SPECIFY_ORIENT, orientationValue=0.0, plyName='Ply-3', region=Region(
    faces=mdb.models['Model-1'].parts['face_sheet'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), suppressed=False, thickness=0.375, thicknessType=
    SPECIFY_THICKNESS)
mdb.models['Model-1'].parts['face_sheet'].compositeLayups['CompositeLayup-1'].resume(
    )
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='face_sheet-1', 
    part=mdb.models['Model-1'].parts['face_sheet'])
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(-200.0, 
    0.0, 0.0), axisPoint=(100.0, -20.0, 0.0), instanceList=('face_sheet-1', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_sheet-1', ), 
    vector=(100.0, 32.0, 24.0))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='face_sheet-2', 
    part=mdb.models['Model-1'].parts['face_sheet'])
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(200.0, 
    0.0, 0.0), axisPoint=(-100.0, -20.0, 0.0), instanceList=('face_sheet-2', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_sheet-2', ), 
    vector=(100.0, 20.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=-180.0, axisDirection=(200.0, 
    0.0, 0.0), axisPoint=(0.0, 12.0, 0.0), instanceList=('face_sheet-1', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_sheet-1', ), 
    vector=(0.0, 0.0, 24.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_sheet-1', ), 
    vector=(0.0, 8.0, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=180.0, axisDirection=(200.0, 
    0.0, 0.0), axisPoint=(0.0, 20.0, 0.0), instanceList=('face_sheet-1', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('face_sheet-1', ), 
    vector=(0.0, -8.0, 24.0))
mdb.models['Model-1'].BuckleStep(maxIterations=300, name='Step-1', numEigen=1, 
    previous='Initial', vectors=30)
mdb.models['Model-1'].Tie(adjust=ON, main=Region(
    side2Faces=mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), name='top_face', positionToleranceMethod=COMPUTED, 
    secondary=Region(
    vertices=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].vertices.getSequenceFromMask(
    mask=('[#a1224514 #ecd #0 #89228000 #10001d80 #104 #24d00000', 
    ' #181 #a0100000 #48914 #62a5 #7c000000 #8a #0', 
    ' #1812bc #40640000 #0 #55db0000 #81002 #1000000 #501925a4', 
    ' #1000000 #48914a #252c0 #4dc10800 #11 #90000000 #1a012a', 
    ' #80e0000 #0 #55d2a000 #40802 #800000 #a034969 #80400000', 
    ' #122452 #94b0 #53704200 #4 #a4000000 #6804a #2038000', 
    ' #0 #9574a800 #10200 #40200000 #280d25a #a0100000 #48914', 
    ' #252c #14dc1080 #1 #a9000000 #1a012 #80e000 #0', 
    ' #255d2a00 #4080 #90080000 #a03496 #28040000 #12245 #94b', 
    ' #45370420 #0 #aa400000 #6804 #203800 #0 #9574a80', 
    ' #1020 #aa400000 #635072 ]', ), )), thickness=ON, tieRotations=ON)
mdb.models['Model-1'].Tie(adjust=ON, main=Region(
    side2Faces=mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), name='bottom_face', positionToleranceMethod=COMPUTED, 
    secondary=Region(
    vertices=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].vertices.getSequenceFromMask(
    mask=('[#0:2 #10422000 #2509 #4480000 #8b02800 #72969 #11c17000', 
    ' #22 #c9800000 #11248000 #4a10 #0 #a4244111 #44a00002', 
    ' #1590 #21160012 #20e52d #74800000 #40018884 #80000001 #884', 
    ' #88000000 #21250022 #5 #10000000 #4a424411 #4a00000 #40000a92', 
    ' #d2116002 #281e52 #38000000 #5000c442 #30000000 #221 #a2000000', 
    ' #48494008 #1 #44000000 #12909104 #81280000 #900002a4 #b4845800', 
    ' #a0794 #8e000000 #14003110 #4c000000 #88 #28800000 #52125002', 
    ' #0 #11000000 #4a42441 #204a0000 #240000a9 #2d211600 #281e5', 
    ' #23800000 #5000c44 #13000000 #22 #8a200000 #14849400 #0', 
    ' #44400000 #1290910 #48128000 #900002a #4b484580 #a079 #34e00000', 
    ' #10d44 #888200 ]', ), )), thickness=ON, tieRotations=ON)
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0.0, 6.0, 12.0))
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(200.0, 6.0, 12.0))
mdb.models['Model-1'].Coupling(alpha=0.0, controlPoint=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[634], )), couplingType=
    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name='left', 
    surface=Region(
    edges=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].edges.getSequenceFromMask(
    mask=('[#0:246 #20108000 #240001 #0:3 #800009 #4800000 #100001', 
    ' #0 #40000048 #0 #2200020 #440 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].edges.getSequenceFromMask(
    mask=('[#2 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].edges.getSequenceFromMask(
    mask=('[#2 ]', ), ), 
    vertices=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].vertices.getSequenceFromMask(
    mask=('[#0:70 #9e000000 #6003800d #df3f50 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].vertices.getSequenceFromMask(
    mask=('[#6 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].vertices.getSequenceFromMask(
    mask=('[#6 ]', ), )), u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
mdb.models['Model-1'].Coupling(alpha=0.0, controlPoint=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[635], )), couplingType=
    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name='right', 
    surface=Region(
    edges=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].edges.getSequenceFromMask(
    mask=('[#0 #48004000 #26009 #0 #80000 #0:2 #48024240', 
    ' #80 #1000 #0:8 #400200 #100 #0:6 #40 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].edges.getSequenceFromMask(
    mask=('[#8 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].edges.getSequenceFromMask(
    mask=('[#8 ]', ), ), 
    vertices=mdb.models['Model-1'].rootAssembly.instances['core_12mm-1'].vertices.getSequenceFromMask(
    mask=('[#80000104 #2c5da64 #0 #400 #4b6000 #2 #0', 
    ' #13000 #0 #80000000 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-1'].vertices.getSequenceFromMask(
    mask=('[#9 ]', ), )+\
    mdb.models['Model-1'].rootAssembly.instances['face_sheet-2'].vertices.getSequenceFromMask(
    mask=('[#9 ]', ), )), u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='left', 
    region=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[634], )), u1=SET, u2=SET
    , u3=SET, ur1=SET, ur2=UNSET, ur3=UNSET)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='right', 
    region=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[635], )), u1=UNSET, u2=
    SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
mdb.models['Model-1'].ConcentratedForce(cf1=-1.0, createStepName='Step-1', 
    distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
    Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[635], )))
mdb.models['Model-1'].parts['face_sheet'].setElementType(elemTypes=(ElemType(
    elemCode=S8R, elemLibrary=STANDARD), ElemType(elemCode=STRI65, 
    elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['face_sheet'].faces.getSequenceFromMask((
    '[#1 ]', ), ), ))
mdb.models['Model-1'].parts['face_sheet'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=1.0)
mdb.models['Model-1'].parts['face_sheet'].generateMesh()
mdb.models['Model-1'].parts['core_12mm'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=2.0)
mdb.models['Model-1'].parts['core_12mm'].generateMesh()
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models.changeKey(fromName='Model-1', toName='shell_model')
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='shell_model', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='12mm', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, numThreadsPerMpiProcess=1, queue=None, resultsFormat=
    ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, 
    waitMinutes=0)
# Save by anair on 2025_01_29-14.32.59; build 2023 2022_09_28-20.11.55 183150
