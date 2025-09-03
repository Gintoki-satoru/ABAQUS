from abaqus import *
from abaqusConstants import *
import sketch
import os
import numpy as np

import section
import regionToolset
import displayGroupMdbToolset as dgm
import mesh



path_modules = 'N:\\Sachdeva\\MT_Nair\\ABAQUS\\MT\\Macros'

os.chdir(path_modules)

# Further packages:
import coordinateTransformation_ellipse as ct

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

#############################   PARAMETERS    #############################

m_a_inner = 100.0  # Semi-major axis of the inner ellipse
m_b_inner = 100.0  # Semi-minor axis of the inner ellipse

thick = 2.5  # Thickness

m_a_outer = m_a_inner + thick  # Semi-major axis of the outer ellipse
m_b_outer = m_b_inner + thick  # Semi-minor axis of the outer ellipse

N = 4  # Number of partitions

mesh_size = 0.5  # Mesh size
edge_mesh = thick/N  # Edge mesh size

Press = 0.1 # Pressure load

##############################   MODEL CREATION   #############################

# Create a new model
Mdb()
modelName = 'EllipseModel_2D'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
mdb.models[modelName].rootAssembly.clearGeometryCache()
mdb.models[modelName].rootAssembly.regenerate()
myModel = mdb.models[modelName]

# Create a new sketch

s = mdb.models['EllipseModel_2D'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(viewStyle=AXISYM)
s.setPrimaryObject(option=STANDALONE)

s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
s.FixedConstraint(entity=g[2])

center = (0.0, 0.0)
majorAxisPoint = (m_a_inner, 0.0)
minorAxisPoint = (0.0, m_b_inner)

s.EllipseByCenterPerimeter(center=center, axisPoint1=majorAxisPoint, axisPoint2=minorAxisPoint)

outerMajorAxisPoint = (m_a_outer, 0.0)
outerMinorAxisPoint = (0.0, m_b_outer)

s.EllipseByCenterPerimeter(center=center, axisPoint1=outerMajorAxisPoint, axisPoint2=outerMinorAxisPoint)

s.autoTrimCurve(curve1=g.findAt((-m_a_inner, 0.0)), point1=(-m_a_inner, 0.0))
s.autoTrimCurve(curve1=g.findAt((-m_a_outer, 0.0)), point1=(-m_a_outer, 0.0))

phi = np.radians(300)
s.autoTrimCurve(
    curve1=g.findAt((ct.pol2cart_x(m_a_inner, phi), ct.pol2cart_y(m_b_inner, phi))),
    point1=(ct.pol2cart_x(m_a_inner, phi), ct.pol2cart_y(m_b_inner, phi))
)

s.autoTrimCurve(
    curve1=g.findAt((ct.pol2cart_x(m_a_outer, phi), ct.pol2cart_y(m_b_outer, phi))),
    point1=(ct.pol2cart_x(m_a_outer, phi), ct.pol2cart_y(m_b_outer, phi))
)
s.Line(point1=(0.0, m_b_inner), point2=(0.0, m_b_outer))

s.Line(point1=(m_a_inner, 0.0), point2=(m_a_outer, 0.0))

p = mdb.models['EllipseModel_2D'].Part(name='EllipsePart', dimensionality=AXISYMMETRIC, 
    type=DEFORMABLE_BODY)
p.BaseShell(sketch=s)

# Create Partition
p = mdb.models['EllipseModel_2D'].parts['EllipsePart']
f = p.faces

phi = np.radians(45)

for i in range(1, N):
    t = p.MakeSketchTransform(
        sketchPlane=f.findAt(coordinates=(ct.pol2cart_x(m_a_inner + i*thick/N, phi),
                                         ct.pol2cart_y(m_b_inner + i*thick/N, phi), 0.0)),
        sketchPlaneSide=SIDE1, origin=(0,0,0)
    )
    s = mdb.models['EllipseModel_2D'].ConstrainedSketch(name='__partition_%d__'%i, sheetSize=200, transform=t)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    
    s.EllipseByCenterPerimeter(center=(0, 0),
                               axisPoint1=(m_a_inner + i*thick/N, 0),
                               axisPoint2=(0, m_b_inner + i*thick/N))
    
    pickedFace = f.findAt(((ct.pol2cart_x(m_a_inner + i*thick/N, phi),
                            ct.pol2cart_y(m_b_inner + i*thick/N, phi), 0.0), ))
    
    p.PartitionFaceBySketch(faces=pickedFace, sketch=s)
    del mdb.models['EllipseModel_2D'].sketches['__partition_%d__'%i]

############ Material properties
myModel.Material(name='Aluminium')
myModel.materials['Aluminium'].Elastic(table=((70000.0,0.34), ))

############ Section
myModel.HomogeneousSolidSection(name='AL_section', 
    material='Aluminium', thickness=None)

def assign_sections():
    for i in range(1, N+1):
        faces = f.findAt(((ct.pol2cart_x(m_a_inner + i*thick/N, phi),
                           ct.pol2cart_y(m_b_inner + i*thick/N, phi), 0.0), ))
        region = regionToolset.Region(faces=faces)
        p.SectionAssignment(region=region, sectionName='AL_section',
                            offset=0.0, offsetType=MIDDLE_SURFACE,
                            offsetField='', thicknessAssignment=FROM_SECTION)

assign_sections()

############ Assembly
a = mdb.models['EllipseModel_2D'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='EllipsePart-1', part=p, dependent=OFF)
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])

# Crerate Sets
e = a.instances['EllipsePart-1'].edges
for i in range(0, N+1):
    edges = e.findAt(((ct.pol2cart_x(m_a_inner + i*thick/N, phi),
                       ct.pol2cart_y(m_b_inner + i*thick/N, phi), 0.0), ))
    a.Set(edges=edges, name='Curve-%d'%i)

f1 = a.instances['EllipsePart-1'].faces
for i in range(1, N+1):
    face1 = f1.findAt(((ct.pol2cart_x(m_a_inner + i*thick/N - thick/(2*N), phi), ct.pol2cart_y(m_b_inner + i*thick/N - thick/(2*N), phi), 0.0), ))
    a.Set(faces=face1, name='Face-%d'%i)

edges_combined = []
for i in range(1, N+1):
    edges = e.findAt(((m_a_inner + i*thick/N - thick/(2*N), 0 , 0.0), ))
    if edges is not None:
        a.Set(edges=[edges], name='side-%d' % i)
        edges_combined.append(edges)

if edges_combined:
    a.Set(edges=edges_combined, name='bottom_load')

edges_combined = []
for i in range(1, N+1):
    edge = e.findAt(((0, m_b_inner + i*thick/N - thick/(2*N), 0.0), ))
    if edge is not None:
        a.Set(edges=[edge], name='top-%d' % i)
        edges_combined.append(edge)

if edges_combined:
    a.Set(edges=edges_combined, name='top_load')

# Create Surface
side1Edges1 = e.findAt(((ct.pol2cart_x(m_a_inner, phi), ct.pol2cart_y(m_b_inner, phi), 0.0), ))
a.Surface(side1Edges=side1Edges1, name='load')


############ Step
myModel.StaticStep(name='LoadingStep', previous='Initial')
myModel.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'COORD'))

######## Boundary Conditions

region = a.sets['top_load']
mdb.models['EllipseModel_2D'].DisplacementBC(name='top_load_BC', 
    createStepName='LoadingStep', region=region, u1=0.0, u2=UNSET, ur3=UNSET, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

region = a.sets['bottom_load']
mdb.models['EllipseModel_2D'].YsymmBC(name='bottom_load_BC', 
    createStepName='LoadingStep', region=region, localCsys=None)

######## Load
region = a.surfaces['load']
mdb.models['EllipseModel_2D'].Pressure(name='PressureLoad', 
    createStepName='LoadingStep', region=region, distributionType=UNIFORM, 
    field='', magnitude=Press, amplitude=UNSET)

######## Mesh

for i in range(0, N+1):
    pickedEdges = a.sets['Curve-{}'.format(i)].edges
    a.seedEdgeBySize(edges=pickedEdges, size=mesh_size, deviationFactor=0.1, constraint=FINER)

pickedEdges = a.sets['top_load'].edges
a.seedEdgeBySize(edges=pickedEdges, size=edge_mesh, deviationFactor=0.1, constraint=FINER)

pickedEdges = a.sets['bottom_load'].edges
a.seedEdgeBySize(edges=pickedEdges, size=edge_mesh, deviationFactor=0.1, constraint=FINER)

elemType1 = mesh.ElemType(elemCode=CAX4R, elemLibrary=STANDARD, 
    secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
    distortionControl=DEFAULT)

for i in range(1, N+1):
    pickedRegions = (a.sets['Face-{}'.format(i)].faces, )
    a.setElementType(regions=pickedRegions, elemTypes=(elemType1,))

partInstances =(a.instances['EllipsePart-1'], )
a.generateMesh(regions=partInstances)

######## Job
job_dir = r'N:\Sachdeva\MT_Nair\FE'
if not os.path.exists(job_dir):
    os.makedirs(job_dir)

os.chdir(job_dir)

jobName = 'EllipseJob_2D'
myJob = mdb.Job(name=jobName, model='EllipseModel_2D', description='',
                type=ANALYSIS, memory=90, memoryUnits=PERCENTAGE,
                getMemoryFromAnalysis=True,
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
                echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)

myJob.submit(consistencyChecking=OFF)
myJob.waitForCompletion()