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


# path_modules = 'N:\\Sachdeva\\MT_Nair\\ABAQUS\\MT\\Macros'
path_modules = 'D:\\psingh\\MT\\ABAQUS\\MT\\Macros'
# path_modules = r"C:\Users\lenovo\Desktop\Aerospace\Thesis\ABAQUS\MT\Macros"
os.chdir(path_modules)

# Further packages:
import coordinateTransformation_ellipse as ct

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

#############################   PARAMETERS    #############################

# Note: All dimensions are to be given as a float (e.g., 10.0 not 10)

m_a_inner = 100.0  # Semi-major axis of the inner ellipse
m_b_inner = 100.0  # Semi-minor axis of the inner ellipse

thick = 5.0  # Thickness

m_a_outer = m_a_inner + thick
m_b_outer = m_b_inner + thick

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
mdb.models['EllipseModel_2D'].XsymmBC(name='top_load_BC', 
    createStepName='LoadingStep', region=region, localCsys=None)

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

########### Create set for ODB ###################

model_name = 'EllipseModel_2D'
instance_name = 'EllipsePart-1'
a_in = m_a_inner
b_in = m_b_inner
ds = mesh_size
set_name = 'NODE-0'

a_asm = mdb.models[model_name].rootAssembly
inst = a_asm.instances[instance_name]

# Bounding box size slightly smaller than mesh size
dx = ds * 0.8
dy = ds * 0.8

phi = 0
points = []
while phi <= math.pi/2:
    x = a_in * math.cos(phi)
    y = b_in * math.sin(phi)
    points.append((x, y))
    dphi = ds / math.sqrt((a_in*math.sin(phi))**2 + (b_in*math.cos(phi))**2)
    phi += dphi

all_nodes = []
for x, y in points:
    nodes_found = inst.nodes.getByBoundingBox(
        xMin=x-dx, xMax=x+dx,
        yMin=y-dy, yMax=y+dy,
        zMin=-1e-6, zMax=1e-6
    )
    all_nodes.extend(nodes_found)

unique_nodes = {n.label:n for n in all_nodes}
node_ids = [n.label for n in unique_nodes.values()]

a_asm.SetFromNodeLabels(
    name=set_name,
    nodeLabels=((instance_name, tuple(node_ids)),)
)

model_name = 'EllipseModel_2D'
instance_name = 'EllipsePart-1'
a_start = 150.0
b_start = 100.0      
ds = 0.5           

a_asm = mdb.models[model_name].rootAssembly
inst = a_asm.instances[instance_name]

# --- Bounding box size ---
dx = ds * 2
dy = ds * 2

for i in range(N):
    a_in = a_start + i * (thick / N)
    b_in = b_start + i * (thick / N)
    phi = 0.0
    points = []
    while phi <= math.pi/2:
        x = a_in * math.cos(phi)
        y = b_in * math.sin(phi)
        points.append((x, y))
        dphi = ds / math.sqrt((a_in*math.sin(phi))**2 + (b_in*math.cos(phi))**2)
        phi += dphi
    all_elements = []
    for x, y in points:
        elems_found = inst.elements.getByBoundingBox(
            xMin=x-dx/2, xMax=x+dx,
            yMin=y-dy/2, yMax=y+dy,
            zMin=-1e-6, zMax=1e-6 )
        all_elements.extend(elems_found)
    unique_elements = {e.label: e for e in all_elements}
    element_list = list(unique_elements.values())
    element_ids = [e.label for e in element_list]
    elem_set_name = 'ELEM-{}'.format(i+1)
    a_asm.SetFromElementLabels(
        name=elem_set_name,
        elementLabels=((instance_name, tuple(element_ids)),)
    )
    print('Created element set: {} with {} elements'.format(elem_set_name, len(element_ids)))

######## Job submission #################
job_dir = r'N:\Sachdeva\MT_Nair\FE'
# job_dir = r"C:\Users\lenovo\Desktop\Aerospace\Thesis\FE_temp"
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

##################### POST-PROCESSING ##########################

odb_path = 'EllipseJob_2D.odb'

if N % 2 == 0:
    middle_index = (N // 2) + 1
else:
    middle_index = (N // 2) + 1

elem_set_name = "ELEM-{}".format(middle_index)

step_name = 'LoadingStep'
output_csv = 'element_stress_ip_coords.csv'

odb = openOdb(odb_path)
assembly = odb.rootAssembly

try:
    elemset = assembly.elementSets[elem_set_name]
except KeyError:
    print('Element set "{}" not found in ODB.'.format(elem_set_name))
    odb.close()
    raise

# --- Last frame ---
frame = odb.steps[step_name].frames[-1]

stress_field = frame.fieldOutputs['S']
stress_subset = stress_field.getSubset(region=elemset, position=INTEGRATION_POINT)

coord_field = frame.fieldOutputs['COORD']
coord_subset = coord_field.getSubset(region=elemset, position=INTEGRATION_POINT)

data = []
for s_val, c_val in zip(stress_subset.values, coord_subset.values):
    data.append([
        s_val.elementLabel,
        s_val.integrationPoint,
        c_val.data[0],  # X (radial)
        c_val.data[1],  # Y (axial)
        s_val.data[0],  # S11
        s_val.data[1],  # S22
        s_val.data[2],  # S33
        s_val.data[3]   # S12
    ])

with open(output_csv, 'wb') as f:
    writer = csv.writer(f)
    writer.writerow(['Element', 'IntegrationPoint', 'X', 'Y', 'S11', 'S22', 'S33', 'S12'])
    for row in data:
        writer.writerow(row)

odb.close()
