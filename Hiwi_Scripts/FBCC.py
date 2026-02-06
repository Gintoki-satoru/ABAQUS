from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

# ----------------------------------------------------------------
# 1. Create a new part for the FBCC unit cell
# ----------------------------------------------------------------
modelName = 'Model-1'
cellName  = 'FBCC_Cell'

model = mdb.models[modelName]
p = model.Part(name=cellName, dimensionality=THREE_D, type=DEFORMABLE_BODY)

# ----------------------------------------------------------------
# 2. Define cube size
# ----------------------------------------------------------------
L = 4.0  # mm, change as needed

corner_coords = [
    (0.0, 0.0, 0.0),
    (L,   0.0, 0.0),  
    (0.0, L,   0.0),  
    (L,   L,   0.0),  
    (0.0, 0.0, L),    
    (L,   0.0, L),    
    (0.0, L,   L),    
    (L,   L,   L)    
]

# ----------------------------------------------------------------
# 3. Create datum points and store IDs
# ----------------------------------------------------------------
dids = {} 
for i, xyz in enumerate(corner_coords):
    dp = p.DatumPointByCoordinate(coords=xyz)
    dids[i] = dp.id

d = p.datums

def add_strut(i, j):
    p.WirePolyLine(points=((d[dids[i]], d[dids[j]]),),
                   mergeType=IMPRINT, meshable=ON)

left_X = [(0,6), (2,4)]
right_X = [(1,7), (3,5)]
bottom_X = [(0,3), (1,2)]
top_X = [(4,7), (5,6)]

x_face_struts = left_X + right_X + bottom_X + top_X

diag_struts = [(0,7), (1,6), (2,5), (3,4)]

for i, j in x_face_struts + diag_struts:
    add_strut(i, j)


mdb.models['Model-1'].Material(name='AL')
mdb.models['Model-1'].materials['AL'].Elastic(table=((70000.0, 0.35), ))

mdb.models['Model-1'].Material(name='carbon')
mdb.models['Model-1'].materials['carbon'].Elastic(type=ENGINEERING_CONSTANTS, 
    table=((181000.0, 10300.0, 10300.0, 0.277, 0.277, 0.0, 7170.0, 7170.0, 
    5960.0), ))

mdb.models['Model-1'].CircularProfile(name='Profile-1', r=0.2)
mdb.models['Model-1'].BeamSection(name='strut', integration=DURING_ANALYSIS, 
    poissonRatio=0.0, profile='Profile-1', material='AL', 
    temperatureVar=LINEAR, consistentMassMatrix=False)
p = mdb.models['Model-1'].parts['FBCC_Cell']

edges_all = p.edges.getByBoundingBox(
    xMin=-1e9, xMax=1e9,
    yMin=-1e9, yMax=1e9,
    zMin=-1e9, zMax=1e9
)
p.Set(edges=edges_all, name='FBCC_AllStruts')
p.SectionAssignment(
    region=p.sets['FBCC_AllStruts'],
    sectionName='strut',
    offset=0.0,
    offsetType=MIDDLE_SURFACE,
    offsetField='',
    thicknessAssignment=FROM_SECTION
)
