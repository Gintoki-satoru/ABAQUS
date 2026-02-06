# -*- coding: mbcs -*-
# Do not delete the following import lines
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
import math
import numpy as np
## all dimensions are in mm
hcore= 2.0 ##10.0
# lf=50.0 # length of the beam
hf=0.1 #facesheets thickness 1.0
#choose the number of cells as a number divisible by core thickness
N_y=2
a_cell=float(hcore/N_y)
N_x= [46]
d_truss=a_cell/10.0 ##Earlier 0.05  
Num_eig_values=1
max_iterations=200
E_module=70000.0  ## MPa
poisson=0.35
N_beamElements=10.0
half_hcore=hcore/2.0
#Boundary Condistions
#case=1 (Default) Fixed from the left side
#case=2 Simply supported from the left side
#case=3 simply supported at the lower corners(Not Useful)
#case=4 Simply supported at the all corners
myBoundaryCondition=1
#Loads Ff=0,4984 Fc=0,00325   pressure=force/area
Es_Exx=1.0/0.0023618305245#1.0/(1.63*(a_cell/d_truss)**(-2.0)) 1.4441287085E-02
Ff=hf/hcore*Es_Exx/(1.0+2.0*hf/hcore*Es_Exx)
Fc=(1.0-2.0*Ff)
 #0.4984/(N_z*a_cell*N_y*a_cell) 
#Fc=#0.00325
fileNamePath='D:/Hatem/3dAR12_SR21/'
LoadFromLeft=0
LoadFromRight=1
MeshSizeFaceSheet=a_cell
MeshSizeCore=a_cell/1.0

truss_cross_area=math.pi*d_truss**2
# Num_cell_x=int(i_forNx)
# Num_cell_y=int(N_y)
# half_hcore=hcore/2.0
results = []
said=0
for i_forNx in N_x:    
    Mdb()
    N_z=i_forNx
    width=float(N_z*a_cell)
    Pf=Ff/(N_z*a_cell*hf)
    lf=float(i_forNx*a_cell)
    SR=float(lf/(hcore+2*hf))
    print('slenderness ratio = %f' % float(lf/(hcore+2*hf)))
    # data=(hcore, hf, lf, a_cell, d_truss,E_module,poisson, Num_eig_values, max_iterations,myBoundaryCondition)
    #results.append(data)
#csv_file_path = 'results.csv'
#combined_data = np.vstack(results)
#np.savetxt(csv_file_path, combined_data, delimiter=',', header='Buckloads,BeamLength,CoreThickness,faceSheetThickness', comments='')
# print('Done!')
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.Line(point1=(0.0, 0.0), point2=(0.0, a_cell))
    s.VerticalConstraint(entity=g[2], addUndoState=False)
    s.Line(point1=(0.0, a_cell), point2=(a_cell/2.0, a_cell/2.0))
    s.Line(point1=(a_cell/2.0, a_cell/2), point2=(0.0, 0.0))
    s.Line(point1=(a_cell/2.0, a_cell/2.0), point2=(a_cell, 0.0))
    s.Line(point1=(a_cell/2.0, a_cell/2.0), point2=(a_cell, a_cell))
    s.Line(point1=(a_cell, 0.0), point2=(a_cell, a_cell))
    s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    # s.linearPattern(geomList=(g[2], g[3], g[4]), vertexList=(), number1=1, 
        # spacing1=20.0, angle1=0.0, number2=N_y, spacing2=a_cell, angle2=90.0)
    p = mdb.models['Model-1'].Part(name='core', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['core']
    p.BaseWire(sketch=s)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['core']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    # creat points and lines for one cell
    p.DatumPointByCoordinate(coords=(0.0, a_cell/2.0, -a_cell/2.0))
    p.DatumPointByCoordinate(coords=(0.0, a_cell, -a_cell))
    p.DatumPointByCoordinate(coords=(0.0, 0.0, -a_cell))
    p.DatumPointByCoordinate(coords=(a_cell/2.0, a_cell/2.0, -a_cell))
    p.DatumPointByCoordinate(coords=(a_cell, a_cell/2.0, -a_cell/2.0))
    p.DatumPointByCoordinate(coords=(a_cell, a_cell, -a_cell))
    p.DatumPointByCoordinate(coords=(a_cell, 0.0, -a_cell))
    ##Lines
    p.WirePolyLine(points=(([0.0, 0.0,0.0 ], [0.0, a_cell/2.0,-a_cell/2.0 ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([0.0, a_cell,0.0 ], [0.0, a_cell/2.0,-a_cell/2.0 ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([0.0, a_cell/2.0,-a_cell/2.0 ], [0.0, a_cell,-a_cell ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([0.0, a_cell/2.0,-a_cell/2.0 ], [0.0, 0.0,-a_cell ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([0.0, 0.0,-a_cell  ], [a_cell/2.0, a_cell/2.0, -a_cell ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([0.0, a_cell,-a_cell  ], [a_cell/2.0, a_cell/2.0, -a_cell ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([0.0, 0.0,-a_cell  ], [0.0, a_cell,-a_cell  ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([a_cell/2.0, a_cell/2.0, -a_cell ], [a_cell, a_cell, -a_cell ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([a_cell/2.0, a_cell/2.0, -a_cell], [a_cell, 0.0, -a_cell]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([a_cell, 0.0,-a_cell  ], [a_cell, a_cell,-a_cell  ]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([a_cell, 0.0, -a_cell], [a_cell, a_cell/2.0, -a_cell/2.0]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([a_cell, a_cell, -a_cell], [a_cell, a_cell/2.0, -a_cell/2.0]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([a_cell, a_cell/2.0, -a_cell/2.0], [a_cell, 0.0, 0.0]), ), mergeType=IMPRINT, meshable=ON)
    p.WirePolyLine(points=(([a_cell, a_cell/2.0, -a_cell/2.0], [a_cell, a_cell, 0.0]), ), mergeType=IMPRINT, meshable=ON)
    #Facesheets
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.rectangle(point1=(-lf/2.0, -hf/2.0), point2=(lf/2.0, hf/2.0))
    p = mdb.models['Model-1'].Part(name='Facesheet', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Facesheet']
    p.BaseSolidExtrude(sketch=s1, depth=width)
    s1.unsetPrimaryObject()
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    #Creat Material
    mdb.models['Model-1'].Material(name='AL')
    mdb.models['Model-1'].materials['AL'].Elastic(table=((E_module, poisson), ))
    #Sections
    mdb.models['Model-1'].HomogeneousSolidSection(name='faceSheetSection', 
        material='AL', thickness=None)
    mdb.models['Model-1'].CircularProfile(name='Profile-1', r=d_truss/2.0)
    mdb.models['Model-1'].BeamSection(name='coreSection', 
        integration=DURING_ANALYSIS, poissonRatio=0.0, profile='Profile-1', 
        material='AL', temperatureVar=LINEAR, consistentMassMatrix=False)
    #Assembly
    #Core
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['core']
    a.Instance(name='core-1', part=p, dependent=ON)
    a.LinearInstancePattern(instanceList=('core-1', ), direction1=(1.0, 0.0, 0.0), 
        direction2=(0.0, 1.0, 0.0), number1=i_forNx, number2=N_y, spacing1=a_cell, 
        spacing2=a_cell)
    a.InstanceFromBooleanMerge(name='CoreXY', instances=a.instances.values(),originalInstances=DELETE, domain=BOTH)
    a.LinearInstancePattern(instanceList=('CoreXY-1', ), direction1=(1.0, 0.0, 0.0), 
        direction2=(0.0, 0.0, -1.0), number1=1, number2=N_z, spacing1=a_cell, 
        spacing2=a_cell)
    a.InstanceFromBooleanMerge(name='Core', instances=a.instances.values(),originalInstances=DELETE, domain=BOTH)
    a.translate(instanceList=('Core-1', ), vector=(-lf/2.0, -hcore/2.0, width/2.0))
    #facesheets
    p = mdb.models['Model-1'].parts['Facesheet']
    a.Instance(name='UpperFacesheet', part=p, dependent=ON)
    a.translate(instanceList=('UpperFacesheet', ), vector=(0.0, half_hcore+hf/2.0, -width/2.0))
    a.Instance(name='BottomFacesheet', part=p, dependent=ON)
    a.translate(instanceList=('BottomFacesheet', ), vector=(0.0, -half_hcore-hf/2.0, -width/2.0))
    ########################################
    a = mdb.models['Model-1'].rootAssembly
    v = a.instances['Core-1'].vertices
    verts = v.findAt(((0.0,0.0,0.0), ))
    a.Set(vertices=verts, name='CoreCenterPoint')
    verts = v.findAt(((lf/2.0,0.0,0.0), ))
    a.Set(vertices=verts, name='CoreCenterPointRight')
    verts = v.findAt(((-lf/2.0,0.0,0.0), ))
    a.Set(vertices=verts, name='CoreCenterPointLeft')
    #Assign sections
    p = mdb.models['Model-1'].parts['Core']
    e = p.edges
    edges = e.getByBoundingCylinder(center1=(-lf, 0.0,0.0), center2=(lf, 0.0,0.0), radius=width+hcore+lf)
    region = regionToolset.Region(edges=edges)
    p = mdb.models['Model-1'].parts['Core']
    p.SectionAssignment(region=region, sectionName='coreSection', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    p1 = mdb.models['Model-1'].parts['Facesheet']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Model-1'].parts['Facesheet']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(cells=cells)
    p = mdb.models['Model-1'].parts['Facesheet']
    p.SectionAssignment(region=region, sectionName='faceSheetSection', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    #Creat Setes
    #Core
    p = mdb.models['Model-1'].parts['Core']
    v = p.vertices
    verts = v.getByBoundingBox(0.0, hcore-a_cell/4.0, -width-a_cell/4.0, lf, hcore+a_cell/4.0,0.0)
    p.Set(vertices=verts, name='UpperPoints')
    verts = v.getByBoundingBox(0.0, -a_cell/4.0, -width-a_cell/4.0, lf, a_cell/4.0,0.0)
    p.Set(vertices=verts, name='BottomPoints')
    verts = v.getByBoundingBox(0.0, half_hcore-a_cell/4.0, -a_cell/4.0, lf, half_hcore+a_cell/4.0,a_cell/4.0)
    p.Set(vertices=verts, name='CoreFaceFrontPoints')
    verts = v.getByBoundingBox(0.0, half_hcore-a_cell/4.0, -width-a_cell/4.0, lf, half_hcore+a_cell/4.0,-width+a_cell/4.0)
    p.Set(vertices=verts, name='CoreFaceBackPoints')
    verts = v.getByBoundingBox(lf-a_cell/4.0, half_hcore-a_cell/4.0, -width-a_cell/4.0, lf+a_cell/4.0, half_hcore+a_cell/4.0,a_cell/4.0)
    p.Set(vertices=verts, name='CoreFaceRightPoints')
    verts = v.getByBoundingBox(-a_cell/4.0, half_hcore-a_cell/4.0, -width-a_cell/4.0, a_cell/4.0, half_hcore+a_cell/4.0,a_cell/4.0)
    p.Set(vertices=verts, name='CoreFaceLeftPoints')
    verts = v.getByBoundingBox(-a_cell/4.0, -a_cell/4.0, -width-a_cell/4.0, a_cell/4.0, hcore+a_cell/4.0,a_cell/4.0)
    p.Set(vertices=verts, name='CoreFaceALLLeftPoints')
    verts = v.getByBoundingBox(lf/2.0-a_cell/4.0, half_hcore-a_cell/4.0, -width-a_cell/4.0, lf/2.0+a_cell/4.0, half_hcore+a_cell/4.0,a_cell/4.0)
    p.Set(vertices=verts, name='CoreFaceCenterPoints')
    verts = v.getByBoundingBox(-lf-a_cell/4.0, half_hcore-a_cell/4.0, -width/2.0-a_cell/4.0, lf+a_cell/4.0, half_hcore+a_cell/4.0,-width/2.0+a_cell/4.0)
    p.Set(vertices=verts, name='CoreFaceCenterPoints2')
    verts = v.getByBoundingBox(lf-a_cell/4.0, -a_cell/4.0, -width-a_cell/4.0, lf+a_cell/4.0, hcore+a_cell/4.0,a_cell/4.0)
    p.Set(vertices=verts, name='CoreFaceALLRightPoints')
    #faceSheetSection
    # s = p.faces
    # side1Faces = s.getSequenceFromMask(mask=('[#2 ]', ), )
    # p.Surface(side1Faces=side1Faces, name='upperface')
    p = mdb.models['Model-1'].parts['Facesheet']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    s = p.faces
    side1Faces=s.findAt(((0.0, -hf/2.0,width/2.0),))
    p.Surface(side1Faces=side1Faces, name='lowerface')
    side1Faces=s.findAt(((0.0, hf/2.0,width/2.0),))
    p.Surface(side1Faces=side1Faces, name='upperface')
    side1Faces=s.findAt(((lf/2.0,0.0,width/2.0 ),))
    p.Surface(side1Faces=side1Faces, name='RightFace')
    side1Faces=s.findAt(((-lf/2.0,0.0,width/2.0 ),))
    p.Surface(side1Faces=side1Faces, name='LeftFace')
    #Constraints
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['UpperFacesheet'].surfaces['lowerface']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.instances['Core-1'].sets['UpperPoints']
    mdb.models['Model-1'].Tie(name='UpperConstraint', master=region1, 
        slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
        tieRotations=ON, constraintEnforcement=NODE_TO_SURFACE, thickness=ON)
    region1=a.instances['BottomFacesheet'].surfaces['upperface']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.instances['Core-1'].sets['BottomPoints']
    mdb.models['Model-1'].Tie(name='BottomConstraint', master=region1, 
        slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
        tieRotations=ON, constraintEnforcement=NODE_TO_SURFACE, thickness=ON)
    #BucklingStep
    mdb.models['Model-1'].BuckleStep(name='Buckling', previous='Initial', 
        numEigen=1, eigensolver=LANCZOS, minEigen=0.0, blockSize=DEFAULT, 
        maxBlocks=DEFAULT)
    # mdb.models['Model-1'].BuckleStep(name='Buckling', previous='Initial', 
        # numEigen=Num_eig_values, vectors=Num_eig_values+8, maxIterations=max_iterations)
    #facesheets sets
    p = mdb.models['Model-1'].parts['Facesheet']
    e = p.edges
    edges = e.findAt(((0.0,hf/2.0,0.0), ),
                 ((lf/2.0,hf/2.0,width/2.0), ),
                 ((-lf/2.0,hf/2.0,width/2.0), ),
                 ((0.0,hf/2.0,width),))
    p.Set(edges=edges, name='upperedges')
    edges = e.findAt(((0.0,-hf/2.0,0.0), ),
                 ((lf/2.0,-hf/2.0,width/2.0), ),
                 ((-lf/2.0,-hf/2.0,width/2.0), ),
                 ((0.0,-hf/2.0,width),))
    p.Set(edges=edges, name='loweredges')
    # facesheet partition
    p = mdb.models['Model-1'].parts['Facesheet']
    plane1=p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=lf/2.0)
    plane2=p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    d1 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d1[plane1.id], cells=pickedCells)
    pickedCells = c.getSequenceFromMask(mask=('[#3 ]', ), )
    p.PartitionCellByDatumPlane(datumPlane=d1[plane2.id], cells=pickedCells)
    #facesheets sets points
    p = mdb.models['Model-1'].parts['Facesheet']
    v = p.vertices
    verts = v.findAt(((0.0,hf/2.0,0.0), ),
                 ((0.0,hf/2.0,width),))
    p.Set(vertices=verts, name='upperpoints_U1')
    verts = v.findAt(((lf/2.0,hf/2.0,width/2.0), ),
                 ((-lf/2.0,hf/2.0,width/2.0),))
    p.Set(vertices=verts, name='upperpoints_U3')
    verts = v.findAt(((0.0,-hf/2.0,0.0), ),
                 ((0.0,-hf/2.0,width),))
    p.Set(vertices=verts, name='lowerpoints_U1')
    verts = v.findAt(((lf/2.0,-hf/2.0,width/2.0), ),
                 ((-lf/2.0,-hf/2.0,width/2.0),))
    p.Set(vertices=verts, name='lowerpoints_U3')
    #BoundaryConditions  
    a = mdb.models['Model-1'].rootAssembly
    region = a.instances['UpperFacesheet'].sets['upperedges']
    mdb.models['Model-1'].DisplacementBC(name='BCedgesUp', createStepName='Initial', 
        region=region, u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    region = a.instances['BottomFacesheet'].sets['loweredges']
    mdb.models['Model-1'].DisplacementBC(name='BCedgesBottom', createStepName='Initial', 
        region=region, u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    region = a.instances['UpperFacesheet'].sets['upperpoints_U1']
    mdb.models['Model-1'].DisplacementBC(name='BCpointUPu1', createStepName='Initial', 
        region=region, u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    region = a.instances['UpperFacesheet'].sets['upperpoints_U3']
    mdb.models['Model-1'].DisplacementBC(name='BCpointUPu3', createStepName='Initial', 
        region=region, u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    region = a.instances['BottomFacesheet'].sets['lowerpoints_U1']
    mdb.models['Model-1'].DisplacementBC(name='BCpointLowU1', createStepName='Initial', 
        region=region, u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    region = a.instances['BottomFacesheet'].sets['lowerpoints_U3']
    mdb.models['Model-1'].DisplacementBC(name='BCpointLowU3', createStepName='Initial', 
        region=region, u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    #Loads
    region = a.instances['Core-1'].sets['CoreFaceALLLeftPoints']
    numbernodeload=len(a.instances['Core-1'].sets['CoreFaceALLLeftPoints'].vertices)
    mdb.models['Model-1'].ConcentratedForce(name='CoreLeftLoad', createStepName='Buckling', 
        region=region, cf1=Fc/numbernodeload, distributionType=UNIFORM, field='', 
        localCsys=None)
    region = a.instances['Core-1'].sets['CoreFaceALLRightPoints']
    mdb.models['Model-1'].ConcentratedForce(name='CoreRightLoad', createStepName='Buckling', 
        region=region, cf1=-Fc/numbernodeload, distributionType=UNIFORM, field='', 
        localCsys=None)
    ######################
    region = a.instances['UpperFacesheet'].surfaces['RightFace']
    mdb.models['Model-1'].Pressure(name='UpperRightFaceLoad', 
        createStepName='Buckling', region=region, distributionType=UNIFORM, 
        field='', magnitude=Pf)
    region = a.instances['UpperFacesheet'].surfaces['LeftFace']
    mdb.models['Model-1'].Pressure(name='UpperLeftFaceLoad', 
        createStepName='Buckling', region=region, distributionType=UNIFORM, 
        field='', magnitude=Pf)
    ################################    
    region = a.instances['BottomFacesheet'].surfaces['LeftFace']
    mdb.models['Model-1'].Pressure(name='BottomLeftFaceLoad', 
        createStepName='Buckling', region=region, distributionType=UNIFORM, 
        field='', magnitude=Pf)  
    region = a.instances['BottomFacesheet'].surfaces['RightFace']
    mdb.models['Model-1'].Pressure(name='BottomRightFaceLoad', 
        createStepName='Buckling', region=region, distributionType=UNIFORM, 
        field='', magnitude=Pf)
    #################################
    #Mesh
    p = mdb.models['Model-1'].parts['Facesheet']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    pickedRegions =(cells, )
    p.seedPart(size=MeshSizeFaceSheet, deviationFactor=0.1, minSizeFactor=0.1)
    elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
    p.generateMesh()  
    #Core
    core = mdb.models['Model-1'].parts['Core'] 
    edges_temp = core.edges                # Create temporary set with all core edges
    region_temp = regionToolset.Region(edges = edges_temp)    
    core.assignBeamSectionOrientation(region=region_temp, method=N1_COSINES, n1=(0.0, 0.0,-1.0)) 
    # Assign beam section orientation. The exact orientation is computed
    # automatically by abaqus somewhere perpendicular to the beam's direction
    core.seedPart(size =MeshSizeCore , deviationFactor = 0.1, minSizeFactor = 0.1)
    core_elemType = mesh.ElemType(elemCode = B31, elemLibrary = STANDARD)
    # Choose element. Element is 3 dimensional linear truss element
    region_temp = regionToolset.Region(edges = core.edges) # Create region with core edges
    core.setElementType(regions = region_temp, elemTypes = (core_elemType, ))
    core.generateMesh()
    a1 = mdb.models['Model-1'].rootAssembly
    a1.regenerate()    
    # elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
    # elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
    # elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD, 
        # secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    #session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
    #Creat Job
    # mdb.Job(name='BucklingJob', model='Model-1', description='', type=ANALYSIS, 
            # atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            # memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            # explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            # modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            # scratch='D:\\Abaqus_scratch', resultsFormat=ODB, 
            # multiprocessingMode=DEFAULT, numCpus=2, numDomains=2, numGPUs=1)
    # mdb.jobs['BucklingJob'].submit(consistencyChecking=OFF)
    # mdb.jobs['BucklingJob'].waitForCompletion()
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=('U',))
    #################
    ################
    a=('{:.0f}'.format(i_forNx*a_cell))
    b=('{:.0f}'.format(N_z*a_cell))
    c=('{:.0f}'.format(hcore+2*hf))
    c11=('{:.0f}'.format(round(d_truss,1)))
    fileNameApp=('Job_L'+str(a)+'_Hc'+str(b)+'_Hs'+str(c)+'_t'+str(c11))
    mdb.saveAs('mdb_'+'Job_L'+str(a)+'_Hc'+str(b)+'_Hs'+str(c)+'_t'+str(c11))
    jobName='Job_L'+str(a)+'_Hc'+str(b)+'_Hs'+str(c)+'_t'+str(c11)
    #################################################################################################
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF,
        multiprocessingMode=DEFAULT, name=jobName, nodalOutputPrecision=SINGLE,
        numCpus=3, numDomains=3, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
        'D:\\Abaqus_scratch', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    #Submit Job: #################################################################
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()
    #########################################################################################################
    #########################################################################################################
    #########################################################################################################
    ###########################      Postprocessing     #####################################################
    #########################################################################################################
    #########################################################################################################
    #########################################################################################################
    myOdb = session.openOdb(name=fileNamePath+jobName+'.odb')
    myFrame = myOdb.steps['Buckling'].frames
    Frame_length=len(myFrame)
    ##########################################################################################################
    ##########################################################################################################
    BuckleLoad=[]
    for i in range(1,Frame_length):
        EigenV=myFrame[i].description
        list=EigenV.split()
        Load=float(list[-1])
        BuckleLoad=BuckleLoad+[Load]
    print(BuckleLoad)
    ##########################################################################################################
    ##########################################################################################################
    if said==0:
        outFile1=open('Eigenvalues.txt','w')
    for count, value in enumerate(BuckleLoad):
        valueinmm=value/lf
        outFile1.write('%10.10E %10.10E %10.10E %10.10E \n'  % (count+1, SR, value, valueinmm))
    if said==len(N_x)-1:
        outFile1.close()
    said=said+1
    #################################################
    # myOdb = session.openOdb(name=fileNamePath+jobName+'.odb')
    # myFrame = myOdb.steps['BucklingJob'].frames
    # Frame_length=len(myFrame)
    # ##########################################################################################################
    # ##########################################################################################################
    # BuckleLoad=[]
    # for i in range(1,Frame_length):
        # EigenV=myFrame[i].description
        # list=EigenV.split()
        # Load=float(list[-1])
        # BuckleLoad=BuckleLoad+[Load]
    # print(BuckleLoad)
    # ##########################################################################################################
    # ##########################################################################################################
    # outFile1=open('Eigenvalues.txt','w')
    # for count, value in enumerate(BuckleLoad):
        # outFile1.write('%10.10E %10.10E \n'  % (count+1, value))
    # outFile1.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    