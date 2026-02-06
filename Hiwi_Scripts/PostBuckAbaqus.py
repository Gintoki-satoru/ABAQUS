# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

def runSimulation(a,b,A11,A12,A22,A66,D11,D12,D16,D22,D26,D66,meshSize,lambdaEnd,evalAtQuater,K11,K22):
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

    nodesX = a/meshSize+1
    nodesY = b/meshSize+1

    ## delete old and open new
    Mdb()
    session.viewports['Viewport: 1'].setValues(displayedObject=None)

    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(a, b))

    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    
    #ABD
    mdb.models['Model-1'].GeneralStiffnessSection(name='ABD', referenceTemperature=None, 
                                                  stiffnessMatrix=(A11, A12, A22, 0.0, 0.0, A66, 
                                                                   0.0, 0.0, 0.0, 
                                                                   D11, 0.0, 0.0, 0.0, 
                                                                   D12, D22, 0.0, 0.0, 0.0, 
                                                                   D16, D26, D66), 
                                                                   applyThermalStress=0, poissonDefinition=DEFAULT, useDensity=OFF)
    mdb.models['Model-1'].sections['ABD'].setValues()
    mdb.models['Model-1'].sections['ABD'].TransverseShearShell(k11=K11, 
    k22=K22, k12=0.0)
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = p.Set(faces=faces, name='Plate')
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName='ABD', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
        engineeringFeatures=OFF, mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    
    ## mesh part
    p.seedPart(size=meshSize, minSizeFactor=0.1)
    elemType1 = mesh.ElemType(elemCode=S8R, elemLibrary=STANDARD)
    faces = p.faces
    p.setElementType(regions=(faces,), elemTypes=(elemType1,))
    p.generateMesh()

    ## create instance
    mbdAssembly = mdb.models['Model-1'].rootAssembly
    mbdAssembly.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['Part-1']
    mbdAssembly.Instance(name='Part-1-1', part=p, dependent=ON)

    ## step solve eigen define
    mdb.models['Model-1'].BuckleStep(name='Step-1', previous='Initial', numEigen=5, 
        vectors=30, maxIterations=1000)

    # ref points
    mbdAssembly = mdb.models['Model-1'].rootAssembly
    mbdAssembly.ReferencePoint(point=(-0.5, b/2, 0.0))        #left
    mbdAssembly.ReferencePoint(point=(a/2, -0.5, 0.0))      #bottom
    mbdAssembly.ReferencePoint(point=(a+0.5, b/2, 0.0))       #right
    mbdAssembly.ReferencePoint(point=(a/2, b+0.5, 0.0))       #top

    ## edgeSets
    e1 = mbdAssembly.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#8 ]', ), )
    mbdAssembly.Set(edges=edges1, name='edgeLeft')
    
    e1 = mbdAssembly.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#4 ]', ), )
    mbdAssembly.Set(edges=edges1, name='edgeTop')

    e1 = mbdAssembly.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    mbdAssembly.Set(edges=edges1, name='edgeRight')

    e1 = mbdAssembly.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#1 ]', ), )
    mbdAssembly.Set(edges=edges1, name='edgeBottom')

    ## refPoints Sets
    r1 = mbdAssembly.referencePoints

    refPoints1=(r1[7], )
    mbdAssembly.Set(referencePoints=refPoints1, name='refPointTop')

    refPoints1=(r1[6], )
    mbdAssembly.Set(referencePoints=refPoints1, name='refPointRight')

    refPoints1=(r1[5], )
    mbdAssembly.Set(referencePoints=refPoints1, name='refPointBottom')

    refPoints1=(r1[4], )
    mbdAssembly.Set(referencePoints=refPoints1, name='refPointLeft')

    ## constaints perpendicular straight edges
    mdb.models['Model-1'].Equation(name='Constraint-Top', terms=((-1.0, 'edgeTop', 
        2), (1.0, 'refPointTop', 2)))
    mdb.models['Model-1'].Equation(name='Constraint-Right', terms=((-1.0, 
        'edgeRight', 1), (1.0, 'refPointRight', 1)))
    mdb.models['Model-1'].Equation(name='Constraint-Bottom', terms=((-1.0, 
        'edgeBottom', 2), (1.0, 'refPointBottom', 2)))
    mdb.models['Model-1'].Equation(name='Constraint-Left', terms=((-1.0, 
        'edgeLeft', 1), (1.0, 'refPointLeft', 1)))

    ## Load
    s1 = mbdAssembly.instances['Part-1-1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
    region = mbdAssembly.Surface(side1Edges=side1Edges1, name='Surf-2')
    mdb.models['Model-1'].ShellEdgeLoad(name='Load-2', createStepName='Step-1', 
        region=region, magnitude=1.0, distributionType=UNIFORM, field='', 
        localCsys=None)
    
    ## BC U3=0
    e1 = mbdAssembly.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#f ]', ), )
    r2 = mbdAssembly.referencePoints
    refPoints2=(r2[4], r2[5], r2[6], r2[7], )
    region = mbdAssembly.Set(edges=edges1, referencePoints=refPoints2, name='Set-U30')
    mdb.models['Model-1'].DisplacementBC(name='BC-U30', createStepName='Step-1', 
        region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, 
        fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    ## BC Left
    r1 =  mbdAssembly.referencePoints
    refPoints1=(r1[4], )
    region =  mbdAssembly.Set(referencePoints=refPoints1, name='Set-10')
    mdb.models['Model-1'].DisplacementBC(name='BC-Left', createStepName='Step-1', 
        region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, 
        fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    ## BC Bottom
    r1 =  mbdAssembly.referencePoints
    refPoints1=(r1[5], )
    region =  mbdAssembly.Set(referencePoints=refPoints1, name='Set-11')
    mdb.models['Model-1'].DisplacementBC(name='BC-Bootm', createStepName='Step-1', 
        region=region, u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, 
        fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    ## create Job
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
        predefinedFields=OFF, connectors=OFF)
    mdb.Job(name='LinBuck', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
    
    
    ## copy LinModel to create postbuck Model
    p1 = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    mdb.Model(name='Model-2', objectToCopy=mdb.models['Model-1'])
    p = mdb.models['Model-2'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)

    import job

    ## ouputfile U of linBuck
    mdb.models['Model-1'].keywordBlock.synchVersions(storeNodesAndElements=False)
    mdb.models['Model-1'].keywordBlock.replace(65, """
    *Output, field, variable=PRESELECT
    *NODE FILE
    U""")

    ## submit
    mdb.jobs['LinBuck'].submit(consistencyChecking=OFF)

    print('waiting')
    mdb.jobs['LinBuck'].waitForCompletion()
    print('not waiting')

    ## extract eigen Values
    r = session.openOdb(name='LinBuck.odb')
    Frames = r.steps['Step-1'].frames
    f = len(Frames)
    F_crit = []
    for i in range(1,f):
        eigV = Frames[i].description
        list = eigV.split()
        Load = float(list[-1])
        F_crit = F_crit +[Load]
    print(F_crit)
    F_crit = F_crit[0]
    print('real F_crit is:')
    print(F_crit)
    ## export Fcrit
    text_file = open("Fcrit.txt", "w")
    text_file.write(str(F_crit))
    text_file.close()

    ###################################################
    ##Post Buck
    ###################################################

    ## replace eigensolver with nonlin riks solver
    mbdAssembly = mdb.models['Model-2'].rootAssembly
    mdb.models['Model-2'].StaticRiksStep(name='Step-1', previous='Initial', 
        maintainAttributes=True, maxNumInc=100, initialArcInc=0.001, 
        minArcInc=1E-36, maxArcInc=1, totalArcLength=1, maxLPF=lambdaEnd, nlgeom=ON) #for modechange: totalArcLength=0.7, else: totalArcLength=1
    
    ## output request
    mdb.models['Model-2'].fieldOutputRequests['F-Output-1'].setValues(variables=(
        'S', 'E', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'SF', 
        'CSTRESS', 'CDISP'))
    
    ## imperfection
    mdb.models['Model-2'].keywordBlock.synchVersions(storeNodesAndElements=False)
    mdb.models['Model-2'].keywordBlock.replace(48, """
    ** ----------------------------------------------------------------
    *IMPERFECTION, File=LinBuck, Step=1
    1, 0.6
    2, 0
    3, 0
    4, 0
    5, 0
    6, 0
    ** ----------------------------------------------------------------
    **
    ** STEP: Step-1
    **""")

    ## create nodeset of center node
    mbdAssembly = mdb.models['Model-2'].rootAssembly
    allNodes = mbdAssembly.instances['Part-1-1'].nodes
    selection = allNodes.getByBoundingBox(a/2-0.000001, b/2-0.000001, -0.000001, a/2+0.000001, b/2+0.000001, 0.000001)
    mbdAssembly.Set(name="midNode", nodes=selection)

    if evalAtQuater:
        ## create nodeset of left quater-center node
        mbdAssembly = mdb.models['Model-2'].rootAssembly
        allNodes = mbdAssembly.instances['Part-1-1'].nodes
        selection = allNodes.getByBoundingBox(a/4-0.000001, b/2-0.000001, -0.000001, a/4+0.000001, b/2+0.000001, 0.000001)
        mbdAssembly.Set(name="quaterMidNode", nodes=selection)

    ## create nodeset of right edge center node
    mbdAssembly = mdb.models['Model-2'].rootAssembly
    allNodes = mbdAssembly.instances['Part-1-1'].nodes
    selection = allNodes.getByBoundingBox(a-0.000001, b/2-0.000001, -0.000001, a+0.000001, b/2+0.000001, 0.000001)
    mbdAssembly.Set(name="midNodeRight", nodes=selection)

    regionDef=mdb.models['Model-2'].rootAssembly.sets['midNode']
    mdb.models['Model-2'].HistoryOutputRequest(name='H-Output-2', 
    createStepName='Step-1', variables=('U3', ), region=regionDef, 
    sectionPoints=DEFAULT, rebar=EXCLUDE)

    ## scale Load
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
    mdb.models['Model-2'].loads['Load-2'].setValues(magnitude=F_crit, resultant=ON)

    ## create Job
    mdb.Job(name='PostBuck', model='Model-2', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
    
    ## submit job
    mdb.jobs['PostBuck'].submit(consistencyChecking=OFF)

    print('waiting')
    mdb.jobs['PostBuck'].waitForCompletion()
    print('Post buck done')

a0 = 200
b0 = 200
A110 = 1.1328e+05
A220 = 1.1328e+05
A120 = 8.4595e+04
A660 = 9.3204e+04
A160 = 0.0
A260 = 0.0

D110 = 7.0797e+05
D220 = 7.0797e+05
D120 = 5.2872e+05
D660 = 5.8252e+05
D160 = 0.0
D260 = 0.0

K110 = 3.4488e+03
K220 = 3.4488e+03

meshSize = 4
lambdaEnd = 2.5
evalAtQuater = True 

runSimulation(a0,b0,
              A110,A120,A220,A660,
              D110,D120,D160,D220,D260,D660,
              meshSize,lambdaEnd,evalAtQuater,K110,K220)