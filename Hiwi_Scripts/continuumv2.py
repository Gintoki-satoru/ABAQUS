# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

min_ply = 5
max_ply = 6
ply_thickness = 0.375
mesh_size = 0.5

def setup_coupling_constraints(num_plies, ply_thickness):
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

    # Delete existing constraints
    # mdb.models['Continuum'].constraints.delete(('Constraint-3', 'Constraint-4', ))

    # Parameters

    a = mdb.models['Continuum'].rootAssembly
    s1 = a.instances['face-1'].faces
    s2 = a.instances['face-2'].faces
    c3 = a.instances['core_12mm-1'].edges

    # Generate ply interface positions
    top_face_start = 12  # Top facesheet starts above the core
    bottom_face_start = 0.0  # Bottom facesheet starts at the core

    # Generate dynamic side1Faces for face-1
    side1Faces1_top = []
    side1Faces1_bottom = []
    for i in range(num_plies):
        y_coord_top = top_face_start + 0.25 + (i * ply_thickness)
        y_coord_bottom = bottom_face_start - 0.25 - (i * ply_thickness)
        face_top = s1.findAt(((0.0, y_coord_top, 8.0),))
        face_bottom = s2.findAt(((0.0, y_coord_bottom, 8.0),))

        # Append the actual face objects
        side1Faces1_top.append(face_top)
        side1Faces1_bottom.append(face_bottom)


    # Generate dynamic side1Faces for face-2
    side1Faces2_top = []
    side1Faces2_bottom = []
    for i in range(num_plies):
        y_coord_top = top_face_start + 0.25  + (i * ply_thickness)
        y_coord_bottom = bottom_face_start - 0.25  - (i * ply_thickness)
        face_top = s1.findAt(((200.0, y_coord_top, 8.0),))
        face_bottom = s2.findAt(((200.0, y_coord_bottom, 8.0),))

        # Append the actual face objects
        side1Faces2_top.append(face_top)
        side1Faces2_bottom.append(face_bottom)

    # Define circumferential edges
    circumEdges3 = c3.findAt(
        ((0.0, 11.0, 24.0),), ((0.0, 7.0, 24.0),), ((0.0, 11.0, 20.0),), ((0.0, 7.0, 20.0),),
        ((0.0, 11.0, 16.0),), ((0.0, 3.0, 24.0),), ((0.0, 3.0, 16.0),), ((0.0, 7.0, 12.0),),
        ((0.0, 11.0, 8.0),), ((0.0, 3.0, 20.0),), ((0.0, 3.0, 0.0),), ((0.0, 7.0, 16.0),),
        ((0.0, 11.0, 12.0),), ((0.0, 3.0, 12.0),), ((0.0, 3.0, 8.0),), ((0.0, 3.0, 4.0),),
        ((0.0, 7.0, 0.0),), ((0.0, 11.0, 4.0),), ((0.0, 7.0, 8.0),), ((0.0, 7.0, 4.0),),
        ((0.0, 11.0, 0.0),)
    )

    # Combine all faces

    region2_left = regionToolset.Region(side1Faces=side1Faces1_top + side1Faces1_bottom,
                                        circumEdges=circumEdges3)

    c3 = a.instances['core_12mm-1'].edges
    circumEdges3 = c3.findAt(((200.0, 11.0, 20.0),), ((200.0, 11.0, 24.0),), ((200.0, 7.0, 20.0),),
                             ((200.0, 7.0, 24.0),), ((200.0, 3.0, 4.0),), ((200.0, 7.0, 4.0),),
                             ((200.0, 11.0, 16.0),), ((200.0, 3.0, 24.0),), ((200.0, 3.0, 12.0),),
                             ((200.0, 7.0, 12.0),), ((200.0, 11.0, 8.0),), ((200.0, 7.0, 8.0),),
                             ((200.0, 3.0, 0.0),), ((200.0, 7.0, 0.0),), ((200.0,3.0, 20.0),), ((200.0, 7.0, 16.0),),
                             ((200.0, 3.0, 8.0),), ((200.0,11.0, 0.0),), ((200.0, 11.0, 4.0),), ((200.0, 3.0, 16.0),),
                             ((200.0, 11.0, 12.0),))
    region2_right = regionToolset.Region(side1Faces=side1Faces2_top + side1Faces2_bottom,
                                         circumEdges=circumEdges3)

    # Reference Points at both ends
    r1 = a.referencePoints
    refPoint_left = (r1[671],)
    refPoint_right = (r1[672],)

    region1_left = regionToolset.Region(referencePoints=refPoint_left)
    region1_right = regionToolset.Region(referencePoints=refPoint_right)

    # Coupling Constraints for both ends
    mdb.models['Continuum'].Coupling(
        name='Constraint-Left',
        controlPoint=region1_left,
        surface=region2_left,
        influenceRadius=WHOLE_SURFACE,
        couplingType=KINEMATIC,
        alpha=0.0,
        localCsys=None,
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
    )

    mdb.models['Continuum'].Coupling(
        name='Constraint-Right',
        controlPoint=region1_right,
        surface=region2_right,
        influenceRadius=WHOLE_SURFACE,
        couplingType=KINEMATIC,
        alpha=0.0,
        localCsys=None,
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
    )


def parametric_study(min_plies, max_plies, thick, m_size):
    import section
    import regionToolset
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import job
    import sketch
    import visualization

    for num_plies in range(min_plies, max_plies):
        total_thickness = num_plies * thick

        # Delete existing offset faces
        p = mdb.models['Continuum'].parts['face']
        p.deleteFeatures(('Offset faces-1', 'Offset faces-2', 'Offset faces-3', 'Offset faces-4', 'Offset faces-5',
                          'Offset faces-6', 'Offset faces-7', 'Offset faces-8'))

        # Increase facesheet thickness
        p.features['Solid extrude-1'].setValues(depth=total_thickness)
        p.regenerate()

        # Offset faces for ply definition using coordinate system
        f = p.faces
        for i in range(num_plies - 1):
            ply_offset = (i + 1) * thick  # Calculate offset for each ply
            p.OffsetFaces(faceList=(f.findAt(coordinates=(-33.333333, -9.0, 0.0)),), distance=ply_offset,
                          trimToReferenceRep=False)

        # Define local coordinate system for layup orientation
        layupOrientation = p.datums[5]
        c = p.cells
        regions = []
        for i in range(num_plies):
            # Define coordinates for the cell at each ply depth
            ply_coords = (-100.0, -17.0, 0.25 + i * thick)
            cells = c.findAt((ply_coords,))
            regions.append(regionToolset.Region(cells=cells))

        # Define composite layup
        compositeLayup = p.CompositeLayup(name='CompositeLayup-1', elementType=CONTINUUM_SHELL, symmetric=False)
        compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, poissonDefinition=DEFAULT)

        # Define reference orientation for the layup
        compositeLayup.ReferenceOrientation(orientationType=SYSTEM, localCsys=layupOrientation, axis=AXIS_3,
                                            stackDirection=STACK_ORIENTATION)

        # Assign plies with alternating orientation (0/90)
        for i in range(num_plies):
            orientation = 0.0 if i % 2 == 0 else 90.0
            compositeLayup.CompositePly(plyName='Ply-{}'.format(i + 1), region=regions[i], material='carbon_epoxy',
                                        thicknessType=SPECIFY_THICKNESS, thickness=0.375,
                                        orientationType=SPECIFY_ORIENT, orientationValue=orientation, axis=AXIS_3,
                                        angle=0.0, numIntPoints=3)

        compositeLayup.resume()

        # Meshing
        p.seedPart(size=m_size, deviationFactor=0.1, minSizeFactor=0.1)

        elemType1 = mesh.ElemType(elemCode=SC8R, elemLibrary=STANDARD, hourglassControl=ENHANCED)
        elemType2 = mesh.ElemType(elemCode=SC6R, elemLibrary=STANDARD)
        elemType3 = mesh.ElemType(elemCode=UNKNOWN_TET, elemLibrary=STANDARD)

        # Set the element types for all cells
        pickedRegions = (p.cells,)
        p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))

        # Remove the edge seeds
        e = p.edges
        for edge in e:
            pickedEdges = (edge,)
            p.deleteSeeds(regions=pickedEdges)

        # Generate the mesh after removing the edge seeds
        p.generateMesh()

        # Apply Coupling Constraints
        setup_coupling_constraints(num_plies, thick)
        p.regenerate()

        # Create job for this ply configuration
        '''job_name = 'CC_{}ply'.format(num_plies)
        mdb.Job(name=job_name, model='Continuum', type=ANALYSIS, memory=90, memoryUnits=PERCENTAGE,
                getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, resultsFormat=ODB,
                numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
        mdb.jobs[job_name].submit(consistencyChecking=OFF)
        mdb.jobs[job_name].waitForCompletion()'''

# Run the parametric study
parametric_study(min_ply, max_ply, ply_thickness, mesh_size)
