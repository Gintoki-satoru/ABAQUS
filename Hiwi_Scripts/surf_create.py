# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
from collections import OrderedDict
from abaqus import session
import __main__

min_ply = 1
max_ply = 5
ply_thickness = 0.375
mesh_size = 0.375
bottom_face_start = 0.0
core_thickness_list = [12]


def create_surface(num_plies, ply_thickness, core_thickness, bottom_face_start):
    import regionToolset
    import part

    core_layer_thickness = 4
    core_layers = core_thickness // core_layer_thickness

    core_name = 'core_{}mm'.format(core_thickness)

    a = mdb.models['Continuum'].rootAssembly
    a.deleteSurfaces(surfaceNames=('Surf-right', 'Surf-left', 'Surf-bottom', 'Surf-top',))
    s1 = a.instances['face-1'].faces
    s2 = a.instances['face-2'].faces
    c3 = a.instances['{}-1'.format(core_name)].edges

    side1Faces1_top = []
    side1Faces1_bottom = []
    side1Faces2_top = []
    side1Faces2_bottom = []

    for i in range(num_plies):
        y_coord_top = core_thickness + 0.25 + (i * ply_thickness)
        y_coord_bottom = bottom_face_start - 0.25 - (i * ply_thickness)

        # Find the faces for the top and bottom sheets
        face_top = s1.findAt(((0.0, y_coord_top, 8.0),))
        face_bottom = s2.findAt(((0.0, y_coord_bottom, 8.0),))

        # Append the face objects for face-1
        side1Faces1_top.append(face_top)
        side1Faces1_bottom.append(face_bottom)


    for i in range(num_plies):
        y_coord_top = core_thickness + 0.25 + (i * ply_thickness)
        y_coord_bottom = bottom_face_start - 0.25 - (i * ply_thickness)

        # Find the faces for the top and bottom sheets
        face_top = s1.findAt(((200.0, y_coord_top, 8.0),))
        face_bottom = s2.findAt(((200.0, y_coord_bottom, 8.0),))

        # Append the face objects for face-2
        side1Faces2_top.append(face_top)
        side1Faces2_bottom.append(face_bottom)


    combined_faces1 = side1Faces1_top + side1Faces1_bottom
    combined_faces2 = side1Faces2_top + side1Faces2_bottom

    circumEdges3 = []
    circumEdges4 = []

    for i in range(core_layers):
        y_coord = 3 + i * 4

        edge1 = c3.findAt(((0.0, y_coord, 24.0),))
        edge2 = c3.findAt(((0.0, y_coord, 20.0),))
        edge3 = c3.findAt(((0.0, y_coord, 16.0),))
        edge4 = c3.findAt(((0.0, y_coord, 12.0),))
        edge5 = c3.findAt(((0.0, y_coord, 8.0),))
        edge6 = c3.findAt(((0.0, y_coord, 4.0),))
        edge7 = c3.findAt(((0.0, y_coord, 0.0),))

        circumEdges3.append(edge1)
        circumEdges3.append(edge2)
        circumEdges3.append(edge3)
        circumEdges3.append(edge4)
        circumEdges3.append(edge5)
        circumEdges3.append(edge6)
        circumEdges3.append(edge7)

    for i in range(core_layers):
        y_coord = 3 + i * 4

        edge1 = c3.findAt(((200.0, y_coord, 24.0),))
        edge2 = c3.findAt(((200.0, y_coord, 20.0),))
        edge3 = c3.findAt(((200.0, y_coord, 16.0),))
        edge4 = c3.findAt(((200.0, y_coord, 12.0),))
        edge5 = c3.findAt(((200.0, y_coord, 8.0),))
        edge6 = c3.findAt(((200.0, y_coord, 4.0),))
        edge7 = c3.findAt(((200.0, y_coord, 0.0),))

        circumEdges4.append(edge1)
        circumEdges4.append(edge2)
        circumEdges4.append(edge3)
        circumEdges4.append(edge4)
        circumEdges4.append(edge5)
        circumEdges4.append(edge6)
        circumEdges4.append(edge7)

    # Create the surface
    a.Surface(side1Faces=combined_faces1, circumEdges=circumEdges3, name='Surf-left')
    a.Surface(side1Faces=combined_faces2, circumEdges=circumEdges4, name='Surf-right')

    p = mdb.models['Continuum'].parts[core_name]
    v = p.vertices

    mdb.models['Continuum'].parts[core_name].deleteSets(setNames=(
        'CoreVertices - 1', 'CoreVertices - 2',))

    vertex_points1 = []
    vertex_points2 = []

    for x in range(0, 201, 4):
        for z in range(0, 25, 4):
            vertex = v.findAt(((x, core_thickness, z),))
            vertex_points1.append(vertex)

    # Create the set
    p.Set(vertices=vertex_points1, name='CoreVertices - 1')

    for x in range(0, 201, 4):
        for z in range(0, 25, 4):
            vertex = v.findAt(((x, bottom_face_start, z),))
            vertex_points2.append(vertex)

    # Create the set
    p.Set(vertices=vertex_points2, name='CoreVertices - 2')

    a = mdb.models['Continuum'].rootAssembly
    s1 = a.instances['face-1'].faces
    side1Faces1 = s1.findAt(((100, core_thickness, 8.0),))
    a.Surface(side1Faces=side1Faces1, name='Surf-top')

    s1 = a.instances['face-2'].faces
    side1Faces1 = s1.findAt(((100, bottom_face_start, 8),))
    a.Surface(side1Faces=side1Faces1, name='Surf-bottom')




def parametric_study(min_plies, max_plies, thick, m_size, core_thickness_list):
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
    import time

    for core_thickness in core_thickness_list:
        core_name = 'core_{}mm'.format(core_thickness)

        # Modify assembly for different core thickness

        a = mdb.models['Continuum'].rootAssembly
        a.suppressFeatures(('face-1', 'face-2',))

        if 'core_12mm-1' in a.features:
            del a.features['core_12mm-1']
        if 'core_16mm-1' in a.features:
            del a.features['core_16mm-1']
        if 'core_20mm-1' in a.features:
            del a.features['core_20mm-1']
        if 'core_24mm-1' in a.features:
            del a.features['core_24mm-1']


        a1 = mdb.models['Continuum'].rootAssembly
        p = mdb.models['Continuum'].parts[core_name]
        new_core =  '{}-1'.format(core_name)
        a1.Instance(name=new_core, part=p, dependent=ON)
        a = mdb.models['Continuum'].rootAssembly
        a.resumeFeatures(('face-1', 'face-2',))

        a = mdb.models['Continuum'].rootAssembly

        del a.features['face-1']
        a1 = mdb.models['Continuum'].rootAssembly
        p = mdb.models['Continuum'].parts['face']
        a1.Instance(name='face-1', part=p, dependent=ON)

        a1 = mdb.models['Continuum'].rootAssembly
        a1.rotate(instanceList=('face-1',), axisPoint=(0.0, 0.0, 0.0), axisDirection=(
            10.0, 0.0, 0.0), angle=-90.0)

        a1.translate(instanceList=('face-1',), vector=(100.0, core_thickness, -1.0))

        boundary_conditions = OrderedDict([
            ('CS', {
                'BC-1': {'u1': SET, 'u2': SET, 'u3': SET, 'ur1': SET, 'ur2': SET, 'ur3': SET},  # Clamped end
                'BC-2': {'u1': UNSET, 'u2': SET, 'u3': SET, 'ur1': UNSET, 'ur2': UNSET, 'ur3': UNSET}
                # Simply supported end
            }),
            ('CC', {
                'BC-1': {'u1': SET, 'u2': SET, 'u3': SET, 'ur1': SET, 'ur2': SET, 'ur3': SET},  # Clamped end
                'BC-2': {'u1': UNSET, 'u2': SET, 'u3': SET, 'ur1': SET, 'ur2': SET, 'ur3': SET}  # Clamped end
            }),
            ('SS', {
                'BC-1': {'u1': SET, 'u2': SET, 'u3': SET, 'ur1': SET, 'ur2': UNSET, 'ur3': UNSET},
                # Simply supported end
                'BC-2': {'u1': UNSET, 'u2': SET, 'u3': SET, 'ur1': UNSET, 'ur2': UNSET, 'ur3': UNSET}
                # Simply supported end
            })
        ])

        for case_name, bc_set in boundary_conditions.items():
            print("Processing {}".format(case_name))
            for num_plies in range(min_plies, max_plies):
                total_thickness = num_plies * thick

                p = mdb.models['Continuum'].parts['face']
                p.deleteFeatures(('Offset faces-1', 'Offset faces-2', 'Offset faces-3', 'Offset faces-4', 'Offset faces-5',
                                  'Offset faces-6', 'Offset faces-7', 'Offset faces-8', 'Offset faces-9', 'Offset faces-10'))

                # Increase facesheet thickness
                p.features['Solid extrude-1'].setValues(depth=total_thickness)
                p.regenerate()

                # Offset faces for ply definition using coordinate system
                f = p.faces
                for i in range(num_plies - 1):
                    ply_offset = (i + 1) * thick
                    p.OffsetFaces(faceList=(f.findAt(coordinates=(-33.333333, -9.0, 0.0)),), distance=ply_offset,
                                  trimToReferenceRep=False)

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

                pickedRegions = (p.cells,)
                p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))

                e = p.edges
                for edge in e:
                    pickedEdges = (edge,)
                    p.deleteSeeds(regions=pickedEdges)

                p.generateMesh()

                # Apply Coupling Constraints
                create_surface(num_plies, ply_thickness, core_thickness, bottom_face_start)

                constraint_names = ('Constraint-leftRP', 'Constraint-rightRP', 'Constraint-top', 'Constraint-bottom')

                for name in constraint_names:
                    if name in mdb.models['Continuum'].constraints.keys():
                        mdb.models['Continuum'].constraints.delete((name,))

                a = mdb.models['Continuum'].rootAssembly
                a.deleteFeatures(('RP-1', 'RP-2',))

                refPoint1 = a.ReferencePoint((0.0, core_thickness / 2.0, 12.0))
                refPoint1_obj = a.referencePoints[refPoint1.id]

                region1 = regionToolset.Region(referencePoints=(refPoint1_obj,))
                a = mdb.models['Continuum'].rootAssembly
                region2 = a.surfaces['Surf-left']
                mdb.models['Continuum'].Coupling(name='Constraint-leftRP', controlPoint=region1,
                                                 surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                                                 alpha=0.0, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
                a = mdb.models['Continuum'].rootAssembly
                refPoint2 = a.ReferencePoint((200.0, core_thickness / 2.0, 12.0))

                refPoint2_obj = a.referencePoints[refPoint2.id]
                region1 = regionToolset.Region(referencePoints=(refPoint2_obj,))
                a = mdb.models['Continuum'].rootAssembly
                region2 = a.surfaces['Surf-right']
                mdb.models['Continuum'].Coupling(name='Constraint-rightRP', controlPoint=region1,
                                                 surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                                                 alpha=0.0, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

                a = mdb.models['Continuum'].rootAssembly
                region1 = a.surfaces['Surf-top']
                a = mdb.models['Continuum'].rootAssembly
                region2 = a.instances['{}-1'.format(core_name)].sets['CoreVertices - 1']
                mdb.models['Continuum'].Tie(name='Constraint-top', main=region1,
                                            secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON,
                                            tieRotations=ON, thickness=ON)
                a = mdb.models['Continuum'].rootAssembly
                region1 = a.surfaces['Surf-bottom']
                a = mdb.models['Continuum'].rootAssembly
                region2 = a.instances['{}-1'.format(core_name)].sets['CoreVertices - 2']
                mdb.models['Continuum'].Tie(name='Constraint-bottom', main=region1,
                                            secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON,
                                            tieRotations=ON, thickness=ON)
                p.regenerate()

                # Apply load
                if 'Load-1' in mdb.models['Continuum'].loads.keys():
                    del mdb.models['Continuum'].loads['Load-1']
                region = regionToolset.Region(referencePoints=(refPoint2_obj,))
                mdb.models['Continuum'].ConcentratedForce(name='Load-1',
                                                          createStepName='Step-1', region=region, cf1=-1.0,
                                                          distributionType=UNIFORM,
                                                          field='', localCsys=None)

                # Apply BC
                for bc_name, bc_values in bc_set.items():
                    region = regionToolset.Region(
                        referencePoints=(refPoint1_obj if bc_name == 'BC-1' else refPoint2_obj,))

                    if bc_name in mdb.models['Continuum'].boundaryConditions.keys():
                        del mdb.models['Continuum'].boundaryConditions[bc_name]

                    mdb.models['Continuum'].DisplacementBC(name=bc_name, createStepName='Initial',
                                                           region=region,
                                                           u1=bc_values['u1'],
                                                           u2=bc_values['u2'],
                                                           u3=bc_values['u3'],
                                                           ur1=bc_values['ur1'],
                                                           ur2=bc_values['ur2'],
                                                           ur3=bc_values['ur3'],
                                                           amplitude=UNSET, distributionType=UNIFORM, fieldName='',
                                                           localCsys=None)

                # Create job for this ply configuration
                '''job_name = 'comp_{}mm_{}_{}ply'.format(core_thickness, case_name, num_plies)
                mdb.Job(name=job_name, model='Continuum', type=ANALYSIS, memory=90, memoryUnits=PERCENTAGE,
                        getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, resultsFormat=ODB,
                        numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
                mdb.jobs[job_name].submit(consistencyChecking=OFF)
                mdb.jobs[job_name].waitForCompletion()'''


# Run the parametric study
parametric_study(min_ply, max_ply, ply_thickness, mesh_size, core_thickness_list)
