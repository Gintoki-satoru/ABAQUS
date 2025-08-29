# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

min_ply = 5
max_ply = 6
ply_thickness = 0.375
mesh_size = 0.375

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
            ply_offset = (i + 1) * thick
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

        elemType1 = mesh.ElemType(elemCode=SC8R, elemLibrary=STANDARD)
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

        # Generate the mesh
        p.generateMesh()

        # Create job
        job_name = 'CC_{}ply'.format(num_plies)
        mdb.Job(name=job_name, model='Continuum', type=ANALYSIS, memory=90, memoryUnits=PERCENTAGE,
                getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, resultsFormat=ODB,
                numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
        mdb.jobs[job_name].submit(consistencyChecking=OFF)
        mdb.jobs[job_name].waitForCompletion()

# Run the parametric study
parametric_study(min_ply, max_ply, ply_thickness, mesh_size)
