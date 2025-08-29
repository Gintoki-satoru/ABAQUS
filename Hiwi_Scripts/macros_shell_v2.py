# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

def macro_shell(min_plies, max_plies, thick):
    import section
    import regionToolset
    import part
    import material
    import assembly
    import step
    import job
    import sketch
    import visualization
    import connectorBehavior

    # Get the part and face region
    p = mdb.models['Continuum'].parts['face']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(faces=faces)

    # Define ply orientations (alternating 0° and 90°)
    orientations = [0.0, 90.0]

    for num_plies in range(min_plies, max_plies + 1):

        if 'CompositeLayup-1' in p.compositeLayups.keys():
            del p.compositeLayups['CompositeLayup-1']

        # Create a new composite layup
        compositeLayup = p.CompositeLayup(name='CompositeLayup-1', description='')
        compositeLayup.orientation.setValues(additionalRotationType=ROTATION_NONE, angle=0.0)
        compositeLayup.setValues(    
            offsetType=BOTTOM_SURFACE)

        compositeLayup.ReferenceOrientation(
            orientationType=GLOBAL, 
            localCsys=None, 
            fieldName='', 
            additionalRotationType=ROTATION_NONE, 
            angle=0.0, 
            axis=AXIS_3
        )

        # Add plies to the new layup
        for i in range(num_plies):
            ply_name = 'Ply-' + str(i + 1)
            orientation_value = orientations[i % 2]  # Alternates between 0° and 90°

            compositeLayup.CompositePly(
                suppressed=False,
                plyName=ply_name,
                region=region,
                material='carbon_epoxy',
                thicknessType=SPECIFY_THICKNESS,
                thickness=thick,
                orientationType=SPECIFY_ORIENT,
                orientationValue=orientation_value,
                additionalRotationType=ROTATION_NONE,
                additionalRotationField='',
                axis=AXIS_3,
                angle=0.0,
                numIntPoints=3
            )

        compositeLayup.resume()
        a = mdb.models['Continuum'].rootAssembly
        a.regenerate()

        # Submit job
        job_name = 'Composite_' + str(num_plies) + 'ply'
        mdb.Job(name=job_name, model='Continuum', description='', type=ANALYSIS,
                atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
                memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
        

        mdb.jobs[job_name].submit()
        mdb.jobs[job_name].waitForCompletion()

        print('Job ' + job_name + ' completed!')

# Define input parameters
min_ply = 4
max_ply = 8
ply_thickness = 0.375

macro_shell(min_ply, max_ply, ply_thickness)