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

    # Get the face region
    p = mdb.models['shell_model'].parts['face_sheet']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(faces=faces)

    # Define ply orientations
    orientations = [0.0, 90.0]


    for num_plies in range(min_plies, max_plies + 1):

        compositeLayup = mdb.models['shell_model'].parts['face_sheet'].compositeLayups['CompositeLayup-1']
        compositeLayup.orientation.setValues(additionalRotationType=ROTATION_NONE, angle=0.0)

        compositeLayup.deletePlies()
        compositeLayup.suppress()


        for i in range(num_plies):
            ply_name = 'Ply-' + str(i+1)
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
                numIntPoints=27
            )

        compositeLayup.resume()
        a = mdb.models['shell_model'].rootAssembly
        a.regenerate()

        # Create and submit job
        job_name = 'Composite_' + str(num_plies) + 'ply'
        mdb.Job(name=job_name, model='shell_model', description='', type=ANALYSIS,
                atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
                memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1,
                multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

        mdb.jobs[job_name].submit()
        mdb.jobs[job_name].waitForCompletion()

        print('Job ' + job_name + ' completed!')



min_ply = 6
max_ply = 8
ply_thickness = 0.375

macro_shell(min_ply, max_ply, ply_thickness)
