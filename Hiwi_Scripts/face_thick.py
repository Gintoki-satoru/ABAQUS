# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import visualization
import csv

min_thick = 0.5
max_thick = 3

def face_thickness(min_thick, max_thick):
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

    # Define the range of thickness values
    lower = int(min_thick / 0.25)
    upper = int(max_thick / 0.25) + 1
    thickness_values = [i * 0.25 for i in range(lower, upper)]

    # CSV file to write results
    with open("Results.csv", mode="wb") as file:
        writer = csv.writer(file)
        writer.writerow(["Thickness (mm)", "Mass (g)", "Eigenvalue"])

        for thickness in thickness_values:

            mdb.models['L_model-4mm-Corrugated'].sections['Al_face'].setValues(
                preIntegrate=OFF,
                material='Al',
                thicknessType=UNIFORM,
                thickness=thickness,
                thicknessField='',
                nodalThicknessField='',
                idealization=NO_IDEALIZATION,
                integrationRule=SIMPSON,
                numIntPts=5
            )

            # Meshing of facesheet
            '''p = mdb.models['L_model-4mm-corrugated'].parts['FACE']
            p.deleteMesh()
            p.seedPart(size=2.0, deviationFactor=0.1, minSizeFactor=0.1)
            p.generateMesh()
            a1 = mdb.models['L_model-4mm-corrugated'].rootAssembly
            a1.regenerate()'''

            # Get the mass of the assembly
            a1 = mdb.models['L_model-4mm-Corrugated'].rootAssembly
            mass_properties = a1.getMassProperties()
            mass = mass_properties['mass']

            '''model = mdb.models['L_model-4mm-twinsheet']
            part = model.parts['twin']

            # Get mass properties
            mass_properties = part.getMassProperties()
            mass = mass_properties['mass']'''



            # Create and submit a job
            job_name = 'core-{:.2f}'.format(thickness).replace('.', '_')
            mdb.Job(name=job_name, model='L_model-4mm-Corrugated', description='', type=ANALYSIS, atTime=None,
                waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE,
                getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0
            )
            mdb.jobs[job_name].submit(consistencyChecking=OFF)
            mdb.jobs[job_name].waitForCompletion()

            # Use the .odb file to extract the eigenvalue
            odb_path = "{}.odb".format(job_name)
            odb = visualization.openOdb(path=odb_path)
            buckling_step = odb.steps['Step-1']
            eigenvalue_description = buckling_step.frames[-1].description
            eigenvalue = float(eigenvalue_description.split("=")[-1].strip())

            writer.writerow([thickness, mass, eigenvalue])

            odb.close()

# Run the macro
face_thickness(min_thick, max_thick)