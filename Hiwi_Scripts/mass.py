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

    # Define the range of thickness values for the face
    lower = int(min_thick / 0.25)
    upper = int(max_thick / 0.25) + 1
    face_thickness_values = [i * 0.25 for i in range(lower, upper)]

    # Define the thickness values for honeycomb_thin
    honeycomb_thin_values = [0.0231, 0.034, 0.0473, 0.0615, 0.075]

    # List of models to run
    model_names = ['L_model-12mm']

    # CSV file to write results
    with open("Results.csv", mode="wb") as file: 
        writer = csv.writer(file)
        writer.writerow(["Model", "Face Thickness (mm)", "Honeycomb Thin (mm)", "Honeycomb Thick (mm)", "Mass (kg)"])  # Removed "Eigenvalue" from headers

        job_name = "reusable_job"

        for model_name in model_names:
            for face_thickness_value in face_thickness_values:
                # Update the facesheet thickness
                mdb.models[model_name].sections['Al_face'].setValues(
                    preIntegrate=OFF,
                    material='Al',
                    thicknessType=UNIFORM,
                    thickness=face_thickness_value,
                    thicknessField='',
                    nodalThicknessField='',
                    idealization=NO_IDEALIZATION,
                    integrationRule=SIMPSON,
                    numIntPts=5
                )

                for honeycomb_thin_value in honeycomb_thin_values:
                    mdb.models[model_name].sections['honeycomb_thin'].setValues(
                        preIntegrate=OFF,
                        material='Al',
                        thicknessType=UNIFORM,
                        thickness=honeycomb_thin_value,
                        thicknessField='',
                        nodalThicknessField='',
                        idealization=NO_IDEALIZATION,
                        integrationRule=SIMPSON,
                        numIntPts=5
                    )
                    honeycomb_thick_value = 2 * honeycomb_thin_value
                    mdb.models[model_name].sections['honeycomb_thick'].setValues(
                        preIntegrate=OFF,
                        material='Al',
                        thicknessType=UNIFORM,
                        thickness=honeycomb_thick_value,
                        thicknessField='',
                        nodalThicknessField='',
                        idealization=NO_IDEALIZATION,
                        integrationRule=SIMPSON,
                        numIntPts=5
                    )
                    a = mdb.models[model_name].rootAssembly
                    a.regenerate()
                    a1 = mdb.models[model_name].rootAssembly
                    mass_properties = a1.getMassProperties()
                    mass = mass_properties['mass']

                    # Meshing of facesheet (commented out)
                    # p = mdb.models[model_name].parts['FACE']
                    # p.deleteMesh()
                    # p.seedPart(size=2.0, deviationFactor=0.1, minSizeFactor=0.1)
                    # p.generateMesh()
                    # a1 = mdb.models[model_name].rootAssembly
                    # a1.regenerate()

                    # Delete the existing job if it exists
                    # if job_name in mdb.jobs:
                      #  del mdb.jobs[job_name]

                    # Create and submit the job with the fixed name
                    # mdb.Job(name=job_name, model=model_name, description='', type=ANALYSIS, atTime=None,
                    #    waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE,
                     #   getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                     #   modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                     #   scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
                    # mdb.jobs[job_name].submit(consistencyChecking=OFF)
                    # mdb.jobs[job_name].waitForCompletion()

                    # Open the .odb file to extract the eigenvalue (commented out)
                    # odb_path = "{}.odb".format(job_name)  # Corrected line
                    # odb = visualization.openOdb(path=odb_path)
                    # buckling_step = odb.steps['BucklingStep']  # Replace 'BucklingStep' with your step name
                    # eigenvalue_description = buckling_step.frames[-1].description  # Extract eigenvalue description

                    # Extract only the numerical value from the eigenvalue description (commented out)
                    # eigenvalue = float(eigenvalue_description.split("=")[-1].strip())  # Extract "800.64" from "Mode 1: EigenValue = 800.64"

                    # Write results to CSV
                    writer.writerow([
                        model_name,
                        face_thickness_value,
                        honeycomb_thin_value,
                        honeycomb_thick_value,
                        mass,


                    # odb.close()

# Run the macro
face_thickness(min_thick, max_thick)