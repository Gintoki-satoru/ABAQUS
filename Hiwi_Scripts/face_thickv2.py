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

    # Range of thickness values for the face
    lower = int(min_thick / 0.25)
    upper = int(max_thick / 0.25) + 1
    face_thickness_values = [i * 0.25 for i in range(lower, upper)]

    honeycomb_thin_values = [0.0231, 0.034, 0.0473, 0.0615, 0.075]

    # List of models to run
    model_names = ['L_model-14mm', 'L_model-16mm']

    # CSV file to write results
    with open("Results.csv", mode="wb") as file:
        writer = csv.writer(file)
        writer.writerow(["Model", "Face Thickness (mm)", "Wall thickness (mm)", "Node thickness (mm)", "Mass (g)", "Eigenvalue"])

        job_name = "L-MODEL"

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
                    # Update honeycomb_thin thickness
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

                    # Update honeycomb_thick thickness
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

                    # Meshing of facesheet
                    p = mdb.models[model_name].parts['FACE']
                    p.deleteMesh()
                    p.seedPart(size=2.0, deviationFactor=0.1, minSizeFactor=0.1)
                    p.generateMesh()
                    a1 = mdb.models[model_name].rootAssembly
                    a1.regenerate()

                    # Get the mass
                    mass_properties = a1.getMassProperties()
                    mass = mass_properties['mass']

                    if job_name in mdb.jobs:
                        del mdb.jobs[job_name]

                    # submit the job
                    mdb.Job(name=job_name, model=model_name, description='', type=ANALYSIS, atTime=None,
                        waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE,
                        getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0
                    )
                    mdb.jobs[job_name].submit(consistencyChecking=OFF)
                    mdb.jobs[job_name].waitForCompletion()

                    # Open the .odb file to extract the eigenvalue
                    odb_path = "{}.odb".format(job_name)
                    odb = visualization.openOdb(path=odb_path)
                    buckling_step = odb.steps['Step-1']
                    eigenvalue_description = buckling_step.frames[-1].description
                    eigenvalue = float(eigenvalue_description.split("=")[-1].strip())

                    # Write results to CSV
                    writer.writerow([
                        model_name,
                        face_thickness_value,
                        honeycomb_thin_value,
                        honeycomb_thick_value,
                        mass,
                        eigenvalue
                    ])
                    odb.close()

# Run the macro
face_thickness(min_thick, max_thick)