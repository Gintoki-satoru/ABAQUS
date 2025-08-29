from abaqus import *
from abaqusConstants import *
import __main__

import section, job, step, material, assembly
import os
import csv

def Macro2():
    import visualization
    import odbAccess

    # Define your model names
    model_names = ['12mm_solid-280']
    thickness_values = [round(x * 0.1, 2) for x in range(5, 31)]  # 0.5 to 3.0

    results = []

    for model_name in model_names:
        if model_name not in mdb.models.keys():
            print("Model {} not found.".format(model_name))
            continue

        section_name = 'face'
        material_name = 'AL'

        # Extract length from model name
        try:
            length = int(model_name.split('-')[-1])
        except:
            print("Could not parse length from model name {}".format(model_name))
            length = None

        csv_file = 'buckling_results.csv'

        if not os.path.exists(csv_file):
            with open(csv_file, 'wb') as f:
                writer = csv.writer(f)
                writer.writerow(['Thickness', 'Length', 'BucklingLoad'])

        for thickness in thickness_values:
            # Set thickness
            mdb.models[model_name].sections[section_name].setValues(
                preIntegrate=OFF,
                material=material_name,
                thicknessType=UNIFORM,
                thickness=thickness,
                thicknessField='',
                nodalThicknessField='',
                idealization=NO_IDEALIZATION,
                integrationRule=SIMPSON,
                numIntPts=27
            )

            a = mdb.models[model_name].rootAssembly
            vp = session.viewports[session.currentViewportName]
            vp.setValues(displayedObject=a)

            job_name = '{}_{}'.format(model_name.split('-')[-1], str(thickness).replace('.', '_'))

            if job_name in mdb.jobs.keys():
                del mdb.jobs[job_name]

            mdb.Job(
                name=job_name,
                model=model_name,
                type=ANALYSIS,
                memory=90,
                memoryUnits=PERCENTAGE,
                getMemoryFromAnalysis=True,
                resultsFormat=ODB,
                multiprocessingMode=DEFAULT,
                numCpus=1
            )

            # Run the job
            mdb.jobs[job_name].submit()
            mdb.jobs[job_name].waitForCompletion()

            # Open the ODB and read the first buckling eigenvalue
            odb_path = job_name + '.odb'
            if os.path.exists(odb_path):
                odb = odbAccess.openOdb(path=odb_path)
                try:
                    step_obj = odb.steps['Step-1']
                    eig_vals = step_obj.frames[1].description
                    eigenvalue = float(eig_vals.split('=')[1])
                except:
                    eigenvalue = 'ERROR'
                odb.close()
            else:
                eigenvalue = 'NO_ODB'

            print('Thickness: {}, Length: {}, Buckling Load: {}'.format(thickness, length, eigenvalue))
            with open(csv_file, 'ab') as f:
                writer = csv.writer(f)
                writer.writerow([thickness, length, eigenvalue])

Macro2()
