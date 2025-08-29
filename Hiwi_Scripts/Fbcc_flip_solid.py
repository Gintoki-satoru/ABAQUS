# -*- coding: mbcs -*-
from abaqus import *
from abaqusConstants import *
import __main__

def Macro2():
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
    import csv
    import odbAccess
    import re
    import os

    model_names = ['12mm_solid', '16mm_solid', '20mm_solid', '24mm_solid']
    depth_values = [round(d, 2) for d in frange(0.5, 3.5, 0.25)]
    bc_types = ['SS', 'CC']
    csv_path = os.path.join(os.getcwd(), "buckling_results.csv")

    with open(csv_path, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Model Name', 'Core Thickness (mm)', 'Facesheet Thickness (mm)', 'BC Type', 'Buckling Load (N)'])

        for model_name in model_names:
            core_thickness = int(re.match(r"(\d+)", model_name).group(1))
            for bc_type in bc_types:
                for depth in depth_values:

                    # Set facesheet extrusion depth
                    p = mdb.models[model_name].parts['facesheet']
                    p.features['Solid extrude-1'].setValues(depth=depth)
                    p.regenerate()

                    a1 = mdb.models[model_name].rootAssembly
                    a1.regenerate()
                    session.viewports['Viewport: 1'].setValues(displayedObject=a1)

                    # Apply boundary conditions based on bc_type
                    if bc_type == 'SS':
                        mdb.models[model_name].boundaryConditions['left'].setValues(
                            u1=SET, u2=SET, u3=SET, ur1=SET, ur2=UNSET, ur3=UNSET)
                        mdb.models[model_name].boundaryConditions['right'].setValues(
                            u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
                    elif bc_type == 'CC':
                        mdb.models[model_name].boundaryConditions['left'].setValues(
                            u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
                        mdb.models[model_name].boundaryConditions['right'].setValues(
                            u1=UNSET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)

                    # Mesh part
                    p = mdb.models[model_name].parts['facesheet']
                    session.viewports['Viewport: 1'].setValues(displayedObject=p)
                    p.generateMesh()
                    a1 = mdb.models[model_name].rootAssembly
                    a1.regenerate()

                    # Submit job
                    job_name = "{}_{}_{}_AK".format(model_name, str(depth).replace('.', '_'), bc_type)
                    mdb.Job(name=job_name, model=model_name, description='', type=ANALYSIS,
                            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
                            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                            scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1,
                            multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

                    mdb.jobs[job_name].submit()
                    mdb.jobs[job_name].waitForCompletion()

                    odb_path = job_name + '.odb'
                    if os.path.exists(odb_path):
                        odb = odbAccess.openOdb(path=odb_path)
                        try:
                            step_obj = odb.steps['Step-1']
                            eig_vals = step_obj.frames[1].description
                            eigenvalue = float(eig_vals.split('=')[1])
                            eigenvalue = eigenvalue / 24.0

                        except:
                            eigenvalue = 'ERROR'
                        odb.close()
                    else:
                        eigenvalue = 'NO_ODB'

                    writer.writerow([model_name, core_thickness, depth, bc_type, eigenvalue])

def frange(start, stop, step):
    while start <= stop + 1e-6:
        yield round(start, 6)
        start += step

Macro2()
