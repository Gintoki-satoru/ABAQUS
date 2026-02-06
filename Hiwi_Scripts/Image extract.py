# -*- coding: utf-8 -*-
# Script to extract first 6 frames as PNG images from 16 ODB files in XY (front) view
# Run this inside Abaqus/CAE or Abaqus/Viewer (File > Run Script)

from abaqus import *
from abaqusConstants import *
import os

# --- User paths ---
base_odb_path = r"U:/ADP grp/ADP LSM/Abaqus/Plate buckling analysis/1-16_M10"
save_img_path = r"U:/ADP grp/ADP LSM/Abaqus/Images/M10"

# --- Settings ---
num_models = 16          # LinBuck_1 to LinBuck_16
num_frames = 6           # First 6 frames per ODB
step_name = 'Step-1'     # Step name in each ODB
var_label = 'U'          # Displacement variable
var_component = 'U3'     # Component to plot

# --- Create folder if not exists ---
if not os.path.exists(save_img_path):
    os.makedirs(save_img_path)

vp = session.viewports['Viewport: 1']

for i in range(1, num_models + 1):
    odb_name = "LinBuck_%d.odb" % i
    odb_path = os.path.join(base_odb_path, odb_name)

    print("Processing: %s" % odb_path)
    if not os.path.exists(odb_path):
        print("  [Warning] File not found: %s" % odb_path)
        continue

    # Open ODB
    odb = session.openOdb(name=odb_path)
    vp.setValues(displayedObject=odb)

    # Set display variables
    vp.odbDisplay.setPrimaryVariable(
        variableLabel=var_label,
        outputPosition=NODAL,
        refinement=(COMPONENT, var_component)
    )
    vp.odbDisplay.display.setValues(plotState=CONTOURS_ON_DEF)

    # Loop through first N frames
    for frame in range(1, num_frames + 1):
        try:
            vp.odbDisplay.setFrame(step=step_name, frame=frame)

            # --- Apply front (XY) view ---
            vp.view.setValues(session.views['Front'])
            vp.view.fitView()  # Optional: zoom to fit the model

            img_name = os.path.join(save_img_path, "%d_%d" % (i, frame))
            session.printToFile(
                fileName=img_name,
                format=PNG,
                canvasObjects=(vp,)
            )
            print("  Saved frame %d -> %s.png" % (frame, img_name))
        except Exception as e:
            print("  [Error] Could not export frame %d: %s" % (frame, str(e)))

    odb.close()

print("=== Export complete ===")
