from abaqus import *
from abaqusConstants import *
from odbAccess import *
import sketch
import os
import numpy as np

import section
import regionToolset
import displayGroupMdbToolset as dgm
import mesh
import math
import csv


session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

Mdb()
modelName = 'SuperEllipse'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
model = mdb.models[modelName]
model.rootAssembly.clearGeometryCache()
model.rootAssembly.regenerate()

################ Parameters #####################
a, b, c = 150.0, 150.0, 100.0   # inner semi-axes
total_length = c
thick = 2.5                     # total thickness
n1, n2 = 0.2, 0.1               # shape exponents
num_points = 30                 # points per curve
num_layers = 4                  # number of layers through thickness
num_theta_sections = 3          # number of θ sections(min 2): For even number, the number of partitions created will be (num_theta_sections + 1)
num_partitions = 4              # number of partitions for face BC
pressure_value = 0.1            # pressure magnitude (MPa)
mesh_size = 2.0                 # mesh size

a_out, b_out, c_out = a + thick, b + thick, c + thick  # outer semi-axes

################ Helper Functions #####################
def signed_power(base, exp):
    return math.copysign(abs(base)**exp, base)

def superellipsoid_point_3d(phi, theta, a, b, c, n1, n2):
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    x = a * signed_power(cos_phi, n1) * signed_power(cos_theta, n2)
    y = b * signed_power(cos_phi, n1) * signed_power(sin_theta, n2)
    z = c * signed_power(sin_phi, n1)
    return x, y, z

def create_datums(p, point_list):
    for pt in point_list:
        p.DatumPointByCoordinate(coords=pt)

def create_wire_splines(part, datum_points, smoothClosed=True):
    wire = part.WireSpline(points=datum_points,
                           mergeType=IMPRINT,
                           meshable=ON,
                           smoothClosedSpline=ON if smoothClosed else OFF)
    return wire

def create_layer_part(model, layer_index, num_theta_sections):
    # --- Layer fraction and scaling ---
    frac1 = float(layer_index) / num_layers
    frac2 = float(layer_index + 1) / num_layers
    # Semi-axes for inner and outer surfaces of this layer
    a_i, b_i, c_i = a + frac1*(a_out - a), b + frac1*(b_out - b), c + frac1*(c_out - c)
    a_o, b_o, c_o = a + frac2*(a_out - a), b + frac2*(b_out - b), c + frac2*(c_out - c)
    # Angular values
    phi_vals = [math.radians(i * 90 / num_points) for i in range(num_points + 1)]
    theta_vals = [math.radians(i * 90 / num_points) for i in range(num_points + 1)]
    # --- Determine θ angles for cross-sections ---
    if num_theta_sections == 2:
        theta_sections = [0, 90]
    elif num_theta_sections >= 3:
        theta_sections = [0, 45, 90]
        if num_theta_sections > 3:
            step = 90.0 / (num_theta_sections - 1)
            extra_angles = [
                round(i * step, 2)
                for i in range(1, num_theta_sections - 1)
                if abs(i * step - 45) > 1e-3
            ]
            theta_sections = sorted(set(theta_sections + extra_angles))
        if num_theta_sections > 3:
            if n2 <= 1.0/6:
                theta_sections[1] = 0.01
                theta_sections[2] = 1
                theta_sections[-3] = 90 - 1
                theta_sections[-2] = 90 - 0.01
            elif n2 <= 1.0/3:
                theta_sections[1] = 1
                theta_sections[-2] = 90 - 1
            elif n2 <= 1.0/2:
                theta_sections[1] = 5
                theta_sections[-2] = 90 - 5
            else:
                print("No modification made to theta_sections")
    else:
        raise ValueError("num_theta_sections must be >= 2")
    theta_vals_deg = [math.degrees(t) for t in theta_vals]
    merged_angles = sorted(set(theta_vals_deg + theta_sections))
    theta_vals = [math.radians(r) for r in merged_angles]
    # --- Case 1 (phi=0°) ---
    inner_case1 = [superellipsoid_point_3d(0.0, th, a_i, b_i, c_i, n1, n2) for th in theta_vals]
    outer_case1 = [superellipsoid_point_3d(0.0, th, a_o, b_o, c_o, n1, n2) for th in theta_vals]
    # --- Case 2 (phi=45°) ---
    phi_forfiv = math.radians(45)
    inner_case2 = [superellipsoid_point_3d(phi_forfiv, th, a_i, b_i, c_i, n1, n2) for th in theta_vals]
    outer_case2 = [superellipsoid_point_3d(phi_forfiv, th, a_o, b_o, c_o, n1, n2) for th in theta_vals]
    part_name = "Layer_" + str(layer_index + 1)
    p = model.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    theta_data = []
    for theta_deg in theta_sections:
        print("Creating section at θ =", theta_deg)
        theta = math.radians(theta_deg)
        # Inner and outer curves for this θ
        inner_curve = [superellipsoid_point_3d(ph, theta, a_i, b_i, c_i, n1, n2) for ph in phi_vals]
        outer_curve = [superellipsoid_point_3d(ph, theta, a_o, b_o, c_o, n1, n2) for ph in phi_vals]
        # Apply backward shift to last point
        inner_curve[-1] = (0.0, 0.0, c_i)
        outer_curve[-1] = (0.0, 0.0, c_o)
        print(inner_curve[-1])
        print(inner_curve[-2])
        print(inner_curve[-3])
        # Create datums and wires
        create_datums(p, inner_curve + outer_curve)
        create_wire_splines(p, inner_curve)
        create_wire_splines(p, outer_curve)
        e1 = p.edges
        inner_ref = inner_curve[-2]
        outer_ref = outer_curve[-2]
        # Create datum points on the edges
        dp_inner = p.DatumPointByEdgeParam(edge=e1.findAt(coordinates=inner_ref), parameter=0.0001)
        dp_outer = p.DatumPointByEdgeParam(edge=e1.findAt(coordinates=outer_ref), parameter=0.0001)
        dp_mid = p.DatumPointByMidPoint(point1=p.datums[dp_inner.id], point2=p.datums[dp_outer.id])
        # Create spline through these datum points
        p.WireSpline(points=(p.datums[dp_inner.id], p.datums[dp_mid.id], p.datums[dp_outer.id]))
        # --- Create a wire spline through these three points (inner → mid → outer)
        dp_inner_start = p.DatumPointByCoordinate(coords=inner_curve[0])
        dp_outer_start = p.DatumPointByCoordinate(coords=outer_curve[0])
        dp_mid_start = p.DatumPointByMidPoint(point1=p.datums[dp_inner_start.id], point2=p.datums[dp_outer_start.id])
        p.WireSpline(points=(p.datums[dp_inner_start.id], p.datums[dp_mid_start.id], p.datums[dp_outer_start.id]))
        theta_data.append({
            "theta_deg": theta_deg,
            "inner_quat": inner_curve[len(inner_curve)//4],
            "inner_threequat": inner_curve[len(inner_curve)*3//4],
            "outer_quat": outer_curve[len(outer_curve)//4],
            "outer_threequat": outer_curve[len(outer_curve)*3//4],
            "dp_mid_coord": p.datums[dp_mid.id].pointOn,
            "dp_mid_start_coord": p.datums[dp_mid_start.id].pointOn
        })
    # --- Case 1 (phi=0°) wires ---
    create_datums(p, inner_case1 + outer_case1)
    create_wire_splines(p, inner_case1)
    create_wire_splines(p, outer_case1)
    # --- Case 2 (phi=45°) wires ---
    create_datums(p, inner_case2 + outer_case2)
    create_wire_splines(p, inner_case2)
    create_wire_splines(p, outer_case2)
    # --- Remove extra wires ---
    center = (0.0, 0.0, c_o)
    radius = 0.5
    edges_near_top = p.edges.getByBoundingSphere(center=center, radius=radius)
    center_bottom = (0.0, 0.0, c_i)
    edges_near_bottom = p.edges.getByBoundingSphere(center=center_bottom, radius=radius)
    RemoveWireEdges = list(edges_near_top) + list(edges_near_bottom)
    if RemoveWireEdges:
        p.RemoveWireEdges(wireEdgeList=RemoveWireEdges)
    num_theta_sections = len(theta_sections)
    inner_a, inner_b, inner_c = a_i, b_i, c_i
    outer_a, outer_b, outer_c = a_o, b_o, c_o
    # --- Create loft sections (same as before) ---
    e1 = p.edges
    loftsections = []
    for section in theta_data:
        pts = (
            e1.findAt(coordinates=section["inner_quat"]),
            e1.findAt(coordinates=section["dp_mid_coord"]),
            e1.findAt(coordinates=section["inner_threequat"]),
            e1.findAt(coordinates=section["outer_quat"]),
            e1.findAt(coordinates=section["outer_threequat"]),
            e1.findAt(coordinates=section["dp_mid_start_coord"]),
        )
        loftsections.append(pts)
    theta_mid = [
        (theta_sections[i] + theta_sections[i + 1]) / 2.0
        for i in range(len(theta_sections) - 1)
    ]
    search_tol = 0.1  # desired tolerance
    inner_path_edges_zero = []
    outer_path_edges_zero = []
    for t_deg in theta_mid:
        theta_rad = math.radians(t_deg)
        # Compute points
        phi = 0.0
        inner_pt = superellipsoid_point_3d(phi, theta_rad, inner_a, inner_b, inner_c, n1, n2)
        outer_pt = superellipsoid_point_3d(phi, theta_rad, outer_a, outer_b, outer_c, n1, n2)
        # Get closest edge for inner point
        found_inner = e1.getClosest(coordinates=(tuple(inner_pt),), searchTolerance=search_tol)
        if found_inner:
            closest_inner = found_inner[0][1]
            edge_inner = e1.findAt((closest_inner,))
            if edge_inner:
                inner_path_edges_zero.append(edge_inner[0])
        # Get closest edge for outer point
        found_outer = e1.getClosest(coordinates=(tuple(outer_pt),), searchTolerance=search_tol)
        if found_outer:
            closest_outer = found_outer[0][1]
            edge_outer = e1.findAt((closest_outer,))
            if edge_outer:
                outer_path_edges_zero.append(edge_outer[0])
    # Convert to tuple for SolidLoft
    inner_path_edges_zero = tuple(inner_path_edges_zero)
    outer_path_edges_zero = tuple(outer_path_edges_zero)
    inner_path_edges_fourfive = []
    outer_path_edges_fourfive = []
    for t_deg in theta_mid:
        theta_rad = math.radians(t_deg)
        # Compute points
        phi = math.radians(45)
        inner_pt = superellipsoid_point_3d(phi, theta_rad, inner_a, inner_b, inner_c, n1, n2)
        outer_pt = superellipsoid_point_3d(phi, theta_rad, outer_a, outer_b, outer_c, n1, n2)
        # Get closest edge for inner point
        found_inner = e1.getClosest(coordinates=(tuple(inner_pt),), searchTolerance=search_tol)
        if found_inner:
            closest_inner = found_inner[0][1]
            edge_inner = e1.findAt((closest_inner,))
            if edge_inner:
                inner_path_edges_fourfive.append(edge_inner[0])
        # Get closest edge for outer point
        found_outer = e1.getClosest(coordinates=(tuple(outer_pt),), searchTolerance=search_tol)
        if found_outer:
            closest_outer = found_outer[0][1]
            edge_outer = e1.findAt((closest_outer,))
            if edge_outer:
                outer_path_edges_fourfive.append(edge_outer[0])
    # Convert to tuple for SolidLoft
    inner_path_edges_fourfive = tuple(inner_path_edges_fourfive)
    outer_path_edges_fourfive = tuple(outer_path_edges_fourfive)
    # --- Define SolidLoft ---
    p.SolidLoft(
        loftsections=loftsections,
        paths=(outer_path_edges_zero, inner_path_edges_zero, inner_path_edges_fourfive, outer_path_edges_fourfive),
        globalSmoothing=ON
    )
    return p

model = mdb.models['SuperEllipse']
for i in range(num_layers):
    create_layer_part(model, i, num_theta_sections)

############ Create SuperEllipsoid by merging layers ############

def assemble_and_merge_layers(num_layers):
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    instance_names = []
    for i in range(1, num_layers + 1):
        part_name = "Layer_" + str(i)
        inst_name = "Layer_" + str(i) + "-1"
        part = model.parts[part_name]
        assembly.Instance(name=inst_name, part=part, dependent=ON)
        instance_names.append(inst_name)
    if num_layers > 1:
        instances_to_merge = tuple(assembly.instances[name] for name in instance_names)
        assembly.InstanceFromBooleanMerge(
            name='SuperEllipsoid',
            instances=instances_to_merge,
            keepIntersections=ON,
            originalInstances=DELETE,
            domain=GEOMETRY
        )
    else:
        mdb.models['SuperEllipse'].parts.changeKey(fromName='Layer_1',toName='SuperEllipsoid')
        mdb.models['SuperEllipse'].rootAssembly.features.changeKey(fromName='Layer_1-1', toName='SuperEllipsoid-1')

assemble_and_merge_layers(num_layers)

############ Material properties ############
model.Material(name='Aluminium')
model.materials['Aluminium'].Elastic(table=((70000.0,0.34), ))

############ Section assignment ############
model.HomogeneousSolidSection(name='AL_section', 
    material='Aluminium', thickness=None)
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']

current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
)
region = regionToolset.Region(cells=current_cells)
p.SectionAssignment(region=region, sectionName='AL_section', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
mdb.models['SuperEllipse'].parts['SuperEllipsoid'].setValues(
    geometryRefinement=FINE)

############ Partitioning ############

#### Strat - 4 ####
def superellipsoid_normal(phi, theta, a, b, c, n1, n2):
    """Analytical unit normal vector on superellipsoid surface."""
    cphi, sphi = math.cos(phi), math.sin(phi)
    ctheta, stheta = math.cos(theta), math.sin(theta)
    # Use your signed_power(base, exp)
    nx = (1.0 / a) * signed_power(cphi, 2 - n1) * signed_power(ctheta, 2 - n2)
    ny = (1.0 / b) * signed_power(cphi, 2 - n1) * signed_power(stheta, 2 - n2)
    nz = (1.0 / c) * signed_power(sphi, 2 - n1)
    # Normalize
    norm = math.sqrt(nx**2 + ny**2 + nz**2)
    return (nx / norm, ny / norm, nz / norm)

def offset_point_along_normal(point, normal, t):
    return tuple(point[i] + t*normal[i] for i in range(3))

theta_case = math.radians(0)
phi_vals = [math.radians(15), math.radians(75)]

for phi_example in phi_vals:
    # Compute inner and outer points on superellipsoid surface
    pt_inner = superellipsoid_point_3d(phi_example, theta_case, a, b, c, n1, n2)
    n_vec = superellipsoid_normal(phi_example, theta_case, a, b, c, n1, n2)
    pt_inner_offset = offset_point_along_normal(pt_inner, n_vec, 1e-3)
    # Create datum points and plane for partition
    datum_inner = p.DatumPointByCoordinate(coords=pt_inner)
    datum_outer = p.DatumPointByCoordinate(coords=pt_inner_offset)
    datum_xy = p.DatumPointByCoordinate(coords=(pt_inner[1], pt_inner[0], pt_inner[2]))
    plane_datum = p.DatumPlaneByThreePoints(
        point1=p.datums[datum_inner.id],
        point2=p.datums[datum_outer.id],
        point3=p.datums[datum_xy.id]
    )  
    cell_to_partition = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )  
    p.PartitionCellByDatumPlane(
        cells=cell_to_partition,
        datumPlane=p.datums[plane_datum.id]
    )

# ---- Partition for face ----#

num_partitions = 6  # number of φ partitions (between 15° and 75°)
phi_min = math.radians(15)
phi_max = math.radians(75)
phi_step = (phi_max - phi_min) / (num_partitions + 1)
phi_face = [phi_min + i * phi_step for i in range(1, num_partitions + 1)]
print("φ partitions (deg):", [math.degrees(ph) for ph in phi_face])  # avoid shadowing

# --- Access model and geometry ---
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']  # part
pf = p.faces      
layer_thickness = float(thick / num_layers)
theta_spl = [math.radians(0), math.radians(90)]
# --- Loop through φ partitions ---
for theta_case in theta_spl:
    for phi_example in phi_face:
        # Create one datum plane per φ
        a_i = a  # use base geometry for datum construction
        b_i = b
        c_i = c
        # Reference point and its normal
        pt_inner = superellipsoid_point_3d(phi_example, theta_case, a_i, b_i, c_i, n1, n2)
        n_vec = superellipsoid_normal(phi_example, theta_case, a_i, b_i, c_i, n1, n2)
        pt_inner_offset = offset_point_along_normal(pt_inner, n_vec, 1e-3)
        pt_xy = (pt_inner[1], pt_inner[0], pt_inner[2])
        # Create datum plane once
        datum_inner = p.DatumPointByCoordinate(coords=pt_inner)
        datum_outer = p.DatumPointByCoordinate(coords=pt_inner_offset)
        datum_xy = p.DatumPointByCoordinate(coords=pt_xy)
        plane_datum = p.DatumPlaneByThreePoints(
            point1=p.datums[datum_inner.id],
            point2=p.datums[datum_outer.id],
            point3=p.datums[datum_xy.id]
        )
        # --- Collect all faces across layers for this φ ---
        all_faces = []
        for layer_index in range(1, num_layers + 1):
            frac = float(layer_thickness / 2)
            a_i = a + frac + (layer_index - 1) * layer_thickness
            b_i = b + frac + (layer_index - 1) * layer_thickness
            c_i = c + frac + (layer_index - 1) * layer_thickness
            # find midpoint between partitions (use φ between faces)
            phi_mid = phi_example - (phi_step / 2.0)
            if phi_mid < phi_min:
                phi_mid = phi_min + (phi_step / 2.0)
            face_pt = superellipsoid_point_3d(phi_mid, theta_case, a_i, b_i, c_i, n1, n2)
            try:
                face_ref = pf.findAt((tuple(face_pt),))  # <-- use FaceArray.findAt
                if face_ref:
                    all_faces.append(face_ref[0])       # append the Face object
            except:
                pass
        # --- Partition all faces at once ---
        if all_faces:
            print("Partitioning φ =", round(math.degrees(phi_example), 2),
                "with", len(all_faces), "faces.")
            p.PartitionFaceByDatumPlane(
                faces=tuple(all_faces),                 # tuple of Face objects
                datumPlane=p.datums[plane_datum.id]
            )
        else:
            print("No faces found for φ =", round(math.degrees(phi_example), 2))


############ Step ############
model.StaticStep(name='LoadingStep', previous='Initial', 
    initialInc=1, minInc=1e-05, maxInc=1.0)

############ Create Surface ############
### Strat - 1 ###
phi_surf = []

for i in range(0, len(phi_vals)):
    phi_surf.append(phi_vals[i] - math.radians(10))

ph_o = math.radians(90)
if n2 and n1 >= 1.0:
    x = math.radians(1)
else:
    x = math.radians(0.01)

phi_surf.append(ph_o - x)
phi_surf_deg = [math.degrees(t) for t in phi_surf]
print("Partition angles (degrees): ", phi_surf_deg)

a1 = mdb.models['SuperEllipse'].rootAssembly
s1 = a1.instances['SuperEllipsoid-1'].faces
search_tol = 0.1
faces_combined = []

for i, phi_i in enumerate(phi_surf):
    theta_i = math.radians(45)
    pt = superellipsoid_point_3d(phi_i, theta_i, a, b, c, n1, n2)
    a1.DatumPointByCoordinate(coords=pt)
    found = s1.getClosest(coordinates=(tuple(pt),), searchTolerance=search_tol)
    if found:
        face_pt = found[0][1]
        seq_face = s1.findAt((face_pt,),)
        if seq_face:
            faces_combined.append(seq_face)

if faces_combined:
    a1.Surface(side1Faces=faces_combined, name='Surf_Load')

########### Create Set ############

layer_thickness = float(thick / num_layers)
for theta_val, set_suffix in zip([0.0, math.radians(90)], ['y', 'x']):
    layer_face_points = []
    for layer_index in range(1, num_layers + 1):
        frac = float((layer_thickness / 2))
        a_i = a + frac + (layer_index - 1) * layer_thickness
        b_i = b + frac + (layer_index - 1) * layer_thickness
        c_i = c + frac + (layer_index - 1) * layer_thickness
        phi_set = [
            phi_vals[0] / 2,
            (phi_vals[0] + phi_vals[1]) / 2 + phi_vals[0],
            math.radians(90) - phi_vals[0] / 2
        ]
        for phi_mid in phi_set:
            pt = superellipsoid_point_3d(phi_mid, theta_val, a_i, b_i, c_i, n1, n2)
            layer_face_points.append(pt)
    layer_face_objs = []
    for pt in layer_face_points:
        found = s1.findAt((pt,))
        if found:
            layer_face_objs.append(found)
    if layer_face_objs:
        a1.Set(faces=layer_face_objs, name='Set-Layer-' + set_suffix)

layer_face_points_45 = []
theta_45 = math.radians(45)
phi_0 = 0.0
for layer_index in range(1, num_layers + 1):
    frac = float((layer_thickness / 2))
    a_i = a + frac + (layer_index - 1) * layer_thickness
    b_i = b + frac + (layer_index - 1) * layer_thickness
    c_i = c + frac + (layer_index - 1) * layer_thickness
    pt = superellipsoid_point_3d(phi_0, theta_45, a_i, b_i, c_i, n1, n2)
    layer_face_points_45.append(pt)

layer_face_objs_45 = []
for pt in layer_face_points_45:
    found = s1.findAt((pt,))
    if found:
        layer_face_objs_45.append(found)

if layer_face_objs_45:
    a1.Set(faces=layer_face_objs_45, name='Set-Layer-z')

############ Load ############
model.Pressure(name='SurfLoad', createStepName='LoadingStep', 
    region=a1.surfaces['Surf_Load'], magnitude=pressure_value)

############ BC ############
region = a1.sets['Set-Layer-x']
mdb.models['SuperEllipse'].XsymmBC(name='BC-2', createStepName='Initial', 
    region=region, localCsys=None)

region = a1.sets['Set-Layer-y']
mdb.models['SuperEllipse'].YsymmBC(name='BC-1', createStepName='Initial', 
    region=region, localCsys=None)

region = a1.sets['Set-Layer-z']
mdb.models['SuperEllipse'].ZsymmBC(name='BC-3', createStepName='Initial', 
    region=region, localCsys=None)

############ Mesh ############
p1 = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
p1.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
p1.generateMesh()
