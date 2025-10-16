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
a, b, c = 150.0, 150.0, 150.0   # inner semi-axes
total_length = c
t = 2.5                         # total thickness
n1, n2 = 20, 20                   # shape exponents
num_points = 30                 # points per curve
num_layers = 1                  # number of layers through thickness
num_partitions = 6              # number of partitions
num_theta_sections = 6          # number of θ sections(min 2): For even number, the number of partitions created will be (num_theta_sections + 1)

a_out, b_out, c_out = a + t, b + t, c + t  # outer semi-axes

################ Helper Functions #####################
def signed_power(base, exp):
    return math.copysign(abs(base)**exp, base)

def superellipsoid_point_3d(phi, theta, a, b, c, n1, n2):
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    x = a * signed_power(cos_phi, 2.0/n1) * signed_power(cos_theta, 2.0/n2)
    y = b * signed_power(cos_phi, 2.0/n1) * signed_power(sin_theta, 2.0/n2)
    z = c * signed_power(sin_phi, 2.0/n1)
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
    if n2 > 4 and num_points > 2:
        # Replace second and second-to-last values
        theta_vals[1] = math.radians(0.1)
        theta_vals[-2] = math.radians(89.1) 
    # --- Case 1 (phi=0°) ---
    inner_case1 = [superellipsoid_point_3d(0.0, th, a_i, b_i, c_i, n1, n2) for th in theta_vals]
    outer_case1 = [superellipsoid_point_3d(0.0, th, a_o, b_o, c_o, n1, n2) for th in theta_vals]
    # --- Determine θ angles for cross-sections ---
    if num_theta_sections == 2:
        theta_sections = [0, 90]
    elif num_theta_sections >= 3:
        # Start with required angles
        theta_sections = [0, 45, 90]
        # Add extra evenly spaced angles if more than 3 sections
        if num_theta_sections > 3:
            step = 90.0 / (num_theta_sections - 1)
            extra_angles = [
                round(i * step, 2)
                for i in range(1, num_theta_sections - 1)
                if abs(i * step - 45) > 1e-3
            ]
            theta_sections = sorted(set(theta_sections + extra_angles))
        # Special adjustment if n2 > 4 and num_theta_sections > 3
        if n2 > 4 and num_theta_sections > 3:
            theta_sections[1] = 0.1                     # second angle
            theta_sections[-2] = 89.1                   # second-to-last angle
    else:
        raise ValueError("num_theta_sections must be >= 2")
    part_name = "Layer_" + str(layer_index + 1)
    p = model.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    theta_data = []
    for theta_deg in theta_sections:
        theta = math.radians(theta_deg)
        # Inner and outer curves for this θ
        inner_curve = [superellipsoid_point_3d(ph, theta, a_i, b_i, c_i, n1, n2) for ph in phi_vals]
        outer_curve = [superellipsoid_point_3d(ph, theta, a_o, b_o, c_o, n1, n2) for ph in phi_vals]
        # Apply backward shift to last point
        inner_curve[-1] = (0.0, 0.0, c_i)
        outer_curve[-1] = (0.0, 0.0, c_o)
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
            "inner_mid": inner_curve[len(inner_curve)//2],
            "outer_mid": outer_curve[len(outer_curve)//2],
            "dp_mid_coord": p.datums[dp_mid.id].pointOn,
            "dp_mid_start_coord": p.datums[dp_mid_start.id].pointOn
        })
    # --- Case 1 (phi=0°) wires ---
    create_datums(p, inner_case1 + outer_case1)
    create_wire_splines(p, inner_case1)
    create_wire_splines(p, outer_case1)
    # --- Remove extra wires ---
    center = (0.0, 0.0, c_o)
    radius = 1
    edges_near_top = p.edges.getByBoundingSphere(center=center, radius=radius)
    center_bottom = (0.0, 0.0, c_i)
    edges_near_bottom = p.edges.getByBoundingSphere(center=center_bottom, radius=radius)
    RemoveWireEdges = list(edges_near_top) + list(edges_near_bottom)
    if RemoveWireEdges:
        p.RemoveWireEdges(wireEdgeList=RemoveWireEdges)
    phi = 0.0
    num_theta_sections = len(theta_sections)
    inner_a, inner_b, inner_c = a_i, b_i, c_i
    outer_a, outer_b, outer_c = a_o, b_o, c_o
    # --- Create loft sections (same as before) ---
    e1 = p.edges
    loftsections = []
    for section in theta_data:
        pts = (
            e1.findAt(coordinates=section["inner_mid"]),
            e1.findAt(coordinates=section["dp_mid_coord"]),
            e1.findAt(coordinates=section["outer_mid"]),
            e1.findAt(coordinates=section["dp_mid_start_coord"]),
        )
        loftsections.append(pts)
    # Generate evenly distributed theta midpoints between each theta section
    theta_mid = [
        (theta_sections[i] + theta_sections[i + 1]) / 2.0
        for i in range(len(theta_sections) - 1)
    ]
    search_tol = 0.1  # or your desired tolerance
    inner_path_edges = []
    outer_path_edges = []
    for t_deg in theta_mid:
        theta_rad = math.radians(t_deg)
        # Compute points
        inner_pt = superellipsoid_point_3d(phi, theta_rad, inner_a, inner_b, inner_c, n1, n2)
        outer_pt = superellipsoid_point_3d(phi, theta_rad, outer_a, outer_b, outer_c, n1, n2)
        # Get closest edge for inner point
        found_inner = e1.getClosest(coordinates=(tuple(inner_pt),), searchTolerance=search_tol)
        if found_inner:
            closest_inner = found_inner[0][1]
            edge_inner = e1.findAt((closest_inner,))
            if edge_inner:
                inner_path_edges.append(edge_inner[0])
        # Get closest edge for outer point
        found_outer = e1.getClosest(coordinates=(tuple(outer_pt),), searchTolerance=search_tol)
        if found_outer:
            closest_outer = found_outer[0][1]
            edge_outer = e1.findAt((closest_outer,))
            if edge_outer:
                outer_path_edges.append(edge_outer[0])
    # Convert to tuple for SolidLoft
    inner_path_edges = tuple(inner_path_edges)
    outer_path_edges = tuple(outer_path_edges)
    # --- Define SolidLoft ---
    p.SolidLoft(
        loftsections=loftsections,
        paths=(outer_path_edges, inner_path_edges),
        globalSmoothing=ON
    )
    return p

model = mdb.models['SuperEllipse']
for i in range(num_layers):
    create_layer_part(model, i, num_theta_sections)

"""p.SolidLoft(loftsections=((e1.findAt(coordinates=(inner_case1[mid_idx])), e1.findAt(coordinates=(inner_case2[mid_idx])), e1.findAt(coordinates=( inner_case3[mid_idx]))), 
                          (e1.findAt(coordinates=(outer_case1[mid_idx])), e1.findAt(coordinates=(outer_case2[mid_idx])), e1.findAt( coordinates=(outer_case3[mid_idx])))), 
                          startCondition=NONE, endCondition=NONE)"""

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
    instances_to_merge = tuple(assembly.instances[name] for name in instance_names)
    assembly.InstanceFromBooleanMerge(
        name='SuperEllipsoid',
        instances=instances_to_merge,
        keepIntersections=ON,
        originalInstances=DELETE,
        domain=GEOMETRY
    )

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

#### Strat - 1 ####
'''
n_long = 4  # number of longitudinal partitions
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
csys = p.DatumCsysByThreePoints(
    name='Datum csys-1',
    coordSysType=CARTESIAN,
    origin=(0.0, 0.0, 0.0),
    point1=(1.0, 0.0, 1.0),
    point2=(0.0, 1.0, 0.0)
)
csys_id = csys.id

dp = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
datumPlane0 = p.datums[dp.id]

angle_increment = 90.0 / n_long
d = p.datums
for i in range(1, n_long):
    angle = i * angle_increment
    dp_rot = p.DatumPlaneByRotation(plane=datumPlane0, axis=d[csys_id].axis3, angle=angle)
    dp_rot_obj = p.datums[dp_rot.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_rot_obj, cells=current_cells)

dp = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
datumPlane0 = p.datums[dp.id]
angle_increment = 90.0 / n_long
d = p.datums
for i in range(1, n_long):
    angle = i * angle_increment
    dp_rot = p.DatumPlaneByRotation(plane=datumPlane0, axis=d[csys_id].axis2, angle=angle)
    dp_rot_obj = p.datums[dp_rot.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_rot_obj, cells=current_cells)'''

#### Strat - 2 ####
'''
p = mdb.models['SuperEllipse'].parts['SuperEllipsoid']
dp_xy = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
datumPlane0 = p.datums[dp_xy.id]

n_partitions = 4
total_length = c
dz = total_length / n_partitions

for i in range(1, n_partitions):
    offset = i * dz
    dp_offset = p.DatumPlaneByOffset(plane=datumPlane0, flip=SIDE1, offset=offset)
    dp_offset_obj = p.datums[dp_offset.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_offset_obj, cells=current_cells)

dp_yz = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
datumPlane1 = p.datums[dp_yz.id]

for i in range(1, n_partitions):
    offset = i * dz
    dp_offset = p.DatumPlaneByOffset(plane=datumPlane1, flip=SIDE1, offset=offset)
    dp_offset_obj = p.datums[dp_offset.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_offset_obj, cells=current_cells)'''

#### Strat - 3 ####
'''
def surface_normal(phi, theta, a, b, c, n1, n2, h=1e-6):
    """Approximate normal via finite differences."""
    p0 = superellipsoid_point_3d(phi, theta, a, b, c, n1, n2)
    p_phi = superellipsoid_point_3d(phi + h, theta, a, b, c, n1, n2)
    p_theta = superellipsoid_point_3d(phi, theta + h, a, b, c, n1, n2)
    t_phi = [p_phi[i] - p0[i] for i in range(3)]
    t_theta = [p_theta[i] - p0[i] for i in range(3)]
    nx = t_phi[1]*t_theta[2] - t_phi[2]*t_theta[1]
    ny = t_phi[2]*t_theta[0] - t_phi[0]*t_theta[2]
    nz = t_phi[0]*t_theta[1] - t_phi[1]*t_theta[0]
    norm = math.sqrt(nx**2 + ny**2 + nz**2)
    return (nx/norm, ny/norm, nz/norm)

def offset_point_along_normal(point, normal, t):
    return tuple(point[i] + t*normal[i] for i in range(3))

theta_case = math.radians(90)
phi_vals = [i * math.pi/2 / num_partitions for i in range(1, num_partitions)]

for phi_example in phi_vals:
    # Compute inner and outer points
    pt_inner = superellipsoid_point_3d(phi_example, theta_case, a, b, c, n1, n2)
    pt_outer = superellipsoid_point_3d(phi_example, theta_case, a_out, b_out, c_out, n1, n2)
    # Compute normal and adjust points slightly
    n_vec = surface_normal(phi_example, theta_case, a, b, c, n1, n2)
    pt_inner_adj = offset_point_along_normal(pt_inner, n_vec, -1e-4)
    pt_outer_adj = offset_point_along_normal(pt_outer, n_vec, +1e-4)
    datum_inner = p.DatumPointByCoordinate(coords=pt_inner_adj)
    datum_outer = p.DatumPointByCoordinate(coords=pt_outer_adj)
    datum_xy = p.DatumPointByCoordinate(coords=(pt_inner_adj[2], pt_inner_adj[1], pt_inner_adj[0]))
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

# Final partition near the tip

n = max(n1, n2)
phi_tip_deg = 57.7 * math.exp(-1.68 * n)
phi_tip = math.radians(phi_tip_deg)

pt_inner = superellipsoid_point_3d(phi_tip, theta_case, a, b, c, n1, n2)
pt_outer = superellipsoid_point_3d(phi_tip, theta_case, a_out, b_out, c_out, n1, n2)

n_vec = surface_normal(phi_tip, theta_case, a, b, c, n1, n2)
pt_inner_adj = offset_point_along_normal(pt_inner, n_vec, -1e-4)
pt_outer_adj = offset_point_along_normal(pt_outer, n_vec, +1e-4)

datum_inner = p.DatumPointByCoordinate(coords=pt_inner_adj)
datum_outer = p.DatumPointByCoordinate(coords=pt_outer_adj)
datum_xy = p.DatumPointByCoordinate(coords=(pt_inner_adj[2], pt_inner_adj[1], pt_inner_adj[0]))

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
'''
#### Strat - 4 ####
def superellipsoid_normal(phi, theta, a, b, c, n1, n2):
    """Analytical unit normal vector on superellipsoid surface."""
    cphi, sphi = math.cos(phi), math.sin(phi)
    ctheta, stheta = math.cos(theta), math.sin(theta)
    # Use your signed_power(base, exp)
    nx = (1.0 / a) * signed_power(cphi, 2 - 2/n1) * signed_power(ctheta, 2 - 2/n2)
    ny = (1.0 / b) * signed_power(cphi, 2 - 2/n1) * signed_power(stheta, 2 - 2/n2)
    nz = (1.0 / c) * signed_power(sphi, 2 - 2/n1)
    # Normalize
    norm = math.sqrt(nx**2 + ny**2 + nz**2)
    return (nx / norm, ny / norm, nz / norm)

def offset_point_along_normal(point, normal, t):
    return tuple(point[i] + t*normal[i] for i in range(3))

theta_case = math.radians(0)
n = max(n1, n2)
phi_tip_deg = 87.605 * math.exp(0.0022474 * n)
phi_tip = math.radians(phi_tip_deg)
phi_vals = [i * math.pi/2 / num_partitions for i in range(1, num_partitions)]
# phi_vals.append(phi_tip)
phi_vals = sorted(phi_vals)

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

############ Step ############
model.StaticStep(name='LoadingStep', previous='Initial', 
    initialInc=1, minInc=1e-05, maxInc=1.0)

############ Create Surface ############

phi_surf = []

for i in range(0, len(phi_vals)):
    phi_surf.append(phi_vals[i] - math.radians(0.5))

ph_o = math.radians(90)
x = (ph_o - phi_vals[-1])/2
phi_surf.append(ph_o - x)

a1 = mdb.models['SuperEllipse'].rootAssembly
s1 = a1.instances['SuperEllipsoid-1'].faces
search_tol = 0.1
faces_combined = []

for i, phi_i in enumerate(phi_surf):
    theta_i = math.radians(45)
    pt = superellipsoid_point_3d(phi_i, theta_i, a, b, c, n1, n2)
    p.DatumPointByCoordinate(coords=pt)
    found = s1.getClosest(coordinates=(tuple(pt),), searchTolerance=search_tol)
    if found:
        face_pt = found[0][1]
        seq_face = s1.findAt((face_pt,),)
        if seq_face:
            # a1.Surface(side1Faces=seq_face, name='Surf_%d' % (i+1))
            faces_combined.append(seq_face)

if faces_combined:
    a1.Surface(side1Faces=faces_combined, name='Surf_Load')


"""phi_i = math.radians(5)
theta_i = math.radians(45)
pt = superellipsoid_point_3d(phi_i, theta_i, a, b, c, n1, n2)
a1.DatumPointByCoordinate(coords=pt)

phi_i = math.radians(45)
theta_i = math.radians(45)
pt = superellipsoid_point_3d(phi_i, theta_i, a, b, c, n1, n2)
a1.DatumPointByCoordinate(coords=pt)

phi_i = 1.1693705988362
theta_i = math.radians(45)
pt = superellipsoid_point_3d(phi_i, theta_i, a, b, c, n1, n2)
p.DatumPointByCoordinate(coords=pt)"""