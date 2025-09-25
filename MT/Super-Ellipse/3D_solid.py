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

# ------------------------
# Parameters
# ------------------------
a_in = 150.0
b_in = 100.0
n = 2.0
thick = 2.5
num_points = 100
n_long = 4  # number of longitudinal partitions
n_lat  = 4  # number of latitudinal partitions

################### Create model ###################
Mdb()
modelName = 'SuperEllipse'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
model = mdb.models[modelName]
model.rootAssembly.clearGeometryCache()
model.rootAssembly.regenerate()

################### Create sketch ###################
s = model.ConstrainedSketch(name='__profile__', sheetSize=2*max(a_in+thick,b_in+thick))

axis_line = s.Line(point1=(0.0, 0.0), point2=(0.0, b_in+thick))
s.setAsConstruction(objectList=(axis_line,))

# --- Inner profile ---
inner_points = []
for i in range(num_points + 1):
    y = b_in * i / num_points
    x = a_in * (1 - (y / b_in) ** n) ** (1.0 / n)
    inner_points.append((x, y))

s.Spline(points=inner_points)

# --- Outer profile ---
outer_points = []
for i in range(num_points + 1):
    y = (b_in + thick) * i / num_points
    x = (a_in + thick) * (1 - (y / (b_in + thick)) ** n) ** (1.0 / n)
    outer_points.append((x, y))

outer_points.reverse()
s.Spline(points=outer_points)

s.Line(point1=inner_points[-1], point2=outer_points[0])
s.Line(point1=outer_points[-1], point2=inner_points[0])

################### Create solid part ###################
part_name = 'Dome'
myPart = model.Part(name=part_name, dimensionality=THREE_D,
                    type=DEFORMABLE_BODY)

myPart.BaseSolidRevolve(sketch=s, angle=90.0, flipRevolveDirection=ON)
s.unsetPrimaryObject()

############ Function for selection in parametric form ############

def superellipse_point(u):
    """Parametric 2D superellipse curve point."""
    cos_u = math.cos(u)
    sin_u = math.sin(u)
    x = a_in * math.copysign(abs(cos_u)**(2.0/n), cos_u)
    y = b_in * math.copysign(abs(sin_u)**(2.0/n), sin_u)
    return x, y

def revolved_surface_point(phi, theta):
    """Returns 3D point on revolved surface."""
    x_curve, y_curve = superellipse_point(phi)
    x = x_curve * math.cos(theta)
    y = y_curve
    z = x_curve * math.sin(theta)
    return (x, y, z)

############ Material properties ############
model.Material(name='Aluminium')
model.materials['Aluminium'].Elastic(table=((70000.0,0.34), ))

############ Section assignment ############
mdb.models['SuperEllipse'].HomogeneousSolidSection(name='AL_section', 
    material='Aluminium', thickness=None)

p = mdb.models['SuperEllipse'].parts['Dome']
long = math.pi / 4      # longitude
lat = -math.pi / 4  # latitude
x, y, z = revolved_surface_point(long, lat)
picked_point = (x, y, z)
c = p.cells
cells = c.findAt(((x,y,z), ))
region = regionToolset.Region(cells=cells)
p = mdb.models['SuperEllipse'].parts['Dome']
p.SectionAssignment(region=region, sectionName='AL_section', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

########### Assembly ##################
a = mdb.models['SuperEllipse'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['SuperEllipse'].parts['Dome']
a.Instance(name='Dome-1', part=p, dependent=OFF)
session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)

########### Create longitudinal partition ##################

'''p = mdb.models['SuperEllipse'].parts['Dome']

dp = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
datumPlane0 = p.datums[dp.id]

axis_line = p.datums[1]

angle_increment = 90.0 / n_long

for i in range(1, n_long):
    angle = i * angle_increment
    dp_rot = p.DatumPlaneByRotation(plane=datumPlane0, axis=axis_line, angle=angle)
    dp_rot_obj = p.datums[dp_rot.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_rot_obj, cells=current_cells)

########### Create latitudinal partition ##################

dp_zx = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
datumPlaneZX = p.datums[dp_zx.id]

dy = b / n_lat

for i in range(1, n_lat):
    offset = i * dy
    dp_offset = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=offset)
    dp_offset_obj = p.datums[dp_offset.id]
    current_cells = p.cells.getByBoundingBox(
        xMin=-1e6, xMax=1e6,
        yMin=-1e6, yMax=1e6,
        zMin=-1e6, zMax=1e6
    )
    p.PartitionCellByDatumPlane(datumPlane=dp_offset_obj, cells=current_cells)'''

'''N_lat = 6  # number of latitude partitions
N_lon = 6 # number of longitude partitions

# ==== SUPERELLIPSE ARC-LENGTH DIVISION ====
def arc_length(u1, u2, steps=100):
    """Compute arc length between u1 and u2 numerically."""
    du = (u2 - u1) / steps
    s = 0.0
    for i in range(steps):
        u = u1 + i * du
        x1, y1 = superellipse_point(u)
        x2, y2 = superellipse_point(u+du)
        ds = math.sqrt((x2-x1)**2 + (y2-y1)**2)
        s += ds
    return s

# Compute full arc from u=0 (top) to u=pi/2 (equator)
u_min, u_max = 0.0, math.pi/2
total_arc = arc_length(u_min, u_max)

# Find u-values such that arc length is equally divided
arc_step = total_arc / (N_lat+1)
u_values = [u_min]
s = 0.0
u = u_min
while len(u_values) < N_lat+1:
    # Increment u until next arc step reached
    u += 0.001
    s = arc_length(u_min, u)
    if s >= arc_step * len(u_values):
        u_values.append(u)

# ==== CREATE DATUM PLANES FOR LATITUDES ====
p = mdb.models['SuperEllipse'].parts['Dome']
datum_planes = []

for u in u_values[1:]:  # skip the first (u=0) which is pole
    x, y = superellipse_point(u)
    dp = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=y)
    datum_planes.append(dp)

# ==== CREATE DATUM PLANES FOR LONGITUDES ====
for i in range(1, N_lon):
    theta = i * (90 / N_lon)
    dp = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    dp = p.DatumPlaneByRotation(plane=p.datums[dp.id], axis=p.datums[1],
                                angle=theta)
    datum_planes.append(dp)

# ==== SELECT ENTIRE CELL AND PARTITION ====

for dp in datum_planes:
    cells = p.cells.getByBoundingBox(-1e3, -1e3, -1e3, 1e3, 1e3, 1e3)
    p.PartitionCellByDatumPlane(datumPlane=p.datums[dp.id], cells=cells)'''

# Abaqus script: equal-arc latitudes with normal-direction datum planes
# Run inside Abaqus/CAE (File -> Run script...)

from abaqus import *
from abaqusConstants import *
import math

# -------------------------
# User parameters (edit)
# -------------------------
modelName = 'SuperEllipse'    # model name
partName  = 'Dome'       # part name (solid dome created earlier)
a_radius  = 100.0             # inner semi-axis x (if you made inner/outer use inner)
b_radius  = 100.0             # inner semi-axis y
thickness = 2.5               # thickness used when creating solid (outer = inner + thickness)
n_exp     = 2.0               # superellipse exponent used for profile (implicit)
N_lat     = 12                 # number of latitudinal partitions (number of bands)
num_profile_samples = 1000    # resolution to compute arc length (increase for more accuracy)

# -------------------------
# Helper math functions
# -------------------------
def superellipse_xy(a, b, n, u):
    """Return (x, y) on superellipse profile for parametric angle u in [0, pi/2]."""
    cu = math.cos(u)
    su = math.sin(u)
    x = a * math.copysign(abs(cu)**(2.0 / n), cu)
    y = b * math.copysign(abs(su)**(2.0 / n), su)
    return x, y

def derivative_superellipse(a, b, n, u, h=1e-6):
    """Numerical derivative dx/du, dy/du at u using central difference."""
    u1 = max(0.0, u - h)
    u2 = min(math.pi / 2.0, u + h)
    x1, y1 = superellipse_xy(a, b, n, u1)
    x2, y2 = superellipse_xy(a, b, n, u2)
    dx = (x2 - x1) / (u2 - u1)
    dy = (y2 - y1) / (u2 - u1)
    return dx, dy

def compute_profile_samples(a, b, n, samples):
    """Return arrays of (u_list, x_list, y_list, s_list) for u in [0, pi/2]."""
    u_list = [i * (math.pi / 2.0) / samples for i in range(samples + 1)]
    x_list = []
    y_list = []
    for u in u_list:
        x, y = superellipse_xy(a, b, n, u)
        x_list.append(x)
        y_list.append(y)
    # cumulative arc length along profile from u=0 to u
    s_list = [0.0]
    for i in range(samples):
        dx = x_list[i+1] - x_list[i]
        dy = y_list[i+1] - y_list[i]
        ds = math.hypot(dx, dy)
        s_list.append(s_list[-1] + ds)
    return u_list, x_list, y_list, s_list

def find_u_at_arc(s_target, u_list, s_list):
    """Interpolate to find u such that arc length s(u) = s_target."""
    # assume s_list is monotonic increasing
    if s_target <= 0.0:
        return u_list[0]
    if s_target >= s_list[-1]:
        return u_list[-1]
    # binary search / linear interpolation
    lo = 0
    hi = len(s_list) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if s_list[mid] < s_target:
            lo = mid
        else:
            hi = mid
    # linear interp between lo and hi
    s_lo = s_list[lo]
    s_hi = s_list[hi]
    u_lo = u_list[lo]
    u_hi = u_list[hi]
    if s_hi == s_lo:
        return u_lo
    t = (s_target - s_lo) / (s_hi - s_lo)
    return u_lo + t * (u_hi - u_lo)

# -------------------------
# Main: compute equal-arc latitudes (including endpoints)
# -------------------------
u_samples, x_samples, y_samples, s_samples = compute_profile_samples(a_radius, b_radius, n_exp, num_profile_samples)
total_arc = s_samples[-1]

# we want N_lat bands -> there are N_lat+1 parallels including top and rim; 
# we'll create datum planes for internal parallels only (exclude pole and rim)
arc_targets = [k * total_arc / (N_lat + 1) for k in range(1, N_lat + 1)]  # internal only

u_lats = [(math.pi / 2.0) * math.sin(math.pi * k / (2 * (N_lat + 1))) 
          for k in range(1, N_lat + 1)]
# compute coordinates on inner profile (and outer if thickness > 0)
lat_points_inner = [superellipse_xy(a_radius, b_radius, n_exp, u) for u in u_lats]
outer_a = a_radius + thickness
outer_b = b_radius + thickness
lat_points_outer = [superellipse_xy(outer_a, outer_b, n_exp, u) for u in u_lats]

# We will use the mid-surface point between inner and outer to define plane position
lat_points_mid = [((xi + xo) / 2.0, (yi + yo) / 2.0) 
                  for (xi, yi), (xo, yo) in zip(lat_points_inner, lat_points_outer)]

# -------------------------
# Work inside Part module
# -------------------------
model = mdb.models[modelName]
p = model.parts[partName]

# safety: delete existing datum planes whose names start with 'LatPlane' to avoid duplicates
# (This deletion is optional and assumes those names were used earlier)
try:
    to_delete = []
    for dId, datum in p.datums.items():
        name = getattr(datum, 'name', '')
        if name and 'LatPlane' in name:
            to_delete.append('Datum plane-%s' % str(dId))
    if to_delete:
        p.deleteFeatures(tuple(to_delete))
except Exception:
    pass

# ensure we have a big bounding box helper function to reselect cells
def select_all_cells(part_obj):
    return part_obj.cells.getByBoundingBox(xMin=-1e9, yMin=-1e9, zMin=-1e9,
                                          xMax=1e9, yMax=1e9, zMax=1e9)

# We will create datum planes by three points for each latitude:
# point_on_mid (xmid, ymid, z=0), point offset in normal direction (xmid + nx, ymid + ny, 0),
# and axis reference point at same y: (0, ymid, 0) to ensure plane spans through axis.
# After creating plane, partition the current cells and reselect cells for next iteration.

created_datums = []
for idx, (u_val, (x_mid, y_mid)) in enumerate(zip(u_lats, lat_points_mid), start=1):
    # compute meridional tangent and normal at u_val using central difference
    dx_du, dy_du = derivative_superellipse(a_radius, b_radius, n_exp, u_val)
    # tangent vector t = (dx_du, dy_du); normal in 2D = (-dy, dx)
    nx2 = -dy_du
    ny2 = dx_du
    norm_len = math.hypot(nx2, ny2)
    if norm_len == 0.0:
        # fallback to radial direction
        nx2, ny2 = x_mid, y_mid
        norm_len = math.hypot(nx2, ny2)
        if norm_len == 0.0:
            nx2, ny2 = 1.0, 0.0
            norm_len = 1.0
    nx2 /= norm_len
    ny2 /= norm_len
    # --- Create datum points ---
    dp1 = p.DatumPointByCoordinate(coords=(x_mid, y_mid, 0.0))                 # mid-surface point
    dp2 = p.DatumPointByCoordinate(coords=(x_mid + nx2, y_mid + ny2, 0.0))     # offset along normal
    dp3 = p.DatumPointByCoordinate(coords=(0, y_mid, -x_mid))                   # axis reference point
    # --- Create datum plane using the three datum points ---
    plane = p.DatumPlaneByThreePoints(point1=p.datums[dp1.id],
                                      point2=p.datums[dp2.id],
                                      point3=p.datums[dp3.id])
    plane_obj = p.datums[plane.id]
    # Name datum for bookkeeping
    try:
        plane_obj.name = 'LatPlane-%d' % idx
    except Exception:
        pass
    created_datums.append(plane_obj)
    # --- Partition all current cells with this datum plane ---
    current_cells = select_all_cells(p)
    p.PartitionCellByDatumPlane(datumPlane=plane_obj, cells=current_cells)


# -------------------------
# Optional: create longitudinal planes (even-angle) if desired
# -------------------------
# create base YZ plane at 0, then rotate about axis (datum line) to get N_lon planes
# you can choose N_lon to match circumferential spacing ~ arc spacing if desired
N_lon = None  # set to an integer to enable; if None, skip longitudinal creation
if isinstance(N_lon, int) and N_lon > 0:
    # create base YZ plane through axis
    dp_yz = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    base_plane = p.datums[dp_yz.id]
    # find axis datum (construction line) — we assume datum 1 is the revolve axis line from sketch
    # you might need to adjust if different in your part
    axis_datum = p.datums[1]
    for i in range(1, N_lon):
        angle = i * (90.0 / N_lon)  # for a 90-degree revolved dome; change if full revolution
        dp_rot = p.DatumPlaneByRotation(plane=base_plane, axis=axis_datum, angle=angle)
        dp_rot_obj = p.datums[dp_rot.id]
        current_cells = select_all_cells(p)
        p.PartitionCellByDatumPlane(datumPlane=dp_rot_obj, cells=current_cells)