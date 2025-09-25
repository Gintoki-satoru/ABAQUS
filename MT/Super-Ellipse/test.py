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


################ Parameters #####################
a, b, c = 150.0, 100.0, 120.0   # semi-axes
t = 2.5                         # thickness
n1, n2 = 4, 4             # shape exponents
num_points = 25                 # resolution along curve

a_out, b_out, c_out = a + t, b + t, c + t

# Get model and create a part to hold datum points
model = mdb.models['Model-1']
p = model.Part(name='DatumSuperEllipsoid', dimensionality=THREE_D, type=DEFORMABLE_BODY)

# ------------------------
# Helper function
# ------------------------
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

# ------------------------
# Case 1: phi = 0, theta = 0 → 90 deg
# ------------------------
theta_vals = [math.radians(i*90/num_points) for i in range(num_points+1)]
phi_case1 = 0.0

inner_points_case1 = []
outer_points_case1 = []
for theta in theta_vals:
    inner_points_case1.append(superellipsoid_point_3d(phi_case1, theta, a, b, c, n1, n2))
    outer_points_case1.append(superellipsoid_point_3d(phi_case1, theta, a_out, b_out, c_out, n1, n2))

# ------------------------
# Case 2: phi = 0 → 90 deg, theta = 0
# ------------------------
phi_vals = [math.radians(i*90/num_points) for i in range(num_points+1)]
theta_case2 = 0.0

inner_points_case2 = []
outer_points_case2 = []
for phi in phi_vals:
    inner_points_case2.append(superellipsoid_point_3d(phi, theta_case2, a, b, c, n1, n2))
    outer_points_case2.append(superellipsoid_point_3d(phi, theta_case2, a_out, b_out, c_out, n1, n2))

# ------------------------
# Case 3: phi = 0 → 90 deg, theta = 90 deg
# ------------------------
theta_case3 = math.radians(90)

inner_points_case3 = []
outer_points_case3 = []
for phi in phi_vals:
    inner_points_case3.append(superellipsoid_point_3d(phi, theta_case3, a, b, c, n1, n2))
    outer_points_case3.append(superellipsoid_point_3d(phi, theta_case3, a_out, b_out, c_out, n1, n2))

# ------------------------
# Create datum points for all three cases
# ------------------------
def create_datums(point_list):
    for pt in point_list:
        p.DatumPointByCoordinate(coords=pt)

# Case 1
create_datums(inner_points_case1 + outer_points_case1)

# Case 2
create_datums(inner_points_case2 + outer_points_case2)

# Case 3
create_datums(inner_points_case3 + outer_points_case3)



datum_inner_case1 = []
for pt in inner_points_case1:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

# ------------------------
# Create a wire spline using these datum points
# ------------------------
wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in outer_points_case1:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

# ------------------------
# Create a wire spline using these datum points
# ------------------------
wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in inner_points_case2:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

# ------------------------
# Create a wire spline using these datum points
# ------------------------
wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in outer_points_case2:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

# ------------------------
# Create a wire spline using these datum points
# ------------------------
wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in inner_points_case3:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

# ------------------------
# Create a wire spline using these datum points
# ------------------------
wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)

datum_inner_case1 = []
for pt in outer_points_case3:
    dp = p.DatumPointByCoordinate(coords=pt)
    datum_inner_case1.append(p.datums[dp.id])

# ------------------------
# Create a wire spline using these datum points
# ------------------------
wire_inner_case1 = p.WireSpline(points=datum_inner_case1,
                                mergeType=IMPRINT,
                                meshable=ON,
                                smoothClosedSpline=ON)
p.regenerate()

dp_inner_case1 = inner_points_case1  # list of (x, y, z)
dp_inner_case2 = inner_points_case2
dp_inner_case3 = inner_points_case3

dp_outer_case1 = outer_points_case1
dp_outer_case2 = outer_points_case2
dp_outer_case3 = outer_points_case3

mid_idx = len(dp_inner_case1) // 2

e1 = p.edges
p.SolidLoft(loftsections=((e1.findAt(coordinates=(dp_outer_case1[mid_idx])), e1.findAt(coordinates=(dp_outer_case2[mid_idx])), e1.findAt(coordinates=( dp_outer_case3[mid_idx]))), 
                          (e1.findAt(coordinates=(dp_inner_case1[mid_idx])), e1.findAt(coordinates=(dp_inner_case2[mid_idx])), e1.findAt( coordinates=(dp_inner_case3[mid_idx])))), 
                          startCondition=NONE, endCondition=NONE)