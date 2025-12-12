import numpy as np
import csv
from math import gamma

# =====================================================
# FUNCTIONS (your provided ones, unchanged)
# =====================================================

def superellipsoid_area(a, b, c, n1, n2, N=5250):
    phi = np.linspace(-np.pi/2, np.pi/2, N)
    theta = np.linspace(-np.pi, np.pi, N)
    PHI, THETA = np.meshgrid(phi, theta)
    def spow(x, p): return np.sign(x) * np.abs(x)**p
    X = a * spow(np.cos(PHI), n1) * spow(np.cos(THETA), n2)
    Y = b * spow(np.cos(PHI), n1) * spow(np.sin(THETA), n2)
    Z = c * spow(np.sin(PHI), n1)
    dX_dphi, dX_dtheta = np.gradient(X, phi, theta)
    dY_dphi, dY_dtheta = np.gradient(Y, phi, theta)
    dZ_dphi, dZ_dtheta = np.gradient(Z, phi, theta)
    r_phi = np.stack([dX_dphi, dY_dphi, dZ_dphi], axis=2)
    r_theta = np.stack([dX_dtheta, dY_dtheta, dZ_dtheta], axis=2)
    cross_prod = np.cross(r_phi, r_theta)
    dA = np.sqrt(np.sum(cross_prod**2, axis=2))
    return np.sum(dA) * (phi[1]-phi[0]) * (theta[1]-theta[0])

def superellipsoid_volume(a, b, c, n1, n2):
    e1 = 2/n1
    e2 = 2/n2
    num = 8 * a * b * c * (gamma(1+1/e1))**2 * gamma(1+1/e2)
    den = gamma(1+2/e1) * gamma(1+(1/e2 + 2/e1))
    return num / den

def shape_factor(A, V_inner, thick, a, b, c):
    S_inf = 3.51 * np.sqrt(A)
    S_0 = A / thick
    ls = 2 * max(a, b, c)
    expr = 1.26 - (2 - np.sqrt(A)/ls) / (9 * np.sqrt(1 - 4.79 * (V_inner**(2/3)) / A))
    n = max(expr, 1.0)
    return (S_0**n + S_inf**n)**(1/n)

def equivalent_heat_coeff(T_air, T_LH2, A_outer, hc, S_ins, k_ins, S_liner, k_liner, S_outer, k_outer, A_outer_liner):
    hc_mm = hc
    R2 = 1 / (A_outer * hc_mm)
    R1 = 1 / (k_ins * S_ins)
    Rliner = 1 / (k_liner * S_liner)
    Router = 1 / (k_outer * S_outer)
    Aeq = np.array([
        [R2,     1,     0,    0],
        [Router,-1,     1,    0],
        [R1,     0,    -1,    1],
        [Rliner, 0,     0,   -1]
    ])
    beq = np.array([T_air, 0, 0, -T_LH2])
    Q, T3, T2, T1 = np.linalg.solve(Aeq, beq)
    Q_total = 1.13 * Q * 1.3
    h_eq = Q_total / ((T_air - T_LH2) * A_outer_liner)
    return Q_total, T3, T2, T1, h_eq

# =====================================================
# WORKFLOW PARAMETERS
# =====================================================

t_ins = 16
t_outer = 2
k_liner = 0.0306
k_outer = 0.0306
k_ins = 3.03e-08
T_air = 300
T_LH2 = 20
hc = 10 / 1e6

input_csv = r"U:\Sachdeva\MT_Nair\ABAQUS\Matlab\superellipsoid_parametric_162.csv"
output_csv = "results_Q_total.csv"

# =====================================================
# READ INPUT CSV
# =====================================================

rows = []
with open(input_csv, newline='') as f:
    reader = csv.reader(f, delimiter=';')
    next(reader)  # skip header
    for r in reader:
        # Skip completely empty rows
        if not r or all(x.strip() == "" for x in r):
            continue
        
        # Take only first 6 columns (pad with '' if row is short)
        vals = (r + [""]*6)[:6]
        
        # Skip if any of the required columns are empty
        if any(v.strip() == "" for v in vals):
            continue

        # Convert to float
        a, b, c, t, n1, n2 = map(float, vals)
        rows.append((a, b, c, t, n1, n2))


# =====================================================
# PROCESS EACH GEOMETRY
# =====================================================

results = []
for a, b, c, thick, n1, n2 in rows:
    A_liner = superellipsoid_area(a, b, c, n1, n2)
    A_ins = superellipsoid_area(a+thick, b+thick, c+thick, n1, n2)
    A_outer = superellipsoid_area(a+thick+t_ins, b+thick+t_ins, c+thick+t_ins, n1, n2)
    A_outer_total = superellipsoid_area(a+thick+t_ins+t_outer, b+thick+t_ins+t_outer, c+thick+t_ins+t_outer, n1, n2)

    V_inner = superellipsoid_volume(a+thick, b+thick, c+thick, n1, n2) - superellipsoid_volume(a, b, c, n1, n2)
    V_ins = superellipsoid_volume(a+thick+t_ins, b+thick+t_ins, c+thick+t_ins, n1, n2) - superellipsoid_volume(a+thick, b+thick, c+thick, n1, n2)
    V_outer = superellipsoid_volume(a+thick+t_ins+t_outer, b+thick+t_ins+t_outer, c+thick+t_ins+t_outer, n1, n2) - superellipsoid_volume(a+thick+t_ins, b+thick+t_ins, c+thick+t_ins, n1, n2)

    S_liner = shape_factor(A_liner, V_inner, thick, a, b, c)
    S_ins = shape_factor(A_ins, V_ins, t_ins, a, b, c)
    S_outer = shape_factor(A_outer, V_outer, t_outer, a, b, c)

    Q_total, T3, T2, T1, h_eq = equivalent_heat_coeff(
        T_air, T_LH2, A_outer_total, hc,
        S_ins, k_ins,
        S_liner, k_liner,
        S_outer, k_outer,
        A_ins
    )

    results.append([a, b, c, thick, n1, n2, Q_total, h_eq, T3, T2, T1])

# =====================================================
# WRITE OUTPUT CSV
# =====================================================

with open(output_csv, "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["a_mm","b_mm","c_mm","t_mm","n1","n2","Q_total","h_eq","T3","T2","T1"])
    writer.writerows(results)

print("Finished. Results written to", output_csv)
