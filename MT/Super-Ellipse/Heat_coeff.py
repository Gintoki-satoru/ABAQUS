import numpy as np
from math import gamma

def superellipsoid_area(a, b, c, n1, n2, N=5250):
    """
    Computes the surface area of a superellipsoid in the first octant,
    then multiplies by 4 to get full area (because phi,theta ∈ [0,π/2]).
    Returned units: same as input (mm^2 if a,b,c in mm).
    """

    phi = np.linspace(-np.pi/2, np.pi/2, N)
    theta = np.linspace(-np.pi, np.pi, N)
    PHI, THETA = np.meshgrid(phi, theta)

    # signed power
    def spow(x, p):
        return np.sign(x) * np.abs(x)**p

    # parametric surface
    X = a * spow(np.cos(PHI), n1) * spow(np.cos(THETA), n2)
    Y = b * spow(np.cos(PHI), n1) * spow(np.sin(THETA), n2)
    Z = c * spow(np.sin(PHI), n1)

    # derivatives
    dX_dphi, dX_dtheta = np.gradient(X, phi, theta)
    dY_dphi, dY_dtheta = np.gradient(Y, phi, theta)
    dZ_dphi, dZ_dtheta = np.gradient(Z, phi, theta)

    r_phi   = np.stack([dX_dphi,   dY_dphi,   dZ_dphi], axis=2)
    r_theta = np.stack([dX_dtheta, dY_dtheta, dZ_dtheta], axis=2)

    cross_prod = np.cross(r_phi, r_theta)
    dA = np.sqrt(np.sum(cross_prod**2, axis=2))

    dphi   = phi[1] - phi[0]
    dtheta = theta[1] - theta[0]

    return np.sum(dA) * dphi * dtheta

def superellipsoid_volume(a, b, c, n1, n2):
    e1 = 2/n1
    e2 = 2/n2
    num = 8 * a * b * c * (gamma(1+1/e1))**2 * gamma(1+1/e2)
    den = gamma(1+2/e1) * gamma(1+(1/e2 + 2/e1))
    return num / den

def shape_factor(A, V_inner, thick, a, b, c):
    """
    Computes shape factor S as in your MATLAB code.
    A in mm^2, V_inner in mm^3, thick in mm
    """
    S_inf = 3.51 * np.sqrt(A)
    S_0 = A / thick
    ls = 2 * max(a, b, c)

    expr = 1.26 - (2 - np.sqrt(A)/ls) / (9 * np.sqrt(1 - 4.79 * (V_inner**(2/3)) / A))
    n = max(expr, 1.0)

    S = (S_0**n + S_inf**n)**(1/n)
    return S

import numpy as np

def equivalent_heat_coeff(
    T_air, T_LH2,
    A_outer,              # mm^2
    hc,                   # W/m2-K
    S_ins, k_ins,         # insulation shape factor, conductivity (W/mm-K)
    S_liner, k_liner,     # liner shape factor, conductivity
    S_outer, k_outer,     # outer wall shape factor, conductivity
    A_outer_liner         # mm^2 reference area for h_eq
):
    """
    Returns: Q_total, T3, T2, T1, h_eq
    """

    # Outer convection
    hc_mm = hc
    R2 = 1 / (A_outer * hc_mm)

    # Insulation
    R1 = 1 / (k_ins * S_ins)

    # Liner
    Rliner = 1 / (k_liner * S_liner)

    # Outer wall metal
    Router = 1 / (k_outer * S_outer)

    # Build linear system for unknowns: Q, T3, T2, T1
    # Equations:
    # R2*Q + T3                 = T_air
    # Router*Q - T3 + T2        = 0
    # R1*Q     - T2 + T1        = 0
    # Rliner*Q         - T1     = -T_LH2

    Aeq = np.array([
        [ R2,      1,     0,    0],
        [Router,  -1,     1,    0],
        [ R1,      0,    -1,    1],
        [Rliner,   0,     0,   -1]
    ])

    beq = np.array([
        T_air,
        0,
        0,
        -T_LH2
    ])

    x = np.linalg.solve(Aeq, beq)
    Q, T3, T2, T1 = x

    # Apply your LH2 scaling
    Q_total = 1.13 * Q * 1.3

    # Equivalent heat transfer coefficient
    h_eq = Q_total / ((T_air - T_LH2) * A_outer_liner)

    return Q_total, T3, T2, T1, h_eq

# 1. Compute area
# a, b, c = 141.0, 141.0, 705.0
# thick = 2                     
# n1, n2 = 1.0, 1.0
# t_ins = 16
# t_outer = 2
# k_liner = 0.0306
# k_outer = 0.0306
# k_ins = 3.0300e-08
# A_liner = superellipsoid_area(a, b, c, n1, n2)
# A_ins = superellipsoid_area(a + thick, b + thick, c + thick, n1, n2)
# A_outer = superellipsoid_area(a + thick + t_ins, b + thick + t_ins, c + thick + t_ins, n1, n2)#
# A_outer_total = superellipsoid_area(a + thick + t_ins + t_outer, b + thick + t_ins + t_outer, c + thick + t_ins + t_outer, n1, n2)

# # 2. Compute volume
# V_inner = superellipsoid_volume(a + thick, b + thick, c + thick, n1, n2) - superellipsoid_volume(a, b, c, n1, n2)
# V_ins = superellipsoid_volume(a + thick + t_ins, b + thick + t_ins, c + thick + t_ins, n1, n2) - superellipsoid_volume(a + thick, b + thick, c + thick, n1, n2)
# V_outer = superellipsoid_volume(a + thick + t_ins + t_outer, b + thick + t_ins + t_outer, c + thick + t_ins + t_outer, n1, n2) - superellipsoid_volume(a + thick + t_ins, b + thick + t_ins, c + thick + t_ins, n1, n2)

# # 3. Compute shape factors
# S_liner = shape_factor(A_liner, V_inner, thick, a, b, c)
# S_ins   = shape_factor(A_ins,   V_ins, t_ins,   a, b, c)
# S_outer = shape_factor(A_outer, V_outer, t_outer, a, b, c)

# # 4. Compute heat transfer coefficient and temperatures
# Q_total, T3, T2, T1, h_eq = equivalent_heat_coeff(
#     T_air=300, T_LH2=20,
#     A_outer=A_outer_total,
#     hc=10 / 1e6,
#     S_ins=S_ins, k_ins=k_ins,
#     S_liner=S_liner, k_liner=k_liner,
#     S_outer=S_outer, k_outer=k_outer,
#     A_outer_liner=A_ins
# )

# print(Q_total)