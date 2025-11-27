import numpy as np
import pandas as pd
from scipy.optimize import root_scalar
from math import gamma


# --- SUPERELLIPSOID VOLUME FUNCTION ---
def superellipsoid_volume(a, b, c, n1, n2):
    e1 = 2.0 / n1
    e2 = 2.0 / n2

    V = 8.0 * a * b * c * \
        (gamma(1 + 1/e1)**2 * gamma(1 + 1/e2)) / \
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)))
    return V


# --- MAIN FUNCTION ---
def generate_superellipsoid_abc_n_csv(filename="superellipsoid_grid.csv"):

    # ---- INPUTS ----
    n_vals      = [1.0, 0.8, 0.6, 0.4]
    ratio_vals  = [1, 3, 5]

    V_inner_target = 0.0587e9
    V_wall_target  = 0.0017496e9
    V_outer_target = V_inner_target + V_wall_target

    # ---- Generate unique (a,b) pairs enforcing a <= b ----
    ab_pairs = []
    for a in ratio_vals:
        for b in ratio_vals:
            if a <= b:               # remove duplicate permutations
                ab_pairs.append((a, b))

    c_vals = ratio_vals
    total_rows = len(n_vals) * len(ab_pairs) * len(c_vals)

    data = []

    # ---- MAIN LOOP ----
    for n in n_vals:
        for (a_ratio, b_ratio) in ab_pairs:
            for c_ratio in c_vals:

                a0 = a_ratio
                b0 = b_ratio
                c0 = c_ratio

                # --- Step 1: scale to match inner volume ---
                V0 = superellipsoid_volume(a0, b0, c0, n, n)
                s = (V_inner_target / V0) ** (1.0 / 3.0)

                a_in = s * a0
                b_in = s * b0
                c_in = s * c0

                # --- Step 2: solve thickness t using root finding ---
                def f(t):
                    return superellipsoid_volume(a_in + t, b_in + t, c_in + t, n, n) - V_outer_target

                sol = root_scalar(f, bracket=[0, 5], method="bisect")
                t = sol.root

                V_out = superellipsoid_volume(a_in + t, b_in + t, c_in + t, n, n)

                # --- Apply rounding ---
                a_in_r = round(a_in, 2)
                b_in_r = round(b_in, 2)
                c_in_r = round(c_in, 2)
                t_r    = round(t, 3)
                V_out_r = round(V_out, 2)

                # Store row
                data.append([
                    n, a_ratio, b_ratio, c_ratio,
                    a_in_r, b_in_r, c_in_r, t_r, V_out_r
                ])

    # ---- Create DataFrame ----
    df = pd.DataFrame(data, columns=[
        "n", "ratio_a", "ratio_b", "ratio_c",
        "a", "b", "c", "t", "V_outer"
    ])

    # ---- Save CSV ----
    df.to_csv(filename, index=False, sep=";")

    print(f"Saved {len(df)} rows to: {filename}")
    return df



# ---- RUN THE FUNCTION ----
if __name__ == "__main__":
    df = generate_superellipsoid_abc_n_csv("superellipsoid_grid.csv")
