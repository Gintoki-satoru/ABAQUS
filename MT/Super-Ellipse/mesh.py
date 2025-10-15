import numpy as np
import trimesh

# --- Superellipsoid tank parameters ---
a, b, c = 50, 50, 30       # outer semi-axes
thickness = 5.0            # wall thickness
n1, n2 = 2.5, 2.5          # shape exponents
n_phi, n_theta = 50, 50    # mesh resolution

# --- Parametric grid (1/8th) ---
phi = np.linspace(0, np.pi/2, n_phi)
theta = np.linspace(0, np.pi/2, n_theta)
phi, theta = np.meshgrid(phi, theta)

def superellipsoid_coords(a, b, c, n1, n2, phi, theta):
    x = a * np.sign(np.cos(phi)) * np.abs(np.cos(phi))**n1 * np.sign(np.cos(theta)) * np.abs(np.cos(theta))**n2
    y = b * np.sign(np.cos(phi)) * np.abs(np.cos(phi))**n1 * np.sign(np.sin(theta)) * np.abs(np.sin(theta))**n2
    z = c * np.sign(np.sin(phi)) * np.abs(np.sin(phi))**n1
    return np.column_stack((x.flatten(), y.flatten(), z.flatten()))

# --- Outer and inner vertices ---
outer_vertices = superellipsoid_coords(a, b, c, n1, n2, phi, theta)
inner_vertices = superellipsoid_coords(a - thickness, b - thickness, c - thickness, n1, n2, phi, theta)

# --- Create faces for outer and inner surfaces ---
faces = []
for i in range(n_phi - 1):
    for j in range(n_theta - 1):
        idx = i * n_theta + j
        offset = outer_vertices.shape[0]

        # outer surface
        faces.append([idx, idx + 1, idx + n_theta])
        faces.append([idx + 1, idx + n_theta + 1, idx + n_theta])

        # inner surface (reversed winding for correct normals)
        faces.append([idx + offset, idx + n_theta + offset, idx + 1 + offset])
        faces.append([idx + 1 + offset, idx + n_theta + offset, idx + n_theta + 1 + offset])

# --- Connect outer and inner surfaces (walls) ---
for i in range(n_phi - 1):
    for j in range(n_theta - 1):
        idx = i * n_theta + j
        offset = outer_vertices.shape[0]

        # vertical walls between outer and inner surface (two triangles per quad)
        faces.append([idx, idx + offset, idx + n_theta])
        faces.append([idx + n_theta, idx + offset, idx + n_theta + offset])
        faces.append([idx + 1, idx + 1 + offset, idx + n_theta + 1])
        faces.append([idx + 1 + offset, idx + n_theta + 1 + offset, idx + n_theta + 1])

faces = np.array(faces, dtype=np.int64)

# --- Combine vertices ---
all_vertices = np.vstack((outer_vertices, inner_vertices))

# --- Create mesh ---
mesh = trimesh.Trimesh(vertices=all_vertices, faces=faces)
mesh.export('superellipsoid_1_8th_tank.stl')
print("✅ Hollow 1/8 superellipsoid tank STL created")
