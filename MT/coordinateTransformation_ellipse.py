import numpy as np
#a is radius of cylinder
# Coordinate transformation: Cylindrical - Cartesian
def cyl2cart_x(a, theta):
    return(a*np.cos(theta))

def cyl2cart_z(a, theta):
    return(a*np.sin(theta))

# Coordinate transformation: Cartesian - Cylindrical
def cart2cyl_radius(x,z):
	return (np.sqrt(x**2+z**2))

def cart2cyl_theta(x,z):
    try:
        return(np.array([np.arctan2(z,x)+2*np.pi if z<0 else np.arctan2(z,x) for x,z in zip(x,z)]))
    except:
        if z<0.0:
            return(np.arctan2(z, x) + 2*np.pi)
        else:
            return(np.arctan2(z, x))

# --- Elliptical versions below ---
# Coordinate transformation: elliptical - Cartesian
def pol2cart_x(a, phi):
    return a * np.cos(phi)

def pol2cart_y(b, phi):
    return b * np.sin(phi)

# def pol2cart3D_x(a, phi,theta):
#     return a * np.cos(phi)*np.cos(theta)

# def pol2cart3D_y(b, phi):
#     return b * np.sin(phi)

# def pol2cart3D_z(a, phi,theta):
#     return a* np.sin(theta)*np.cos(phi)

# def pol2cart3D_x(a, phi,theta):
#     return a * np.cos(phi)*np.cos(theta)

# def pol2cart3D_y(b, phi):
#     return b * np.sin(phi)

# def pol2cart3D_z(a, phi,theta):
#     return a* np.sin(theta)*np.cos(phi)

def pol2cart3D_x(a, t):
    return a * np.cos(t)

def pol2cart3D_z(a, b, alpha, t):
    # d is the semi-axis in the z direction of the intersection ellipse
    d = 1 / np.sqrt(np.tan(alpha)**2 / b**2 + 1 / a**2)
    return d * np.sin(t)

def pol2cart3D_y(a, b, alpha, t):
    # y = z * tan(alpha)
    d = 1 / np.sqrt(np.tan(alpha)**2 / b**2 + 1 / a**2)
    return d * np.tan(alpha) * np.sin(t)

def intersection_ellipse_x(a, b, alpha, t):
    # lambda as defined above
    lam = np.sin(alpha) / (np.cos(alpha) * np.sqrt(2))
    D = 1 / np.sqrt(1 / a**2 + 2 * lam**2 / b**2)
    return (a * np.cos(t) + D * np.sin(t)) / np.sqrt(2)

def intersection_ellipse_z(a, b, alpha, t):
    lam = np.sin(alpha) / (np.cos(alpha) * np.sqrt(2))
    D = 1 / np.sqrt(1 / a**2 + 2 * lam**2 / b**2)
    return (a * np.cos(t) - D * np.sin(t)) / np.sqrt(2)

def intersection_ellipse_y(a, b, alpha, t):
    lam = np.sin(alpha) / (np.cos(alpha) * np.sqrt(2))
    D = 1 / np.sqrt(1 / a**2 + 2 * lam**2 / b**2)
    return np.sqrt(2) * lam * D * np.sin(t)

def cart2pol_phi(x, y, a, b):
    # Works for scalar or np.ndarray
    x = np.asarray(x)
    y = np.asarray(y)
    return(np.array([np.arctan2(a * y, b * x)+2*np.pi if y<0 else np.arctan2(a * y, b * x) for x,y in zip(x,y)]))

def cart2pol_radius(x, y, a, b):
    x = np.asarray(x)
    y = np.asarray(y)
    return np.sqrt((x / a) ** 2 + (y / b) ** 2)

def pol2cart3D_x_ax(a, alpha, t):
    return a * np.cos(t) * np.cos(alpha)

def pol2cart3D_y_ax(b, alpha):
    return b * np.sin(alpha)

def pol2cart3D_z_ax(a, alpha, t):
    return a * np.sin(t) * np.cos(alpha)

def pol2cart3D_x_intersect(a, b, alpha, t):
    phi_star = np.arctan(b / (a * np.tan(alpha) * np.sin(t)))
    return a * np.sin(phi_star) * np.cos(t)

def pol2cart3D_y_intersect(a, b, alpha, t):
    phi_star = np.arctan(b / (a * np.tan(alpha) * np.sin(t)))
    return b * np.cos(phi_star)

def pol2cart3D_z_intersect(a, b, alpha, t):
    phi_star = np.arctan(b / (a * np.tan(alpha) * np.sin(t)))
    return a * np.sin(phi_star) * np.sin(t)