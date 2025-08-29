import numpy as np

# Coordinate transformation: Cylindrical - Cartesian
def pol2cart_x(radius, theta):
    return(radius*np.cos(theta))

def pol2cart_y(radius, theta):
    return(radius*np.sin(theta))

# Coordinate transformation: Cartesian - Cylindrical
def cart2pol_radius(x,y):
	return (np.sqrt(x**2+y**2))

def cart2pol_theta(x,y):
    try:
        return(np.array([np.arctan2(y,x)+2*np.pi if y<0 else np.arctan2(y,x) for x,y in zip(x,y)]))
    except:
        if y<0.0:
            return(np.arctan2(y, x) + 2*np.pi)
        else:
            return(np.arctan2(y, x))

