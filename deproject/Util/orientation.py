import numpy as np

def Inclination(ellipsoid_type, theta, phi):
    """Calculate inclination angle of axisymmetric ellipsoid

    Args:
        ellipsoid_type (str): 'oblate' or 'prolate'
        theta (float): projection angle theta [rad]
        phi (float): projection angle phi [rad]

    Returns:
        float: inclination angle [deg]
    """
    if ellipsoid_type == 'oblate':
        return theta * 180 / np.pi 
    elif ellipsoid_type == 'prolate':
        inc = np.arccos(np.cos(phi) * np.sin(theta))
        inc = inc * 180 / np.pi
        inc = np.where(inc > 90., 180 - inc, inc)
        return inc
    else:
        raise ValueError("Supported ellipsoid_type: ['oblate', 'prolate'].")


def Sphere_random_point(num):
    """Generate random orientations on a sphere, characterized by polar angle theta and azimuthal angle phi

    Args:
        num (int): number of pairs of (theta, phi) to return

    Returns:
        _type_: (theta, phi)
    """
    u, v = np.random.uniform(low=0, high=1, size = (num, 2)).T 
    theta = np.arccos(u)
    phi = np.pi * v
    return theta, phi