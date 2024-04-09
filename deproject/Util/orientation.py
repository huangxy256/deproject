import numpy as np

def _Inclination(ellipsoid_type, theta, phi):
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


def Inclination(oblate, theta, phi, deg = 1):
    """Calculate inclination angle of axisymmetric ellipsoid

    Args:
        oblate (bool): oblate or prolate
        theta (float): projection angle theta [rad]
        phi (float): projection angle phi [rad]

    Returns:
        float: inclination angle [deg]
    """
    if oblate == 1:
        inc = theta
    else:
        inc = np.arccos(np.cos(phi) * np.sin(theta))
        inc = np.where(inc > np.pi/2, np.pi - inc, inc)
    
    if deg == 1:
        inc = np.degrees(inc)
    else:
        pass
    
    return inc
    


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


def Isotropic_inclination(num, oblate=1, deg=0):
    """Sample inclination angle isotropically

    Args:
        num (_type_): number of inclination angles 
        oblate (bool, optional): oblate or prolate
        deg (int, optional): whether to return the inclination angle in degrees or in radians. Defaults to 0 (radians).

    Returns:
        array of float: array of isotropic inclination angle
    """
    theta, phi = Sphere_random_point(num)
    inc = Inclination(oblate=oblate, theta=theta, phi=phi, deg=deg)
    return inc