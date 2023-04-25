import numpy as np

def Cap_A(zeta, xi, theta, phi):
    """Calculate parameter A from Binney 1985

    Args:
        zeta (float): Aixs ratio b/a, 0<=zeta<=1
        xi (float): Axis ratio c/a
        theta (float): Poalr angle in spherical coordinate, in Radian
        phi (float): Azimuthal angle in spherical coordinate, in Radian

    Returns:
        float: Parameter A from Binney 1985
    """
    return (np.cos(2 * phi) * (2 - zeta**2 - xi**2) + (2 + zeta**2 + xi**2) + np.sin(phi)**2 * np.cos(2 * theta) * 2 * (zeta**2 - xi**2)) / (4 * zeta**2 * xi**2)

def Cap_B(zeta, xi, theta, phi):
    """Calculate parameter B from Binney 1985

    Args:
        zeta (float): Aixs ratio b/a, 0<=zeta<=1
        xi (float): Axis ratio c/a
        theta (float): Poalr angle in spherical coordinate, in Radian
        phi (float): Azimuthal angle in spherical coordinate, in Radian

    Returns:
        float: Parameter B from Binney 1985
    """
    return -1 * (np.sin(2 * phi) * (zeta**2 + xi**2 - 2 + np.cos(2 * theta) * (zeta**2 - xi**2))) / (2 * zeta**2 * xi**2)

def Cap_C(zeta, xi, theta, phi):
    """Calculate parameter C from Binney 1985

    Args:
        zeta (float): Aixs ratio b/a, 0<=zeta<=1
        xi (float): Axis ratio c/a
        theta (float): Poalr angle in spherical coordinate, in Radian
        phi (flsoat): Azimuthal angle in spherical coordinate, in Radian

    Returns:
        float: Parameter C from Binney 1985
    """
    return (np.cos(2 * phi) * (zeta**2 + xi**2 - 2) + (2 + zeta**2 + xi**2) + np.cos(phi)**2 * np.cos(2 * theta) * 2 * (zeta**2 - xi**2)) / (4 * zeta**2 * xi**2)

def Small_f(zeta, xi, theta, phi):
    """Calculate parameter f from Binney 1985

    Args:
        zeta (float): Aixs ratio b/a, 0<=zeta<=1
        xi (float): Axis ratio c/a
        theta (float): Poalr angle in spherical coordinate, in Radian
        phi (float): Azimuthal angle in spherical coordinate, in Radian

    Returns:
        float: Paramater f from Binney 1985
    """
    return np.sin(theta)**2 / zeta**2 + np.cos(theta)**2 / xi**2

