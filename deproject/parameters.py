import numpy as np

__all__ = ['Cap_A', 'Cap_B', 'Cap_C', 'Small_f']

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
    norm_factor = _Norm(zeta=zeta, xi=xi)
    return (np.sin(theta)**2 / zeta**2 + np.cos(theta)**2 / xi**2 * (np.sin(phi)**2 + np.cos(phi)**2 / zeta**2)) * norm_factor**2

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
    norm_factor = _Norm(zeta=zeta, xi=xi)
    return (np.cos(theta) * np.sin(2 * phi) * (1 - 1/zeta**2) / xi**2) * norm_factor**2

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
    norm_factor = _Norm(zeta=zeta, xi=xi)
    return ((np.cos(phi)**2 + np.sin(phi)**2 / zeta**2) / xi**2) * norm_factor**2

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
    norm_factor = _Norm(zeta=zeta, xi=xi)
    return ((np.cos(theta)**2 / xi**2) + np.sin(theta)**2 * (np.cos(phi)**2 + (np.sin(phi)**2 / zeta**2))) * norm_factor

def _Norm(zeta, xi):
    return (zeta * xi)**(2.0/3.0)