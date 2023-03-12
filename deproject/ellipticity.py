from deproject.parameters import *

# print(Cap_A(1, 1, 0, 0))

def Ellipticity(zeta, xi, theta, phi):
    """Calculate ellipticity from axis ratios and angles of projection

    Args:
        zeta (float): Aixs ratio b/a, 0<=zeta<=1
        xi (float): Axis ratio c/a
        theta (float): Poalr angle in spherical coordinate
        phi (float): Azimuthal angle in spherical coordinate

    Returns:
        float: Ellipticity
    """
    A = Cap_A(zeta, xi, theta, phi)
    B = Cap_B(zeta, xi, theta, phi)
    C = Cap_C(zeta, xi, theta, phi)
    
    numerator = A + C - np.sqrt((A - C)**2 + B**2)
    denominator = A + C + np.sqrt((A - C)**2 + B**2)
    e = 1 - np.sqrt(numerator / denominator)
    
    return e

# print(Ellipticity(1, 1, 0, 0))