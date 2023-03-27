from deproject.parameters import *
import numpy as np


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

def Ellipticity_12(xx, yy, Sigma, weight = None):
    """Calculate ellipticities with second moments

    Args:
        xx (float): 2d x coordinate
        yy (float): 2d y coordinate
        Sigma (float): 2d projected image
        weight(float): 2d weight map in calculating second moments

    Returns:
        float, float, float: Ellipticities e1, e2 and e
    """
    if weight == None:
        weight = np.ones_like(Sigma)
    sum_Sigma = np.sum(Sigma * weight)
    center_x = np.sum(Sigma * xx * weight) / sum_Sigma
    center_y = np.sum(Sigma * yy * weight) / sum_Sigma

    q11 = np.sum(Sigma * (xx - center_x)**2 * weight) / sum_Sigma
    q22 = np.sum(Sigma * (yy - center_y)**2 * weight) / sum_Sigma  
    q12 = np.sum(Sigma * (xx - center_x)* (yy - center_y) * weight) / sum_Sigma

    e1 = (q11 - q22) / (q11 + q22)
    e2 = 2 * q12 / (q11 + q22)
    e = np.hypot(e1, e2)

    return e1, e2, e