from deproject.parameters import *
import numpy as np


def AxisRatio(zeta, xi, theta, phi):
    """Calculate apparent ratio of the projected ellipse from 3d axis ratios and angles of projection

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
    q = np.sqrt(numerator / denominator)

    return q

def Ellipticity(zeta, xi, theta, phi):
    """Calculate elliticity defined via $e = \frac{1-q}{1+q}$

    Args:
        zeta (float): Aixs ratio b/a, 0<=zeta<=1
        xi (float): Axis ratio c/a
        theta (float): Poalr angle in spherical coordinate
        phi (float): Azimuthal angle in spherical coordinate

    Returns:
        _type_: _description_
    """
    q = AxisRatio(zeta = zeta, xi = xi, theta = theta, phi = phi)
    e = (1 - q) / (1 + q)

    return e

def Ellipticity_12(I_xy, x_grid, y_grid, weight_map = None):
    """_summary_

    Args:
        I_xy (_type_): _description_
        x_grid (_type_): _description_
        y_grid (_type_): _description_
        weight_map (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    if weight_map == None:
        weight_map = np.ones_like(I_xy)

    sum_I_xy = np.sum(I_xy * weight_map)
    center_x = np.sum(I_xy * x_grid * weight_map) / sum_I_xy
    center_y = np.sum(I_xy * y_grid * weight_map) / sum_I_xy

    q11 = np.sum(I_xy * (x_grid - center_x)**2 * weight_map) / sum_I_xy
    q22 = np.sum(I_xy * (y_grid - center_y)**2 * weight_map) / sum_I_xy  
    q12 = np.sum(I_xy * (xx_gridx - center_x)* (y_grid - center_y) * weight_map) / sum_I_xy

    e1 = (q11 - q22) / (q11 + q22)
    e2 = 2 * q12 / (q11 + q22)
    e = np.hypot(e1, e2)

    return e1, e2, e

def center_xy(I_xy, x_grid, y_grid, weight_map = None):
    """_summary_

    Args:
        I_xy (_type_): _description_
        x_grid (_type_): _description_
        y_grid (_type_): _description_
        weight_map (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    if weight_map == None:
        weight_map = np.ones_like(I_xy)
    sum_Ixy = np.sum(I_xy * weight_map)
    center_x = np.sum(I_xy * x_grid * weight_map) / sum_Ixy
    center_y = np.sum(I_xy * y_grid * weight_map) / sum_Ixy
    coord = np.array([center_x, center_y])
    return coord