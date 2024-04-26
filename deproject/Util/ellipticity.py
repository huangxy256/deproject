import numpy as np

def Ellipticity_from_grid(I_xy, x_grid, y_grid, weight_map = None):
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

def Center_xy(I_xy, x_grid, y_grid, weight_map = None):
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

def Axis_ratio2ellipticity(q):
    """Compute ellipticity of an ellipse from its axis ratio, e = (1-q) / (1+q)

    Args:
        q (float): minor axis / major axis

    Returns:
        float: ellipticity (1-q)/(1+q)
    """
    return (1 - q) / (1 + q)


def Ellipticity2axis_ratio(e):
    """Compute the axis ratio of an ellipse from its axis ratio

    Args:
        e (float): ellipticity = (1-q)/(1 + q)

    Returns:
        float: minor axis / major axis
    """
    return (1 - e) / (1 + e)