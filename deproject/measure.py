import numpy as np

from scipy.interpolate import interp1d

def Radial2Image(radial_profile, radius, radius_type, axis_ratio, grid_size=None, grid_spacing=100):
    xy_upbound = (radius[-1] - radius[0]) / np.sqrt(2)
    if grid_size == None:
        x = np.linspace(- xy_upbound, xy_upbound, num=grid_spacing)
        y = np.linspace(- xy_upbound, xy_upbound, num=grid_spacing)
    elif grid_size > xy_upbound:
        print('Radius outside interpolation bound!')
        return 0
    elif grid_size < 0:
        print('grid_size has to be positive number!')
        return 0
    else:
        x = np.linspace(- grid_size, grid_size, num=grid_spacing)
        y = np.linspace(- grid_size, gird_size, num=grid_spacing)

    xx, yy = np.meshgrid(x, y)
    sigma = np.zeros_like(xx)
    radius_xy = np.sqrt(xx**2 + yy**2)

    if radius_type == 'log':
        radial_func = interp1d(np.log10(radius * np.sqrt(axis_ratio)), radial_profile, kind='linear', fill_value='extrapolate')
        sigma = radial_func(np.log10(radius_xy))
    elif radius_type == 'linear':
        radial_func = interp1d(radius * np.sqrt(axis_ratio), radial_profile, kind='linear', fill_value='extrapolate')
        sigma = radial_func(radius_xy)
    else:
        print('Unknown radius type!')

    return sigma, xx, yy