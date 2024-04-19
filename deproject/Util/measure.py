import numpy as np

from scipy.interpolate import interp1d

def _Radial2Image(radial_profile, radius, radius_type, axis_ratio, grid_size=None, grid_spacing=100):
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

def Radial2Image(radial_profile, R_average, radius_type, grid_limit=None, grid_num=200):
    """Convert radial profile to 2d image

    Args:
        radial_profile (float): azimuthally averaged radial profile 
        R_average (float): radial variable corresponding to the radial profile
        radius_type (string): 'log' or 'linear', radial variable linearly defined in log or linear space
        grid_limit (float, optional): maximum radius from center [arcsec]. Defaults to None. If None, will be max(R_average).
        grid_num (int, optional): total number of grid points in 1 dimension. Defaults to 200.

    Returns:
        float: a 2d image corresponding to the input radial profile.
    """
    xy_upbound = (R_average[-1] - R_average[0]) / np.sqrt(2)
    if grid_limit == None:
        x = np.linspace(- xy_upbound, xy_upbound, num=grid_num // 2)
        y = np.linspace(- xy_upbound, xy_upbound, num=grid_num // 2)
    elif grid_limit > xy_upbound:
        print('Radius outside interpolation bound!')
        return 0
    elif grid_limit < 0:
        print('grid_size has to be positive number!')
        return 0
    else:
        x = np.linspace(- grid_limit, grid_limit, num=grid_num // 2)
        y = np.linspace(- grid_limit, grid_limit, num=grid_num // 2)

    xx, yy = np.meshgrid(x, y)
    sigma = np.zeros_like(xx)
    radius_xy = np.sqrt(xx**2 + yy**2)

    if radius_type == 'log':
        radial_func = interp1d(np.log10(R_average), radial_profile, kind='linear', fill_value='extrapolate')
        sigma = radial_func(np.log10(radius_xy))
    elif radius_type == 'linear':
        radial_func = interp1d(R_average, radial_profile, kind='linear', fill_value='extrapolate')
        sigma = radial_func(radius_xy)
    else:
        print('Available radius type: [log, linear].')

    return sigma, xx, yy

def Slope_at_thetaE(radial_profile, R_average, radius_type, thetaE):
    """Calculate the log-log slope of the radial profile at theta_E

    Args:
        radial_profile (float): azimuthally averaged radial profile
        R_average (float): radial variable corresponding to the radial profile
        radius_type (string): 'log'. 
        thetaE (float): Einstein radius

    Returns:
        float: log-log slope of the radial profile at theta_E.
    """
    if radius_type == 'log':
        derivative_func = _Deivative(np.log10(radial_profile), np.log10(R_average), x_type='linear', y_type='linear')
    return derivative_func(thetaE)


def _Deivative(y, x, y_type='log', x_type='log'):
    """_summary_

    Args:
        y (_type_): _description_
        x (_type_): _description_
        y_type (str, optional): _description_. Defaults to 'log'.
        x_type (str, optional): _description_. Defaults to 'log'.

    Returns:
        _type_: _description_
    """
    derivative = np.zeros_like(y)
    if (y_type == 'log' and x_type == 'log'):
        derivative[0] = np.log10(y[1] / y[0]) / np.log10(x[1] / x[0])
        for i in range(1, len(derivative) - 1):
            derivative[i] = np.log10(y[i+1] / y[i-1]) / np.log10(x[i+1] / x[i-1])
        derivative[-1] = np.log10(y[-1] / y[-2]) / np.log10(x[-1] / x[-2])

    elif (y_type == 'linear' and x_type == 'linear'):
        derivative[0] = (y[1] - y[0]) / (x[1] - x[0])
        for i in range(1, len(derivative) - 1):
            derivative[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
        derivative[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])

    return interp1d(x, derivative, kind='linear', fill_value='extrapolate')


def einstein_radius_from_rp(radial_profile, R_circ, R_min=1e-3, R_max=1e1):
    """Compute the Einstein radius from the circularized radial profile

    Args:
        radial_profile (arr of float): circularized radial profile of mass density
        R_circ (arr of float): circularized radius
        R_min (float, optional): minimum radius of the convergence integrand. Defaults to 1e-3.
        R_max (_type_, optional): maximum radius of the convergence integrand (should be larger than the Einstein radius). Defaults to 1e1.

    Returns:
        float: estimate of the Einstein radius
    """
    # define a finer grid for interpolation
    interp_func = interp1d(np.log10(R_circ), np.log10(radial_profile))
    R_min = np.max([R_min, R_circ.min()])
    R_max = np.min([R_max, R_circ.max()])
    R_finer = np.geomspace(R_min, R_max, 10 * len(R_circ))
    radial_profile = 10**interp_func(np.log10(R_finer))

    # perform integral in log space
    radial_profile_ = (radial_profile[1:] + radial_profile[:-1]) / 2
    R_finer_ = (R_finer[1:] + R_finer[:-1]) / 2
    dlog_r = (np.log10(R_finer[2]) - np.log10(R_finer[1])) * np.log(10)
    # add the mass within the innermost bin and assume it is constant
    kappa_innermost = radial_profile[0] * np.pi * R_finer[0]**2

    # the first part is the logarithmic integrand, the second part the circle integrand
    kappa_slice = radial_profile_ * dlog_r * R_finer_ * (2 * np.pi * R_finer_)
    kappa_slice = np.append(kappa_innermost, kappa_slice)

    kappa_cdf = np.cumsum(kappa_slice)
    # calculate average convergence at radius
    kappa_average = kappa_cdf / (np.pi * R_finer**2)

    # we interpolate as the inverse function and evaluate this function for average kappa = 1
    # (assumes monotonic decline in average convergence)
    inv_interp = interp1d(np.log10(kappa_average), np.log10(R_finer))
    try: 
        theta_E = 10**inv_interp(0)
    except:
        theta_E = np.nan
    return theta_E