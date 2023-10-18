from deproject.parameters import *
import numpy as np
from deproject.ellipticity import AxisRatio

import deproject as dp

from scipy.integrate import quad

from scipy.interpolate import interp1d

def Project_2d(x_grid, y_grid, zeta, xi, theta, phi, profile):
    """_summary_

    Args:
        x_grid (_type_): _description_
        y_grid (_type_): _description_
        zeta (_type_): _description_
        xi (_type_): _description_
        theta (_type_): _description_
        phi (_type_): _description_
        profile (_type_): _description_

    Returns:
        _type_: _description_
    """
    cap_A = Cap_A(zeta=zeta, xi=xi, theta=theta, phi=phi)
    cap_B = Cap_B(zeta=zeta, xi=xi, theta=theta, phi=phi)
    cap_C = Cap_C(zeta=zeta, xi=xi, theta=theta, phi=phi)
    small_f = Small_f(zeta=zeta, xi=xi, theta=theta, phi=phi)
    
    as_sq = (cap_A * x_grid**2 + cap_B * x_grid * y_grid + cap_C * y_grid**2) / small_f

    Sigma = np.zeros_like(as_sq)

    for i in range(as_sq.shape[0]):
        for j in range(as_sq.shape[1]):
            Sigma[i, j] = quad(profile.Density_in_project, 0, np.inf, args= (as_sq[i, j]))[0] * 2 / np.sqrt(small_f)

    return Sigma


def RadialProfile(x_major, zeta, xi, theta, phi, profile, interpolate=False):
    """Calculate azimuthally averaged radial profile from a 3d triaxial density profile

    Args:
        x_major (float): array of radial variable along the major axis of the projected ellipse
        zeta (float): b/a
        xi (float): c/a
        theta (float): azimuthal angle of the projection
        phi (float): polar angle of the projection
        profile (object): _description_

    Returns:
        _type_: _description_
    """
    cap_A = Cap_A(zeta=zeta, xi=xi, theta=theta, phi=phi)
    cap_B = Cap_B(zeta=zeta, xi=xi, theta=theta, phi=phi)
    cap_C = Cap_C(zeta=zeta, xi=xi, theta=theta, phi=phi)
    small_f = Small_f(zeta=zeta, xi=xi, theta=theta, phi=phi)

    eigenv_minus = (cap_A + cap_C - np.sqrt((cap_A - cap_C)**2 + cap_B**2)) / 2

    as_sq_diag = eigenv_minus / small_f * x_major**2

    Sigma_rad = np.zeros_like(as_sq_diag)

    for i in range(as_sq_diag.shape[0]):
        Sigma_rad[i] = quad(profile.Density_in_project, 0, np.inf, args= (as_sq_diag[i]))[0] * 2 / np.sqrt(small_f)

    axis_ratio = AxisRatio(zeta=zeta, xi=xi, theta=theta, phi=phi)

    if interpolate:
        print('Lower bound for interpolation: {}' .format(np.sqrt(axis_ratio) * np.min(x_major)))
        print('Upper bound for interpolation: {}' .format(np.sqrt(axis_ratio) * np.max(x_major)))
        return interp1d(np.sqrt(axis_ratio) * x_major, Sigma_rad)
    else: 
        return Sigma_rad, np.sqrt(axis_ratio) * x_major


def RadialProfile_any(x, x_type, zeta, xi, theta, phi, profile, interpolate = False, get_interpolation_range = False):

    cap_A = Cap_A(zeta=zeta, xi=xi, theta=theta, phi=phi)
    cap_B = Cap_B(zeta=zeta, xi=xi, theta=theta, phi=phi)
    cap_C = Cap_C(zeta=zeta, xi=xi, theta=theta, phi=phi)
    small_f = Small_f(zeta=zeta, xi=xi, theta=theta, phi=phi)

    eigenv_minus = (cap_A + cap_C - np.sqrt((cap_A - cap_C)**2 + cap_B**2)) / 2
    eigenv_plus = (cap_A + cap_C + np.sqrt((cap_A - cap_C)**2 + cap_B**2)) / 2

    if x_type == 'major':

        as_sq_diag = eigenv_minus / small_f * x**2

        Sigma_rad = np.zeros_like(as_sq_diag)

        for i in range(as_sq_diag.shape[0]):
            Sigma_rad[i] = quad(profile.Density_in_project, 0, np.inf, args= (as_sq_diag[i]))[0] * 2 / np.sqrt(small_f)

        if interpolate:
            print('Lower bound for interpolation: {}' .format(np.min(x)))
            print('Upper bound for interpolation: {}' .format(np.max(x)))
            if get_interpolation_range == True:
                return interp1d(x, Sigma_rad), {'lower_bound': np.min(x), 'upper_bound': np.max(x)}
            else:
                return interp1d(x, Sigma_rad)
        else:
            return Sigma_rad, x

    elif x_type == 'average':

        axis_ratio = AxisRatio(zeta=zeta, xi=xi, theta=theta, phi=phi)

        as_sq_diag = eigenv_minus / small_f * x**2

        Sigma_rad = np.zeros_like(as_sq_diag)

        for i in range(as_sq_diag.shape[0]):
            Sigma_rad[i] = quad(profile.Density_in_project, 0, np.inf, args= (as_sq_diag[i]))[0] * 2 / np.sqrt(small_f)

        if interpolate:
            print('Lower bound for interpolation: {}' .format(np.sqrt(axis_ratio) * np.min(x)))
            print('Upper bound for interpolation: {}' .format(np.sqrt(axis_ratio) * np.max(x)))
            if get_interpolation_range == True:
                return interp1d(np.sqrt(axis_ratio) * x, Sigma_rad), {'lower_bound': np.sqrt(axis_ratio) * np.min(x), 'upper_bound': np.sqrt(axis_ratio) * np.max(x)}
            else:
                return interp1d(np.sqrt(axis_ratio) * x, Sigma_rad)
        else: 
            return Sigma_rad, np.sqrt(axis_ratio) * x

    elif x_type == 'minor':

        as_sq_diag = eigenv_plus / small_f * x**2

        Sigma_rad = np.zeros_like(as_sq_diag)

        for i in range(as_sq_diag.shape[0]):
            Sigma_rad[i] = quad(profile.Density_in_project, 0, np.inf, args= (as_sq_diag[i]))[0] * 2 / np.sqrt(small_f)

        if interpolate:
            print('Lower bound for interpolation: {}' .format(np.min(x)))
            print('Upper bound for interpolation: {}' .format(np.max(x)))
            if get_interpolation_range == True:
                return interp1d(x, Sigma_rad), {'lower_bound': np.min(x), 'upper_bound': np.max(x)}
            else:
                return interp1d(x, Sigma_rad)
        else:
            return Sigma_rad, x

    else:
        raise ValueError("Supported x_type: 'major', 'average' and 'minor'.")