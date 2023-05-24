from deproject.parameters import *
import numpy as np

import deproject as dp

from scipy.integrate import quad

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


def RadialProfile(x_major, zeta, xi, theta, phi, profile):
    """_summary_

    Args:
        x_major (_type_): _description_
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

    eigenv_minus = (cap_A + cap_C - np.sqrt((cap_A - cap_C)**2 + cap_B**2)) / 2

    as_sq_diag = eigenv_minus / small_f * x_major**2

    Sigma_rad = np.zeros_like(as_sq_diag)

    for i in range(as_sq_diag.shape[0]):
        Sigma_rad[i] = quad(profile.Density_in_project, 0, np.inf, args= (as_sq_diag[i]))[0] * 2 / np.sqrt(small_f)

    return Sigma_rad