from deproject.Profiles.sis_truncated_angular import *  
from deproject.Profiles.hernquist import Hernquist
from deproject.projection import Projection
from mgefit.mge_fit_1d import mge_fit_1d
from deproject.Kinematics import mge_misc
import numpy as np


def get_truncsis_mge(sigma_v, lens_cosmo, projection, shape, rc = None, r_min = 0.001, r_max = 1000, r_num = 100, interpolate_num = 300, get_profile = False, kwargs_mge = {}, plot_mge = False):
    """this function calculate the MGE parameter for a triaxial "truncated" SIS profile projected to asigned LoS, and can also return the radial profile

    Args:
        sigma_v (_type_): _description_
        lens_cosmo (_type_): _description_
        projection (_type_): _description_
        shape (_type_): _description_
        rc (_type_, optional): _description_. Defaults to None.
        r_min (float, optional): _description_. Defaults to 0.001.
        r_max (int, optional): _description_. Defaults to 1000.
        r_num (int, optional): _description_. Defaults to 100.
        interpolate_num (int, optional): _description_. Defaults to 300.
        get_profile (bool, optional): _description_. Defaults to False.
        kwargs_mge (dict, optional): _description_. Defaults to {}.
        plot_mge (bool, optional): _description_. Defaults to False.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    profile_sis = SIS_truncated_angular(sigma_v=sigma_v, rc_arcsec = rc, lens_cosmo=lens_cosmo)

    rad_arcsec = np.logspace(np.log10(r_min), np.log10(r_max), num = r_num)

    if shape == 'oblate':
        radial_avg = projection.RadialProfile(R=rad_arcsec, profile=profile_sis, R_align='major', interpolate=True)
    elif shape == 'prolate':
        radial_avg = projection.RadialProfile(R=rad_arcsec, profile=profile_sis, R_align='minor', interpolate=True)
    else: 
        raise ValueError("3D shape {} not recongnized!".format(shape))

    rad_arcsec = np.logspace(np.log10(r_min), np.log10(r_max), num = interpolate_num)

    radial_avg = radial_avg(rad_arcsec) / 1e6

    ngauss = kwargs_mge.get('ngauss', 20)
    inner_slope = kwargs_mge.get('inner_slope', 3.)
    outer_slope = kwargs_mge.get('outer_slope', 1.)
    quiet_mge = kwargs_mge.get('quiet', True)

    sis_mge = mge_fit_1d(rad_arcsec, radial_avg, ngauss=ngauss, inner_slope=inner_slope, outer_slope=outer_slope, plot = False, quiet = quiet_mge)

    sigma_sis = sis_mge.sol[1]
    amp_sis = sis_mge.sol[0] / (np.sqrt(2 * np.pi) * sigma_sis)

    if plot_mge:
        mge_misc.plot_mge(amp_sis, sigma_sis, rad_arcsec, radial_avg)

    if get_profile:
        return amp_sis, sigma_sis, rad_arcsec, radial_avg
    else:
        return amp_sis, sigma_sis


#######################################################################################


def get_hernquist_mge(m_star, Re_kpc, lens_cosmo, projection, shape, r_min = 0.001, r_max = 500, r_num = 100, interpolate_num = 300, get_profile = False, kwargs_mge = {}, plot_mge = False):

    Rs_star = Re_kpc / 1000 * 0.551
    sigma0, Rs_angle_star = lens_cosmo.hernquist_physical2angle(M = m_star, Rs=Rs_star)
    print({'sigma0': sigma0, 'Rs_angle_star': Rs_angle_star})

    profile_hernquist = Hernquist(Rs=Rs_angle_star, sigma0=sigma0)

    rad_arcsec = np.logspace(np.log10(r_min), np.log10(r_max), num = r_num)

    if shape == 'oblate':
        light_avg = projection.RadialProfile(R=rad_arcsec, profile=profile_hernquist, R_align='major', interpolate=True)
    elif shape == 'prolate':
        light_avg = projection.RadialProfile(R=rad_arcsec, profile=profile_hernquist, R_align='minor', interpolate=True)
    else: 
        raise ValueError("3D shape {} not recongnized!".format(shape))

    rad_arcsec = np.logspace(np.log10(r_min), np.log10(r_max), num = interpolate_num)
    light_avg = light_avg(rad_arcsec)

    ngauss = kwargs_mge.get('ngauss', 20)
    inner_slope = kwargs_mge.get('inner_slope', 2.)
    outer_slope = kwargs_mge.get('outer_slope', 3.)
    quiet_mge = kwargs_mge.get('quiet', True)

    hernquist_mge = mge_fit_1d(rad_arcsec, light_avg, ngauss=ngauss, inner_slope=inner_slope, outer_slope=outer_slope, plot = False, quiet = quiet_mge)

    sigma_hernquist = hernquist_mge.sol[1]
    amp_hernquist = hernquist_mge.sol[0] / (np.sqrt(2 * np.pi) * sigma_hernquist)

    if plot_mge:
        mge_misc.plot_mge(amp_hernquist, sigma_hernquist, rad_arcsec, light_avg)

    if get_profile:
        return amp_hernquist, sigma_hernquist, rad_arcsec, light_avg
    else:
        return amp_hernquist, sigma_hernquist





