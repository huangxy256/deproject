import numpy.testing as npt
import numpy as np

import deproject.Profiles.Hernquist
import deproject.projection

from lenstronomy.LightModel.light_model import LightModel

def test_Project_2d():
    
    rho0 = 8.0
    Rs = 1.0

    profile = deproject.Profiles.Hernquist.Hernquist(rho0=rho0, Rs=Rs)

    zeta = 1.0
    xi = 1.0
    theta = 0.5
    phi = 0.6

    x = np.linspace(-10, 10, num = 100)
    y = np.linspace(-10, 10, num = 100)

    xx, yy = np.meshgrid(x, y)

    sigma = deproject.projection.Project_2d(x_grid=xx, y_grid=yy, zeta=zeta, xi=xi, phi=phi, theta=theta, profile=profile)

    e1_he = 0.
    e2_he = 0.

    light_model = LightModel(light_model_list = ['HERNQUIST_ELLIPSE'])
    kwargs_light = [{'amp': rho0, 'Rs': Rs, 'e1': -e1_he, 'e2': e2_he}]
    sigma_lenstronomy = light_model.surface_brightness(xx, yy, kwargs_light)

    npt.assert_almost_equal(sigma, sigma_lenstronomy, decimal = 5)


def test_RadialProfile():

    rho0 = 8.0
    Rs = 1.0

    profile = deproject.Profiles.Hernquist.Hernquist(rho0=rho0, Rs=Rs)

    zeta = 1.0
    xi = 1.0
    theta = 0.5
    phi = 0.6

    xpp = np.logspace(np.log10(0.01), np.log10(10), num = 100)

    sigma_rad = deproject.projection.RadialProfile(x_major=xpp, zeta=zeta, xi=xi, theta=theta, phi=phi, profile=profile)

    e1_he = 0.
    e2_he = 0.

    light_model = LightModel(light_model_list = ['HERNQUIST_ELLIPSE'])
    kwargs_light = [{'amp': rho0, 'Rs': Rs, 'e1': -e1_he, 'e2': e2_he}]
    sigma_rad_lenstronomy = light_model.surface_brightness(x=xpp, y=0, kwargs_list=kwargs_light)

    npt.assert_almost_equal(sigma_rad, sigma_rad_lenstronomy, decimal = 3)


def test_RadialProfile_any():

    # todo: complete this part

    # convert the following into test functions of RadialProfile_any

    # theta = 1.3
    # phi = 4.5

    # rad = np.logspace(np.log10(0.1), np.log10(10), num = 50)
    # radial_average = RadialProfile(rad, zeta = zeta_single, xi = xi_single, theta = theta, phi = phi, profile = profile_sis)

    # plt.loglog(radial_average[1], radial_average[0], label = 'average')

    # radial_average_new = RadialProfile_any(x = rad, x_type= 'average', zeta = zeta_single, xi = xi_single, theta = theta, phi = phi, profile = profile_sis)

    # plt.loglog(radial_average_new[1], radial_average_new[0], label = 'average_new')

    # radial_major = RadialProfile_any(x = rad, x_type= 'major', zeta = zeta_single, xi = xi_single, theta = theta, phi = phi, profile = profile_sis)

    # radial_minor = RadialProfile_any(x = rad, x_type= 'minor', zeta = zeta_single, xi = xi_single, theta = theta, phi = phi, profile = profile_sis)

    # axis_ratio = AxisRatio(zeta = zeta_single, xi = xi_single, theta = theta, phi = phi)

    # plt.loglog(radial_major[1], radial_major[0], label = 'major')
    # plt.loglog(radial_minor[1], radial_minor[0], label = 'minor')

    # plt.loglog(radial_major[1] * axis_ratio, radial_major[0], label = 'minor from major', ls = '--')

    # plt.loglog(radial_major[1] * np.sqrt(axis_ratio), radial_major[0], ls = '-.', label = 'average from major')

    # plt.loglog(radial_minor[1] / np.sqrt(axis_ratio), radial_minor[0], ls = ':', label = 'average from minor')

    # print(axis_ratio)

    # plt.legend()