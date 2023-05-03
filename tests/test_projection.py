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