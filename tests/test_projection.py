import numpy.testing as npt
import numpy as np
import pytest

from deproject.Profiles.hernquist import Hernquist
from deproject.projection import Projection
from deproject.Cosmo.default_cosmo import get_default_lens_cosmo
from lenstronomy.LensModel.lens_model import LensModel

class TestProjection(object):
    
    def setup_method(self):
        self.projection = Projection(zeta = 0.5, xi = 0.5, theta = 0, phi = 0)
        self.projection_sph = Projection(zeta=1, xi=1, theta=0, phi=0)

    def test_Cap_A(self):
        npt.assert_almost_equal(self.projection.Cap_A, (0.25)**(-2.0/3.0))

    def test_Cap_B(self):
        npt.assert_almost_equal(self.projection.Cap_B, 0)
    
    def test_Cap_C(self):
        npt.assert_almost_equal(self.projection.Cap_C, (0.25)**(1.0/3.0))

    def test_Lcase_f(self):
        npt.assert_almost_equal(self.projection.Lcase_f, (0.25)**(-1.0/3.0))

    def test_AxisRatio(self):
        npt.assert_almost_equal(self.projection.AxisRatio(), 0.5)

    def test_Ellipticity(self):
        npt.assert_almost_equal(self.projection.Ellipticity(), 1.0/3.0)

    def test_Orientation_phi(self):
        pass

    def test_RadialProfile(self):
        R = np.logspace(np.log10(0.001), np.log10(100), num = 50)
        lens_cosmo = get_default_lens_cosmo()
        M_star = 1e11
        Rs_Mpc = 0.008
        sigma0, Rs_angle_star = lens_cosmo.hernquist_physical2angle(M=M_star, Rs=Rs_Mpc)
        hernquist = Hernquist(Rs=Rs_angle_star, sigma0=sigma0)
        radial_profile = self.projection_sph.RadialProfile(R=R, profile=hernquist)[0]

        lens_list = ['HERNQUIST']
        lens_model = LensModel(lens_list)
        kwargs_lens = [{'sigma0': sigma0, 'Rs': Rs_angle_star}]
        radial_kappa = lens_model.kappa(x=R, y=0, kwargs = kwargs_lens)

        npt.assert_almost_equal(radial_profile, radial_kappa, decimal = 6)

    def test_Project_2d(self):
        x = np.linspace(0.01, 80, num = 50)
        y = np.linspace(0.01, 80, num = 50)
        xx, yy = np.meshgrid(x, y)

        lens_cosmo = get_default_lens_cosmo()
        M_star = 1e11
        Rs_Mpc = 0.008
        sigma0, Rs_angle_star = lens_cosmo.hernquist_physical2angle(M=M_star, Rs=Rs_Mpc)
        hernquist = Hernquist(Rs=Rs_angle_star, sigma0=sigma0)
        kappa_numer = self.projection_sph.Project_2d(x_grid=xx, y_grid=yy, profile=hernquist)

        lens_list = ['HERNQUIST']
        lens_model = LensModel(lens_list)
        kwargs_lens = [{'sigma0': sigma0, 'Rs': Rs_angle_star}]
        kappa_lenst = lens_model.kappa(x = xx, y=yy, kwargs = kwargs_lens)

        npt.assert_almost_equal(kappa_numer, kappa_lenst, decimal = 6)


if __name__ == '__main__':
    pytest.main()

