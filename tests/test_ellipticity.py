from deproject.ellipticity import *
import numpy.testing as npt


def test_Ellipticity():

    test_ellipticity = Ellipticity(zeta = 1, xi = 1, theta = 0, phi = 0)

    assert test_ellipticity == 0

# def test_Ellipticity_from_grid():

#     x = np.arange(-50, 50, 1.0)
#     y = np.arange(-50, 50, 1.0)
#     xx, yy = np.meshgrid(x, y)

#     # theoretical power-law profile
#     s = 0.0000000001
#     Sigma = np.pi / np.sqrt(xx**2 + yy**2 + s**2)

#     test_ellipticity_12 = Ellipticity_12(xx, yy, Sigma)

#     npt.assert_almost_equal(test_ellipticity_12, np.zeros(len(test_ellipticity_12)), decimal = 3)
