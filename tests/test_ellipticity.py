from deproject.ellipticity import *

def test_Ellipticity():

    test_ellipticity = Ellipticity(zeta = 1, xi = 1, theta = 0, phi = 0)

    assert test_ellipticity == 0