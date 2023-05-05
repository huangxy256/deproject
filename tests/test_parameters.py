from deproject.parameters import *
import numpy.testing as npt

def test_Cap_A():

    cap_A = Cap_A(zeta = 1.0, xi = 1.0, theta = 0.78, phi = 5.6)

    npt.assert_almost_equal(cap_A, 1.0, decimal = 6)


def test_Cap_B():

    cap_B = Cap_B(zeta = 1.0, xi = 1.0, theta = 0.6, phi = 3.4)

    npt.assert_almost_equal(cap_B, 0.0, decimal = 6)

def test_Cap_C():

    cap_C = Cap_C(zeta = 1.0, xi = 1.0, theta = 1.23, phi = 2.45)

    npt.assert_almost_equal(cap_C, 1.0, decimal = 6)

def test_Small_f():

    small_f = Small_f(zeta = 1.0, xi = 1.0, theta = 1.67, phi = 3.8)

    npt.assert_almost_equal(small_f, 1.0, decimal = 6)
