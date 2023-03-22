from deproject.parameters import *

def test_Cap_A():

    cap_A = Cap_A(zeta = 1.0, xi = 1.0, theta = 0, phi = 0)

    assert cap_A == 1.0


def test_Cap_B():

    cap_B = Cap_B(zeta = 1.0, xi = 1.0, theta = 0, phi = 0)

    assert cap_B == 0

def test_Cap_C():

    cap_C = Cap_C(zeta = 1.0, xi = 1.0, theta = 0, phi = 0)

    assert cap_C == 1

def test_Small_f():

    small_f = Small_f(zeta = 1.0, xi = 1.0, theta = 0, phi = 0)

    assert small_f == 1