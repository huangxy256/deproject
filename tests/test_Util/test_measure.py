import numpy as np
import numpy.testing as npt
import pytest

from lenstronomy.LensModel.lens_model import LensModel
from deproject.Util.measure import *

def test_einstein_radius_from_rp():

    lens_list = ['SIS']
    theta_E_true = np.random.randint(1, 10, 1)[0]
    lens_kwargs = [{"theta_E": theta_E_true, "center_x": 0, "center_y": 0}]
    lens_model = LensModel(lens_model_list=lens_list)
    R_circ = np.geomspace(1e-3, 1e2, 50)
    kappa1d = lens_model.kappa(R_circ, 0, kwargs=lens_kwargs)

    # interpolation within range
    theta_E_measure = einstein_radius_from_rp(kappa1d, R_circ, R_min=1e-3, R_max=1e1)
    npt.assert_almost_equal(theta_E_measure, theta_E_true, decimal=3)
    # interpolation out of range
    theta_E_measure = einstein_radius_from_rp(kappa1d, R_circ, R_min=1e-4, R_max=1e3)
    npt.assert_almost_equal(theta_E_measure, theta_E_true, decimal=3)
    # bad lens
    lens_kwargs = [{"theta_E": 100, "center_x": 0, "center_y": 0}]
    kappa1d = lens_model.kappa(R_circ, 0, kwargs=lens_kwargs)
    theta_E_measure = einstein_radius_from_rp(kappa1d, R_circ, R_min=1e-3, R_max=1e1)
    assert np.isnan(theta_E_measure)