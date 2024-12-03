import numpy as np 
import numpy.testing as npt
import pytest

from deproject.MGE_analysis.mge_proj import MGE_Proj

class TestMGE_proj(object):
    def setup_method(self):
        self.peak_intr = np.array([2.20210302e+06, 3.67627885e+05, 8.53917638e+04, 2.10410004e+04,
        4.96276754e+03, 1.17985733e+03, 2.80693649e+02, 6.67230673e+01,
        1.58651606e+01, 3.77209603e+00, 8.96871202e-01, 2.13239793e-01,
        5.06998587e-02, 1.20582122e-02, 2.86249789e-03, 6.82997181e-04,
        1.62970230e-04, 3.61941839e-05, 6.52052159e-06, 4.95437631e-07])
        self.sigma_intr = np.array([9.70963979e-01, 2.16660416e+00, 4.39529445e+00, 8.97480038e+00,
        1.84430923e+01, 3.78111269e+01, 7.75425126e+01, 1.59030404e+02,
        3.26140044e+02, 6.68856955e+02, 1.37170552e+03, 2.81314527e+03,
        5.76910594e+03, 1.18328701e+04, 2.42684669e+04, 4.97269852e+04,
        1.02672783e+05, 2.10859277e+05, 4.16521866e+05, 1.01383345e+06])
        self.qintr_oblate = 0.5
        self.qintr_prolate = 1.2
        self.MGE_proj_oblate = MGE_Proj(self.peak_intr, self.sigma_intr, self.qintr_oblate)
        self.MGE_proj_prolate = MGE_Proj(self.peak_intr, self.sigma_intr, self.qintr_prolate)

    def test_qobs(self):
        npt.assert_almost_equal(self.MGE_proj_oblate.qobs(0.), 1.)
        npt.assert_almost_equal(self.MGE_proj_prolate.qobs(0.), 1.)

        npt.assert_almost_equal(self.MGE_proj_oblate.qobs(np.pi/2), self.qintr_oblate)
        npt.assert_almost_equal(self.MGE_proj_prolate.qobs(np.pi/2), self.qintr_prolate)

    def test_surf(self):
        surf_oblate = self.MGE_proj_oblate.surf(np.pi/2)
        npt.assert_almost_equal(surf_oblate, self.peak_intr * self.sigma_intr * np.sqrt(np.pi*2))

        surf_prolate = self.MGE_proj_prolate.surf(np.pi/2)
        npt.assert_almost_equal(surf_prolate, self.peak_intr * self.sigma_intr * np.sqrt(np.pi*2))

    def test_mtot_proj(self):
        mtot = self.MGE_proj_oblate.mtot() / 1e13
        mtot_proj = self.MGE_proj_oblate.mtot_proj(50.) / 1e13
        npt.assert_almost_equal(mtot_proj, mtot)

        mtot = self.MGE_proj_prolate.mtot() / 1e13
        mtot_proj = self.MGE_proj_prolate.mtot_proj(80.) / 1e13
        npt.assert_almost_equal(mtot_proj, mtot)
        
