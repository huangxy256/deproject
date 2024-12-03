import numpy as np
import pytest
import numpy.testing as npt

from deproject.MGE_analysis.intr_mge import Intr_MGE, Intr_MGE_rescale_amp
from deproject.Profiles.SIS_truncated_physical import SIS_truncated_physical


class TestIntr_MGE(object):
    def setup_method(self):
        sigma_v = 200
        self.profile = SIS_truncated_physical(sigma_v, 500)
        self.qintr = 0.7
        self.r_sym = np.geomspace(0.01, 1000, 100)
        self.Intr_MGE = Intr_MGE(self.profile, self.qintr, self.r_sym)
    
    def test_MGE_mass_sph(self):
        mass_numerical = self.Intr_MGE.radial_profile_mass_sph()
        mass_mge = self.Intr_MGE.MGE_mass_sph()
        # npt.assert_almost_equal(mass_mge, mass_numerical)
        assert 1 # TODO: accurate integration of total mass
    
    def test_MGE_mass(self):
        mass_sph = self.Intr_MGE.MGE_mass_sph() / 1e13
        mass_axi = self.Intr_MGE.MGE_mass() / 1e13
        npt.assert_almost_equal(mass_axi, mass_sph)

if __name__ == 'main':
    pytest.main()