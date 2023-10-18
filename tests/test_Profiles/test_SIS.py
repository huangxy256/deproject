import numpy as np
import numpy.testing as npt
from deproject.Profiles.SIS import SIS

from deproject.Cosmo.lens_cosmo import LensCosmo
from astropy.cosmology import FlatLambdaCDM


def test_Density_3d():
    sigma = 200.

    z_lens = 0.5
    z_source = 1.5
    cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3, Ob0 = 0.05)
    lens_cosmo = LensCosmo(z_lens=z_lens, z_source=z_source, cosmo=cosmo)

    theta_E = lens_cosmo.sis_sigma_v2theta_E(v_sigma=sigma)

    sis_phys = SIS(sigma=sigma, unit = 'physical')
    sis_angular = SIS(theta_E=theta_E, unit = 'arcsec')

    r_kpc = np.logspace(np.log10(0.001), np.log10(20), num = 50)
    r_arcsec = r_kpc / 1000. / lens_cosmo.Dd / ((2 * np.pi) / 360 / 3600)

    den_phys = sis_phys.Density_3d(r = r_kpc)
    den_ang = sis_angular.Density_3d(lens_cosmo=lens_cosmo, r = r_arcsec)

    npt.assert_allclose(den_phys, den_ang, rtol = 1e-3)


