from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3, Ob0 = 0.05)
from deproject.Cosmo.lens_cosmo import LensCosmo
lens_cosmo = LensCosmo(z_lens = 0.5, z_source = 1.5, cosmo = cosmo)

def get_default_lens_cosmo():
    return lens_cosmo

def get_default_cosmo():
    return cosmo