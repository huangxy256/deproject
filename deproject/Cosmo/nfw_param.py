import numpy as np

__all__ = ['NFWParam']

class NFWParam(object):
    rhoc = 2.77536627e11  # critical density [h^2 M_sun Mpc^-3]

    def __init__(self, cosmo=None):
        """_summary_

        Args:
            cosmo (object, optional): astropy.cosmology instance. Defaults to None (Planck 2018 cosmology).
        """
        if cosmo == None:
            from astropy.cosmology import default_cosmology
            self.cosmo = default_cosmology.get()
        else:
            self.cosmo = cosmo

    def rhoc_z(self, z):
        """calculate critical density at z

        Args:
            z (float): redshift

        Returns:
            float: critical density at z [h^2 M_sun Mpc^-3]
        """
        return self.rhoc * (self.cosmo.efunc(z))**2

    def r200_M(self, M, z):
        """calculate R_200 for a halo of mass M/h

        Args:
            M (float): mass of the dark matter halo [M_sun/h]
            z (float): redshift of the halo

        Returns:
            float: R_200 [Mpc/h]
        """
        return (3 * M / (4 * np.pi * 200 * self.rhoc_z(z)))**(1.0/3.0)

    def rho0_c(self, c, z):
        """calculate density normalisation as a function of concentration c

        Args:
            c (float): concentration parameter of the DM halo defined with c = R_200/Rs
            z (float): redshift

        Returns:
            float: density normalisation [h^2 M_sun Mpc^-3]
        """
        return (200 * c**3 * self.rhoc_z(z)) / (3 * (np.log(1.0+c) - c/(1.0+c)))

    