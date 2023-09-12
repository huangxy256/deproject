import numpy as np
import astropy.constants as const
import astropy.units as u

__all__ = ['SIS']

class SIS(object):

    def __init__(self, unit, sigma = None, theta_E = None, lens_cosmo = None):

        self.unit = unit

        if unit == 'physical':
            if not ((sigma != None) and (theta_E == None)):
                raise ValueError('When unit is physical, input only sigma in km/s!')
            self.sigma = sigma
        elif unit == 'arcsec':
            if not ((theta_E != None) and (sigma == None)):
                raise ValueError('When unit is arcsec, input only theta_E in arcsec!')
            self.theta_E = theta_E
            if lens_cosmo != None:
                self.lens_cosmo = lens_cosmo
        else:
            raise ValueError("Available units: ['physical', 'arcsec'].")

    def Density_3d(self, r, lens_cosmo = None):
        """_summary_

        Args:
            r (_type_): _description_
            lens_cosmo (LensCosmo instance, optional): _description_. Defaults to None.

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        if self.unit == 'physical':
            return ((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G) / (r * u.kpc)**2).to(u.M_sun / u.kpc**3).value
        elif self.unit == 'arcsec':

            if hasattr(self, 'lens_cosmo'):
                if lens_cosmo == None:
                    lens_cosmo = self.lens_cosmo
                else:
                    raise ValueError('Already has lens_cosmo attribute in SIS instance. Do not input lens_cosmo in this function.')
            else:
                if lens_cosmo == None:
                    raise ValueError('No lens_cosmo instance provided. Please specify lens_cosmo in either SIS object construction or pass as a variable to this function.')
                else:
                    pass
            
            return (((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G) / (lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * (360 * 3600 / 2 / np.pi)**2 / r**2

    def _Density_in_project(self, z, as_sq):
        """Integrand of density in LoS integration [M_sun/kpc^3]

        Args:
            z (_type_): _description_
            as_sq (_type_): _description_

        Returns:
            _type_: _description_
        """
        return ((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G * (z**2 + as_sq) * u.kpc**2)).to(u.M_sun / u.kpc**3).value

    def Density_in_project(self, z, as_sq, lens_cosmo = None):
        """ 

        Args:
            z (_type_): _description_
            as_sq (_type_): _description_
            lens_cosmo (_type_, optional): _description_. Defaults to None.

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        if self.unit == 'physical':
            return ((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G * (z**2 + as_sq) * u.kpc**2)).to(u.M_sun / u.kpc**3).value
        # elif self.unit == 'arcsec':

        #     if hasattr(self, 'lens_cosmo'):
        #         if lens_cosmo == None:
        #             lens_cosmo = self.lens_cosmo
        #         else:
        #             raise ValueError('Already has lens_cosmo attribute in SIS instance. Do not input lens_cosmo in this function.')
        #     else:
        #         if lens_cosmo == None:
        #             raise ValueError('No lens_cosmo instance provided. Please specify lens_cosmo in either SIS object construction or pass as a variable to this function.')
        #         else:
        #             pass

        #     rho0 = lens_cosmo.sis_angular_density_norm()
        #     return self.theta_E * rho0 / 1e9 * (2 * np.pi / 360 / 3600) / (z**2 + as_sq**2)
        