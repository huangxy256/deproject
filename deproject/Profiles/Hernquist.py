import numpy as np

__all__ = ['Hernquist']

class Hernquist:

    def __init__(self, Rs, sigma0):
        """The Hernquist profile class

        Args:
            Rs (_type_): _description_
            sigma0 (_type_): _description_
        """
        rho0 = self.sigma02rho(sigma0=sigma0, Rs=Rs)
        self.rho0 = rho0
        self.Rs = Rs

    @staticmethod
    def sigma02rho(sigma0, Rs):
        """convert sigma0 to 3d density normalisation

        Args:
            sigma0 (float): projected surface density [Sigma_crit]
            Rs (float): scale radius [arcsec]

        Returns:
            float: 3d density normalisation [Sigma_crit arcsec^-1]
        """
        return sigma0 / Rs


    def Density_3d_spherical(self, r):
        """3d density 

        Args:
            r (_type_): radius [arcsec]

        Returns:
            _type_: _description_
        """
        return self.rho0 / ((r / self.Rs) * (1 + (r / self.Rs))**3)

    def Density_3d_triaxial(self, x, y, z, zeta, xi, get_effective_radius = False):
        """3d density

        Args:
            x (_type_): x coordinate, aligned with intrinsic major axis
            y (_type_): y coordinate, aligned with intrinsic intermediate axis
            z (_type_): z xoordinate, aligned with intrinsic minor axis
            zeta (_type_): intrinsic axis ratio b/a
            xi (_type_): intrinsic axis ratio c/a
            get_effective_radius (bool, optional): whether to return effective radius. Defaults to False.

        Returns:
            _type_: _description_
        """
        av = (zeta * xi)**(1/3) * np.sqrt(x**2 + y**2/zeta**2 + z**2/xi**2)
        density = self.rho0 / ((av / self.Rs) * (1 + (av / self.Rs))**3)

        if get_effective_radius:
            return density, av
        else:
            return density

    def Project_integrand(self, z, as_sq):
        """Integrand in projection integral, z to be integrated 

        Args:
            z (float): line of sight variable
            as_sq (float): square of radial variable in projection

        Returns:
            float: density in projection
        """
        return self.rho0 / ((np.sqrt(z**2 + as_sq)/ self.Rs) * (1 + (np.sqrt(z**2 + as_sq)/ self.Rs))**3)