import numpy as np
import deproject.constant as const
from deproject.Cosmo.nfw_param import NFWParam

class LensCosmo(object):
    def __init__(self, z_lens, z_source, cosmo):
        """class to manage the physical units and distances present in a single plane lens with fixed input cosmology

        Args:
            z_lens (float): redshift of the lens
            z_source (float): redshift of the lensing source
            cosmo (astropy.cosmology instance): _description_
        """
        self.z_lens = z_lens
        self.z_source = z_source
        self.cosmo = cosmo
        self.nfw_param = NFWParam(cosmo=cosmo)

    @property
    def h(self):
        """dimensionless Hubble constant

        Returns:
            float: value of the dimensionless Hubble constant
        """
        return self.cosmo.H0.value / 100.0

    @property
    def Dds(self):
        """angular diameter distance between lens and source

        Returns:
            float: value of angular diameter distance between lens and source in the given cosmology [Mpc]
        """
        return (self.cosmo.angular_diameter_distance_z1z2(self.z_lens, self.z_source)).value

    @property
    def Dd(self):
        """angular diameter distance of the lens

        Returns:
            float: value of angular diameter distance between observer (z=0) and the lens in the given cosmology [Mpc]
        """
        return (self.cosmo.angular_diameter_distance(self.z_lens)).value

    @property
    def Ds(self):
        """angular diameter distance of the source

        Returns:
            float: value of angular diameter distance between observer (z=0) and the source in the given cosmology [Mpc]
        """
        return (self.cosmo.angular_diameter_distance(self.z_source)).value

    @property
    def Sigma_crit(self):
        """calculate Sigma_crit

        Returns:
            float: Sigma_crit [M_sun Mpc^-2]
        """
        if not hasattr(self, '_sigma_crit_mpc'):
            const_SI = const.c**2 / (4 * np.pi * const.G) # [kg m^-1]
            conversion = const.Mpc / const.M_sun # [m M_sun kg^-1 Mpc^-1]
            factor = const_SI * conversion # [M_sun Mpc^-1]
            self._sigma_crit_mpc = self.Ds / (self.Dd * self.Dds) * factor
        return self._sigma_crit_mpc

    def Mpc2arcsec_lens(self, physical_Mpc):
        """convert physical distance in lens plane to arcseconds

        Args:
            physical_Mpc (float): physical distance in lens plane [Mpc]

        Returns:
            float: angular diameter [arcsec]
        """
        return physical_Mpc / self.Dd / const.arcsec

    def arcsec2Mpc_lens(self, arcsec):
        """convert angular diameter to physical distance in Mpc

        Args:
            arcsec (float): angular diameter [arcsec]

        Returns:
            float: physical distance in lens plane [Mpc]
        """
        return arcsec * const.arcsec * self.Dd

    def nfwParam_physical(self, M, c):
        """calculate physical NFW parameters: rho0, Rs, R200

        Args:
            M (float): physical mass [M_sun]
            c (float): concentration parameter
        Returns:
            float, float, float: rho0 [M_sun Mpc^-3], Rs [Mpc], R200 [Mpc]
        """
        r200 = self.nfw_param.r200_M(M * self.h, self.z_lens) / self.h
        rho0 = self.nfw_param.rho0_c(c, self.z_lens) * self.h**2
        Rs = r200 / c
        return rho0, Rs, r200

    def nfw_physical2angle(self, M, c):
        """calculate Rs and alphs_Rs in angle units from physical mass and concentration

        Args:
            M (float): physical mass of NFW halo [M_sun]
            c (float): concentration parameter

        Returns:
            float, float: Rs [arcsec], alpha_Rs [arcsec]
        """
        rho0, Rs, r200 = self.nfwParam_physical(M, c)
        Rs_angle = Rs / self.Dd / const.arcsec # [arcsec]
        alpha_Rs = rho0 * (4 * Rs ** 2 * (1 + np.log(1. / 2.)))
        alpha_Rs = alpha_Rs / self.Sigma_crit / self.Dd / const.arcsec
        return Rs_angle, alpha_Rs

    def hernquist_physical2angle(self, M, Rs):
        """calculate sigma0 and Rs_angle in angular units from mass and physical Rs

        Args:
            M (float): spherical overdensity mass defined with mdef at redshift z [M_sun]
            Rs (float): scale radius [Mpc]

        Returns:
            float, float: sigma0 [Sigma_crit], Rs_angle [arcsec]
        """
        Rs_angle = Rs / self.Dd / const.arcsec # [arcsec]
        rhos = M / (2 * np.pi) / Rs**3 # [M_sun Mpc^-3]
        sigma0 = rhos * Rs / self.Sigma_crit
        return sigma0, Rs_angle

    def sis_sigma_v2theta_E(self, v_sigma):
        """convert velocity dispersion to Einstein radius for a SIS profile

        Args:
            v_sigma (float): velocity dispersion [km s^-1]

        Returns:
            float: Einstein radius [arcsec]
        """
        theta_E = 4 * np.pi * (v_sigma * 1000. / const.c)**2 * self.Dds / self.Ds / const.arcsec
        return theta_E

    def sis_sigma_v_sq2theta_E(self, v_sigma_sq):
        """convert velocity dispersion square to Einstein radius for a SIS profile 

        Args:
            v_sigma_sq (float): square of velocity dispersion [km^2 s^-2]

        Returns:
            float: Einstein radius [arcsec]
        """
        theta_E = 4 * np.pi * v_sigma_sq * 1000**2 / (const.c)**2 * self.Dds / self.Ds / const.arcsec
        return theta_E

    def sis_theta_E2sigma_v(self, theta_E):
        """convert Einstein radius to velocity dispersion for a SIS profile

        Args:
            theta_E (float): Einstein radius [arcsec]

        Returns:
            float: velocity dispersion [km/s]
        """
        sigma_v = np.sqrt(theta_E * const.arcsec / (4 * np.pi) * self.Ds / self.Dds * const.c**2) / 1000
        return sigma_v

