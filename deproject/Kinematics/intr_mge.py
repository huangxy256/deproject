import numpy as np 
from mgefit.mge_fit_1d import mge_fit_1d

class Intr_MGE_1(object):
    
    def __init__(self, profile, qintr, r_sym):
        """class to manage the intrinsic MGE of axisymmetric profile. The mass baseline is the mass of a spherical profile which has the same sigma as the profile on the symmetry axis of the axisymmetric profile

        Args:
            profile (deproject.Profiles subclass instance): profile instance from the deproject.Profiles module
            qintr (float): intrinsic axis ratio of the axisymmetric profile; if oblate, qintr < 1; if prolate, qintr > 1
            r_sym (array of float): coordinate on the symmetry axis 
        """
        self.profile = profile
        self.qintr = qintr
        self._r_sym = r_sym

    def radial_profile_sph(self):
        """Compute the radial profile for a spherical profile with the given parameterization

        Returns:
            array of float: spherical radial profile at r_sym
        """
        # radial profile for a spherical profile
        rho = self.profile.Density_3d_spherical(r=self._r_sym)
        return rho

    def MGE_param_sph(self, kwargs_mge={}, plot_mge=0, fignum=None):
        """Compute the MGE parameters for a spherical intrinsic radial profile

        Args:
            kwargs_mge (dict, optional): _description_. Defaults to {}.
            plot_mge (bool, optional): whether to plot the MGE components. Defaults to 0.
            fignum (int, optional): figure identification number. Defaults to None.

        Returns:
            arrays of float: peaks of the MGEs
            arrays of float: sigmas of the MGEs
        """
        ngauss = kwargs_mge.get('ngauss', 20)
        inner_slope = kwargs_mge.get('inner_slope', 3.)
        outer_slope = kwargs_mge.get('outer_slope', 1.)
        quiet_mge = kwargs_mge.get('quiet', 1)
        rho = self.radial_profile_sph()
        surf0, sigma = mge_fit_1d(self._r_sym, rho, ngauss=ngauss, inner_slope=inner_slope, outer_slope=outer_slope, plot = plot_mge, quiet = quiet_mge, fignum=fignum).sol
        peak = surf0 / np.sqrt(2 * np.pi) / sigma
        return peak, sigma

    def MGE_mass_sph(self):
        """Compute the total mass/luminosity of the spherical profile

        Returns:
            float: total mass/luminiosity
        """
        peak, sigma = self.MGE_param_sph()
        mtot_sph = (1.0 * peak * (np.sqrt(2 * np.pi) * sigma)**3).sum()
        return mtot_sph
        
    def radial_profile_mass_sph(self):
        """Compute the total mass/luminosity from the sampled radial profile using simpson integration. Note: must have large integration range for precision. 

        Returns:
            float: total mass/luminosity
        """
        rho = self.radial_profile_sph()
        x = self._r_sym
        mass = simpson(x=x, y=rho * x**2 * 4 * np.pi) * 1.0
        return mass
    
    def MGE_param(self, kwargs_mge={}, plot_mge=0, fignum=None):
        """Compute the MGE parameters for the intrinsic axisymmetric density/luminosity profile

        Args:
            kwargs_mge (dict, optional): _description_. Defaults to {}.
            plot_mge (bool, optional): whether to plot the MGE components. Defaults to 0.
            fignum (int, optional): figure identification number. Defaults to None.

        Returns:
            array of float: peaks of the MGEs
            array of float: sigmas of the MGEs
        """
        peak, sigma = self.MGE_param_sph()
        mtot_sph = (1.0 * peak * (np.sqrt(2 * np.pi) * sigma)**3).sum()

        mtot_axi = (self.qintr * peak * (np.sqrt(2 * np.pi) * sigma)**3).sum()
        peak = peak * mtot_sph / mtot_axi
        return peak, sigma



class Intr_MGE_2(object):
    
    def __init__(self, profile, qintr, r_sym):
        """class to manage the intrinsic MGE of axisymmetric profile. The mass baseline is the mass of a spherical profile which has the same peak amplitude as the profile on the symmetry axis of the axisymmetric profile. The intrinsic sigma of the profile on the symmetry axis and the intrinsic sigma of the spherical profile satisfy: sigma_sph = sigma_axi * qintr**(1/3)

        Args:
            profile (_type_): profile instance from the deproject.Profiles module
            qintr (_type_): intrinsic axis ratio of the axisymmetric profile. If oblate, qintr < 1; else qintr > 1
            r_sym (_type_): coordinate on the symmetry axis 
        """
        self.profile = profile
        self.qintr = qintr
        self._r_sym = r_sym

    def radial_profile_sph(self):
        # radial profile for a spherical profile
        rho = self.profile.Density_3d_spherical(r=self._r_sym)
        return rho

    def MGE_param_sph(self, kwargs_mge={}, plot_mge=0, fignum=None):
        ngauss = kwargs_mge.get('ngauss', 20)
        inner_slope = kwargs_mge.get('inner_slope', 3.)
        outer_slope = kwargs_mge.get('outer_slope', 1.)
        quiet_mge = kwargs_mge.get('quiet', 1)
        rho = self.radial_profile_sph()
        surf0, sigma = mge_fit_1d(self._r_sym, rho, ngauss=ngauss, inner_slope=inner_slope, outer_slope=outer_slope, plot = plot_mge, quiet = quiet_mge, fignum=fignum).sol
        peak = surf0 / np.sqrt(2 * np.pi) / sigma
        return peak, sigma

    def MGE_mass_sph(self):
        peak, sigma = self.MGE_param_sph()
        mtot_sph = (1.0 * peak * (np.sqrt(2 * np.pi) * sigma)**3).sum()
        return mtot_sph
        
    def radial_profile_mass_sph(self):
        rho = self.radial_profile_sph()
        x = self._r_sym
        mass = simpson(x=x, y=rho * x**2 * 4 * np.pi) * 1.0
        return mass
    
    def MGE_param(self, kwargs_mge={}, plot_mge=0, fignum=None):
        peak, sigma_sph = self.MGE_param_sph()
        sigma = sigma_sph * self.qintr**(-1/3)
        return peak, sigma

    def MGE_mass(self):
        peak, sigma = self.MGE_param()
        mtot = (self.qintr * peak * (np.sqrt(2 * np.pi) * sigma)**3).sum()
        return mtot
