import numpy as np 
from mgefit.mge_fit_1d import mge_fit_1d

class Intr_MGE(object):

    def __init__(self, profile, qintr, oblate, r_sym):
        self.profile = profile
        self.qintr = qintr
        self._oblate = oblate
        self._r_sym = r_sym

    @property
    def oblate(self):
        return self._oblate
    
    @oblate.setter
    def oblate(self, oblate_shape):
        self._oblate = oblate_shape

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
        peak, sigma = self.MGE_param_sph()
        mtot_sph = (1.0 * peak * (np.sqrt(2 * np.pi) * sigma)**3).sum()

        if self._oblate:
            qintr1 = self.qintr
        else:
            qintr1 = 1 / self.qintr
        
        mtot_axi = (qintr1 * peak * (np.sqrt(2 * np.pi) * sigma)**3).sum()
        peak = peak * mtot_sph / mtot_axi
        return peak, sigma