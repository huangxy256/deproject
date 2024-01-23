import numpy as np 
from deproject.Kinematics.mge_misc import sum_gaussian_components

class MGE_Proj(object):

    def __init__(self, peak_intr, sigma, qintr, oblate = 1):

        self.peak_intr = peak_intr
        self.sigma = sigma 
        self.qintr = qintr
        self._oblate = oblate

    @property
    def oblate(self):
        return self._oblate
    
    @oblate.setter
    def oblate(self, oblate_shape):
        self._oblate = oblate_shape

    def qobs(self, inc):
        if self._oblate:
            return np.sqrt((self.qintr * np.sin(inc))**2 + np.cos(inc)**2) # eq.(35) of Cappellari (2020, MNRAS)
        else:
            return self.qintr / np.sqrt(np.sin(inc)**2 + self.qintr**2 * np.cos(inc)**2) # eq.(36) of Cappellari (2020, MNRAS)

    def surf(self, inc):
        qobs = self.qobs(inc)
        surf = self.peak_intr * np.sqrt(2 * np.pi) * self.sigma * self.qintr / qobs
        return surf

    def mtot(self):
        return (self.qintr * self.peak_intr * (np.sqrt(2 * np.pi) * self.sigma)**3).sum()

    def mtot_proj(self, inc):
        surf = self.surf(inc)
        qobs = self.qobs(inc)
        return (2 * np.pi * surf * self.sigma**2 * qobs).sum()

    def radial_profile(self, inc, rad):
        surf = self.surf(inc)
        density = sum_gaussian_components(x=rad, peaks=surf, sigmas=self.sigma)
        return density

    def radial_profile_circularized(self, inc, rad_circularized):
        surf = self.surf(inc)
        qobs = self.qobs(inc)
        if self._oblate:
            density = sum_gaussian_components(x=rad_circularized/np.sqrt(qobs), peaks=surf, sigmas=self.sigma)
        else:
            density = sum_gaussian_components(x=rad_circularized * np.sqrt(qobs), peaks=surf, sigmas=self.sigma)
        return density