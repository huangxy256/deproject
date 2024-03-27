import numpy as np 
from deproject.Kinematics.mge_misc import sum_gaussian_components

class MGE_Proj(object):

    def __init__(self, peak_intr, sigma, qintr):

        self.peak_intr = peak_intr
        self.sigma = sigma 
        self.qintr = qintr

    def qobs(self, inc):
        if self.qintr <= 1:
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