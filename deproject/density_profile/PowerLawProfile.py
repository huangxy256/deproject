

class PowerLawProfile:

    def __init__(self, r, rho0, r0, alpha):
        self.r = r
        self.rho0 = rho0
        self.r0 = r0
        self.alpha = alpha

    def Density(self):
        """calculate density profile

        Returns:
            float: density profile
        """
        return rho0 * (r/r0)**(-alpha)

    