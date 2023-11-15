from scipy.stats import truncnorm
import numpy as np 

def generate_bins_from_uniform(low, high, size_1d):
    return np.random.uniform(low = low, high = high, size = [size_1d, 2]).T


def generate_bins_from_normal(low_bound, up_bound, size_1d, xmean, xscale, ymean, yscale):
    xlow = low_bound[0]
    ylow = low_bound[1]

    xhigh = up_bound[0]
    yhigh = up_bound[1]

    ax, bx = (xlow - xmean) / xscale, (xhigh - xmean) / xscale
    ay, by = (ylow - ymean) / yscale, (yhigh - ymean) / yscale

    xbin = truncnorm.rvs(ax, bx, size = size_1d, loc = xmean, scale = xscale)
    ybin = truncnorm.rvs(ay, by, size = size_1d, loc = ymean, scale = yscale)

    return xbin, ybin