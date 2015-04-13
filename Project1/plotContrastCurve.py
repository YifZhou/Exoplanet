#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d

if __name__ == '__main__':
    inFN = sys.argv[1]
    sep, curve = np.loadtxt(inFN, unpack = True)
    new_sep = np.linspace(sep.min(), sep.max(), 100)
    new_curve = interp1d(sep, curve, kind='slinear')(new_sep)
    plt.close('all')
    plt.plot(new_sep, new_curve)
    plt.xlabel('Angular Separation (arcsec)')
    plt.ylabel('S/N')
    plt.savefig(sys.argv[2])
    