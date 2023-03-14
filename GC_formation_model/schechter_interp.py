# Licensed under BSD-3-Clause License - see LICENSE

# Compute maximum cluster mass, Mmax, as a function of the total mass in clusters, Mgc for a schechter function, dN/dM \propto exp(-M/Mc)M^(alpha)

import sys

import numpy as np
from scipy import interpolate
import mpmath

from . import astro_utils

log_mv = np.linspace(3.9, 8.6, num = 2000)
dlog_mv_inv = 1./(log_mv[1] - log_mv[0])
gamma_arr1, gamma_arr2 = np.zeros(len(log_mv)), np.zeros(len(log_mv))
s = 0.02
log_mmaxt = np.arange(4.01, 8.6, step = s)

def upper_gamma2(log_mvv): # linear interpolation of the upper incomplete gamma function for the case of a -2 power law
	return astro_utils.lininterp(log_mvv, log_mv, gamma_arr2, dlog_mv_inv)

def upper_gamma1(log_mvv): # linear interpolation of the upper incomplete gamma function for the case of a -1 power law
	return astro_utils.lininterp(log_mvv, log_mv, gamma_arr1, dlog_mv_inv)

def init(mc, alpha = -2.0):
	for i in range(len(log_mv)):
		mvv = 10**log_mv[i]
		gamma_arr2[i] = mpmath.gammainc(alpha+1.0, mvv/mc) 
		gamma_arr1[i] = mpmath.gammainc(alpha+2.0, mvv/mc) 
	#print "Schechter function init complete"

def generate(mc, alpha = -2.0, mmin = 1e5): # generates functions to interpolate Mgc(M0) and M0(Mmax); combine to give Mgc(Mmax), 
    ug51 = upper_gamma1(np.log10(mmin))
    mgc = mc*np.array([(ug51 - upper_gamma1(log_mmaxtv))/upper_gamma2(log_mmaxtv) for log_mmaxtv in log_mmaxt])
    mgc_to_mmax = interpolate.interp1d(np.log10(mgc), log_mmaxt)
    return mgc_to_mmax
