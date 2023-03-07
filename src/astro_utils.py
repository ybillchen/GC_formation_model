import numpy as np
from scipy import interpolate
# from astropy.cosmology import WMAP9 as cosmo # Use same cosmology as Illustris
from astropy.cosmology import Planck15 as cosmo # Use same cosmology as TNG
import astropy.units as u

#use WMAP-9 cosmology, as adopted in (original, not TNG) Illustris
# H0 = 7.195e-5 # h = 0.704 in Myr^-1
# omega_matter = 0.2726
H0 = cosmo.H0.to(1/u.Myr).value
omega_matter = cosmo.Om0
omega_lambda = 1. - omega_matter
ratio = omega_matter/omega_lambda
G = 0.00449972292 #pc^3 /(Msun * Myr^2)
#Behroozi+ 13 stellar mass-halo mass relation, plus some other useful functions
def fSMHM(x, alpha, delta, gamma) : 
    if(x > -2.):
        val = -np.log10(10**(alpha*x)+1) + delta*((np.log10(1+np.exp(x)))**gamma)/(1+np.exp(10**(-x)))
    else :
        val = -np.log10(10**(alpha*x)+1)
    return val
def SMHMparameters(z) : #standard Behroozi relation (Mvir)
    a = 1./(1.+z)
    nu = np.exp(-4*(a**2))

    log_eps = -1.777+((-0.006*(a-1))*nu) - (0.119*(a-1))
    eps = 10**log_eps

    log_m1 = 11.514 + (-1.793*(a-1)+ (-0.251*z))*nu
    m1 = 10**log_m1

    alpha = -1.412 + (.713*(a-1))*nu
    delta = 3.508 + (2.608*(a-1) + (-0.043*z))*nu
    gamma = 0.316 + (1.319*(a-1) + (0.279*z))*nu
    xi = np.random.normal(0,.218 - .023*(a-1)) #scatter
    return m1, alpha, delta, gamma, eps, xi
def SMHMparameters2(z): #Kravtsov+ 2014 M200 (no scatter) 
    a = 1./(1.+z)
    nu = np.exp(-4*(a**2))
    log_eps = -1.618
    eps = 10**log_eps

    log_m1 = 11.39
    m1 = 10**log_m1

    alpha = -1.795
    delta = 4.345
    gamma = 0.619

    xi = np.random.normal(0,.218 - .023*(a-1)) #scatter
    return m1, alpha, delta, gamma, eps, xi
def SMHMparameters3(z): #Kravtsov+ 2014 Mvir (no scatter)
    a = 1./(1+z)
    nu = np.exp(-4*a*a)
    log_eps = -1.663
    eps = 10**log_eps

    log_m1 = 11.43
    m1 = 10**log_m1

    alpha = -1.750
    delta = 4.290
    gamma = 0.595

    xi = np.random.normal(0,.218 - .023*(a-1)) #scatter
    return m1, alpha, delta, gamma, eps, xi
def SMHM(Mh, z, k = False, scatter = False, mdef = 'm200') : 
    if(not k):
        m1, alpha, delta, gamma, eps, xi = SMHMparameters(z) 
    else:
        if(mdef == 'm200'):
            m1, alpha, delta, gamma, eps, xi = SMHMparameters2(z) 
        else: #Mvir
            m1, alpha, delta, gamma, eps, xi = SMHMparameters3(z)
    logSM  = np.log10(eps*m1)+fSMHM(np.log10(Mh/m1), alpha, delta, gamma) - fSMHM(0, alpha, delta, gamma) 
    if(scatter):
    	logSM += xi
    SM = 10**logSM
    return SM
def lininterp(xarb,x, y, d_inv):
    f = (xarb - x[0])*d_inv;
    b = int(f);
    f -= b;
    return y[b] + f*(y[b+1] - y[b]);
dz = .0001
zt = np.arange(-0.2, 1000, step = dz)
#compute cosmic time(z)
f = (omega_matter/(1-omega_matter))*((1+zt)**3)
t = 2./3./np.sqrt(1-omega_matter)*np.log((1.+np.sqrt(1.+f))/np.sqrt(f)) #in units of 1/H0
ct = t/H0
#compute E(z)
ezt = (omega_matter*((1+zt)**3) + omega_lambda)**0.5
#compute evolving virial overdensity (Bryan & Norman 1998)
omega = omega_matter*(1+zt)**3 / (ezt*ezt)
x = omega - 1.0
ozt = 18*np.pi*np.pi + 82*x - 39*x*x
def cosmicTime(z, units = 'Myr') : #returns time(z) in Myr
    if(hasattr(z, '__len__')):
        t = np.array([lininterp(zv, zt, ct, 1./dz) for zv in z])
    else:
        t = lininterp(z, zt, ct, 1./dz)
    if(units == 'Gyr'):
        t = t/1e3
    if(units == 'yr'):
        t = t*1e6
    return t
def timeToRedshift(time, units = 'Myr') :    
    if(units == "Gyr"):
        time *= 1e3
    elif(units == "yr"):
        time /= 1e6
    x = 1.5*(omega_lambda)**.5 * H0*time
    a = (ratio)**(1./3) * (np.sinh(x) ** (2./3))
    return 1./a - 1
def thub(z, units = 'yr'): #return 1/H(z) in years
    Ez = E(z)
    H_inv_Myr = H0*Ez
    H_inv_yr = H_inv_Myr/1e6
    thub = 1./H_inv_yr
    if(units == 'yr'):
        return thub
    elif(units == 'Myr'):
        return thub/1e6
    elif(units == 'Gyr'):
        return thub/1e9
def E(z): 
    return lininterp(z, zt, ezt, 1./dz)
def overdensity(z):
    return lininterp(z, zt, ozt, 1./dz)
def bulge_mass(sm, z): #M_bulge(M*, z), sm in log units
    a = 1./(1+z)
    z_mul = 1.0 - 0.5*(1.0-a)
    return sm+np.log10(z_mul/(1.0+np.exp(-1.12882*(sm-10.1993)))) #returns m_bulge in 
def virialRadius(m, z):
    ez = E(z)
    rhocrit = 3*H0*H0*ez*ez/(8*np.pi*G)
    return (3*m/(4*np.pi*overdensity(z)*rhocrit))**(1./3)
def distance(x1, y1, z1, x2, y2, z2): #for r/rvir computation
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) 
def sigmaDM(mh, rvir):
    return np.sqrt(G*mh/(2*rvir))
def vvir(mh, z):
    rvir = virialRadius(mh, z)
    return np.sqrt(G*mh/rvir)
def tvir(mh, z):
    rvir = virialRadius(mh,z)
    vvir = np.sqrt(G*mh/rvir)
    return 3.6e5*(vvir*.9785/100)**2 #from vdB lecture slides
def find_nearest(array,value, return_idx = False):
    if(hasattr(value, '__len__')):
        l = []
        for v in value:
            idx = (np.abs(array-v)).argmin()
            if(return_idx):
                l.append(idx)
            else:
                l.append(array[idx])
        return l
    else:
        idx = (np.abs(array-value)).argmin()
        if(return_idx):
            return idx
        else:
            return array[idx]
#f = '/Users/nchoksi/gdrive/Astro/smbh_mergers/sfr_charlie.txt'
def mag_to_flux(m): #return flux density in Jy given magnitude
    return 10**((23.9 - m)/2.5)/1e6  #https://en.wikipedia.org/wiki/Jansky
def flux_to_mag(f, zp_jy = 3631.0): #return AB magnitude given flux
    return -2.5*np.log10(f/zp_jy) #https://en.wikipedia.org/wiki/AB_magnitude
def sample_power_law(r, xmin, alpha): #assume xmax extends to infinity; alpha is the power-law index
    return xmin*(1-r)**(1./(alpha+1))
def sample_power_law2(r, xmin, xmax, alpha): #more accurate version of above; alpha is the power-law index
    return (r*(xmax**(alpha+1) - xmin**(alpha+1)) + xmin**(alpha+1))**(1./(alpha+1))

def addEnds(x,p25, minx, maxx): #p25 = y values
    slope_p25_low = (p25[1] - p25[0])/(x[1] - x[0])
    
    slope_p25_high = (p25[-1] - p25[-2])/(x[-1] - x[-2])
    
    p25_low = p25[0] + slope_p25_low*(-x[0] + minx)
    
    p25_high = p25[-1] + slope_p25_high*(-x[-1] + maxx)
    
    p25 = np.insert(p25,0,  p25_low)
    p25 = np.append(p25, p25_high)
    x = np.insert(x, 0, minx)
    x = np.append(x, maxx)
    
    return x, p25
def kde_gauss(data, x, h, weights = None):
    if(weights is None):
        weights = np.array([1.0 for v in data])
    y = weights[:,None]*np.exp(-((x-data[:,None])/h)**2/2)
    return y.sum(0)/(np.sqrt(2*np.pi)*h*len(data))
