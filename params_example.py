# Licensed under BSD-3-Clause License - see LICENSE

import numpy as np

# disrupt mode: constant, tidal
disrupt_mode = 'tidal'

# cosmological parameters
h100 = 0.6774
Ob = 0.0486
Om = 0.3089

# adjustable model parameters (run mode)
p2 = 18
p3 = 0.5
kappa = 2.5

# grid length in calculating tidal tensor, in kpc
d_tid = 0.3 

# x and y parameters if disrupt_mode == 'tidal'
disrupt_x = 2./3.
disrupt_y = 4./3.

# Normalized period of rotation for disruption if disrupt_mode == 'constant'
pr = 0.5 

# log(Mcut) - knee of the Schechter function. Anything with log_mc >
# 10 will be treated as a pure power-law, to speed up runtime
# (analytic solution).
log_mc = 7.0 

seed = 0 # random seed

# Output file names: will result in cat_base_mcvalue_p2_p3.txt ,
# allcat_base_mcvalue_p2_p3.txt, merits_base.txt
allcat_base = 'allcat'
merit_name = 'merits.txt'

# Optional model parameters. These are usually taken as fixed, but in
# principle can vary
mpb_only = False

# Scaling relations parameters, defined in CGL18
mmr_slope = 0.35 # Mass slope of galaxy mass-metallicity relation
mmr_pivot = 10.5 # Pivot in the power-law of mass-metallicity relation
mmr_evolution = 0.9 # Redshift slope of the galaxy mass-metallicity relation
max_feh = 0.3 # Max [Fe/H]

# Scaling of the gas depletion time with redshift, tdep \propto
# (1+z)^(-alpha), alpha as here 
tdep = 0.3 

sigma_mg = 0.2 # Galactic MMR scatter, in dex [0]
sigma_mc = 0 # Cluster metallicity scatter (within a single galaxy), in dex [0.3]
sigma_gas = 0.3 # Cold gas fraction ratio scatter, in dex

# Turn on stellar mass scatter (only True or False; if True, use the
# Behroozi 2013 has an evolving stellar mass scatter)
sm_scat = True  
gaussian_process = True # whether to evole it with Gaussian process

gauss_l = 2 # correlation length of Gaussian process, in Gyr
gauss_l_sm = 2 # correlation length of Gaussian process for stellar mass, in Gyr

regen_feh = False
seed_feh = 0 # seed for re-generate feh

log_Mhmin = 8.0 # Min halo mass
log_Mmin = 4.0 # Min cluster mass to draw from CIMF

t_lag = 0.01 # time lag for assigning stellar particles, in Gyr

low_mass = False # whether to impose the low mass GCs < Mmin 
log_Mmin_low_mass = 4.0

low_mass_attempt_N = 1 # number of attempts to for clusters if GC Mass is below Mmin

form_nuclear_cluster = True
no_random_at_formation = False 

base_tree = '/nfs/astro2/ybchen/tng50_halos/'
base_halo = '/nfs/astro2/ybchen/tng50_halos/'

# redshift list
redshift_snap = np.loadtxt('/nfs/astro2/ybchen/tng50_halos/TNG_z_list.txt', dtype=float)

# path of mass loss due to stellar evolution
path_massloss = 'data/massloss.txt'

# subhalo list
subs = [523889]

# full snap list if not all snaps are full
full_snap = np.loadtxt('/nfs/astro2/ybchen/tng50_halos/TNG_full_snap.txt', dtype=int)
# full_snap = [2,3]

# Input: color-metallicity transformations to be used for the Virgo
# Cluster GCs. Options right now are "LG14", "CGL18", "V19"
color_metallicity = 'CGL18' 

UVB_constraint = 'KM22'

skip = None

resultspath = '/nfs/astro2/ybchen/temp/'

verbose = True

params = {
    'disrupt_mode':disrupt_mode,
    'h100':h100,
    'Ob':Ob,
    'Om':Om,
    'p2':p2, 
    'p3':p3, 
    'kappa':kappa, 
    'd_tid':d_tid,
    'disrupt_x':disrupt_x,
    'disrupt_y':disrupt_y,
    'log_mc':log_mc, 
    'seed':seed,
    'seed_feh':seed_feh,
    'mpb_only':mpb_only, 
    'mmr_slope':mmr_slope, 
    'mmr_pivot':mmr_pivot,
    'mmr_evolution':mmr_evolution, 
    'max_feh':max_feh,
    'tdep':tdep,
    'sigma_mg':sigma_mg, 
    'sigma_mc':sigma_mc, 
    'sigma_gas':sigma_gas,
    'sm_scat':sm_scat,
    'log_Mhmin':log_Mhmin,
    'log_Mmin':log_Mmin, 
    'pr':pr, 
    't_lag':t_lag,
    'base_tree':base_tree, 
    'base_halo':base_halo,
    'redshift_snap':redshift_snap, 
    'path_massloss':path_massloss,
    'subs':subs,
    'full_snap':full_snap,
    'color_metallicity':color_metallicity, 
    'resultspath':resultspath,
    'allcat_base':allcat_base, 
    'merit_name':merit_name,
    'low_mass':low_mass, 
    'log_Mmin_low_mass':log_Mmin_low_mass,
    'form_nuclear_cluster':form_nuclear_cluster, 
    'low_mass_attempt_N':low_mass_attempt_N,
    'no_random_at_formation':no_random_at_formation, 
    'gaussian_process':gaussian_process,
    'gauss_l': gauss_l,
    'gauss_l_sm': gauss_l_sm,
    'regen_feh':regen_feh,
    'UVB_constraint':UVB_constraint, 
    'verbose':verbose,
    'skip':skip,
    }