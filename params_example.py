# Licensed under BSD-3-Clause License - see LICENSE

import numpy as np

# you need to specify all paths starting with $ to your local paths

# disrupt mode: constant, tidal
disrupt_mode = 'tidal'

# adjustable model parameters (run mode)
# Chen & Gnedin (2024a) values for TNG50
p2 = 18
p3 = 0.5
kappa = 1.5

# grid length in calculating tidal tensor, in kpc
d_tid = 0.3 

# x and y parameters if disrupt_mode == 'tidal'
disrupt_x = 0.67
disrupt_y = 1.33

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

# Scaling relations parameters, as in Chen & Gnedin 2023
mmr_slope = 0.3 # Mass slope of galaxy mass-metallicity relation
mmr_pivot = 9.0 # Pivot in the power-law of mass-metallicity relation
mmr_evolution = 1.0 # Redshift slope of the galaxy mass-metallicity relation
mmr0 = -0.5 # The constant to shift MMR
max_feh = 0.3 # Max [Fe/H]

# Scaling of the gas depletion time with redshift, tdep \propto
# (1+z)^(-alpha), alpha as here 
tdep = 0.3 

sigma_mg = 0.3 # Galactic MMR scatter, in dex [0]
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

base_tree = '$base_tree'
base_halo = '$base_halo'
resultspath = '$result'

# redshift list
redshift_snap = np.loadtxt('$z_list.txt', dtype=float)

# path of mass loss due to stellar evolution
path_massloss = 'data/massloss.txt' # not used

mu_sev = 0.55 # fraction of stellar mass after stellar evolution

rmax_form = 3.0 # max radius to form GCs, in kpc

# subhalo list
subs = [523889] # TNG50

# Input: color-metallicity transformations to be used for the Virgo
# Cluster GCs. Options right now are "LG14", "CGL18", "V19"
color_metallicity = 'CGL18' 

UVB_constraint = 'KM22'

# output the label whether the GC mass exceeds stellar mass
exceed_stellar = False

# whether to fix stellar mass when the GC mass exceeds stellar mass
fix_stellar = False

skip = None

verbose = True

full_snap_only = False

max_lag_ratio = 0.5

form_nuclear_cluster = True

# Below: simulation specific parameters

# Cosmological parameters
Om = 0.3089; Ob = 0.0486; h100 = 0.6774 # TNG50

# full snap list if not all snaps are full
full_snap = [2, 3, 4, 6, 8, 11, 13, 17, 21, 25, 33, 40, 50, 59, 67, 72, 78, 84, 91, 99] # TNG50

analyse_snap = 99

params = {
    'disrupt_mode':disrupt_mode,
    'h100':h100,
    'Ob':Ob,
    'Om':Om,
    'p2':p2, 
    'p3':p3, 
    'kappa':kappa, 
    'p2_arr':p2_arr,
    'p3_arr':p3_arr,
    'kappa_arr':kappa_arr,
    'd_tid':d_tid,
    'disrupt_x':disrupt_x,
    'disrupt_y':disrupt_y,
    'log_mc':log_mc, 
    'seed':seed,
    'seed_list':seed_list,
    'seed_feh':seed_feh,
    'mpb_only':mpb_only, 
    'test_mmr':test_mmr,
    'mmr_slope':mmr_slope, 
    'mmr_pivot':mmr_pivot,
    'mmr_evolution':mmr_evolution, 
    'max_feh':max_feh,
    'mmr0':mmr0,
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
    'mu_sev':mu_sev,
    'rmax_form':rmax_form,
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
    'skip':skip,
    'gaussian_process':gaussian_process,
    'gauss_l': gauss_l,
    'gauss_l_sm': gauss_l_sm,
    'regen_feh':regen_feh,
    'UVB_constraint':UVB_constraint, 
    'exceed_stellar':exceed_stellar,
    'fix_stellar':fix_stellar,
    'verbose':verbose,
    'analyse_snap':analyse_snap,
    'full_snap_only':full_snap_only,
    'max_lag_ratio':max_lag_ratio,
    }