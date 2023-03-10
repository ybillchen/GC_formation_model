# Licensed under BSD-3-Clause License - see LICENSE

import numpy as np

# mode: run or calibrate
mode = 'run'

# disrupt mode: constant, density, tidal, tidal_centrifugal, Mark_Gieles
disrupt_mode = 'Mark_Gieles'

# cosmological parameters
h100 = 0.6774
Ob = 0.0486
Om = 0.3089

# adjustable model parameters (run mode)
p2 = 8.8
p3 = 0.58
kappa = 9.0

# grid length in calculating tidal tensor, in kpc
d_tid = 0.3 

# adjustable model parameters (calibrate mode) Below, we'll loop
# over all combinations of p2_arr and p3_arr: if you want to just run
# a single set of model params, then just make each array a list of
# length 1
p2_arr = [10, 12, 14, 16, 18, 20] # 6 p2
p3_arr = [0.5, 0.6, 0.7] # 3 p3
kappa_arr = [1, 2, 2.5, 3, 3.5, 4] # more kappa

# log(Mcut) - knee of the Schechter function. Anything with log_mc >
# 10 will be treated as a pure power-law, to speed up runtime
# (analytic solution).
log_mc = 7.0 

# Various things you can play with (solely for ease of use, doesn't
# affect model calculations at all)
run_all = True
log_mh_min, log_mh_max = 6., 16.  
##MW: 11.86, 12.38
N = 9999
seed = 0

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

log_Mmin = 4.0 # Min cluster mass to draw from CIMF
pr = 0.5 # Normalized period of rotation for disruption; t_tid \propto P

t_lag = 0.01 # time lag for assigning stellar particles, in Gyr

low_mass = False # whether to impose the low mass GCs < Mmin 
log_Mmin_low_mass = 4.0

low_mass_attempt_N = 1 # number of attempts to for clusters if GC Mass is below Mmin

form_nuclear_cluster = True
no_random_at_formation = False 

# base path of TNG50-1
base = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L35n2160TNG/output/' 
base_tree = '/n/holyscratch01/vogelsberger/billchen/temp/'
base_halo = '/n/holyscratch01/vogelsberger/billchen/temp/'

# redshift list
redshift_snap = np.loadtxt('/n/holyscratch01/vogelsberger/billchen/temp/TNG_z_list.txt', dtype=float) 

# subhalo list
subs = [523889]

# full snap list if not all snaps are full
full_snap = [2,3,4,6,8,11,13,17,21,25,33,40,50,59,67,72,78,84,91,99]

# Input: color-metallicity transformations to be used for the Virgo
# Cluster GCs. Options right now are "LG14", "CGL18", "V19"
color_metallicity = 'CGL18' 

UVB_constraint = 'KM22'

resultspath = '/n/holyscratch01/vogelsberger/billchen/temp/'

verbose = True

params = {
    'mode':mode, 
    'disrupt_mode':disrupt_mode,
    'h100':h100,
    'Ob':Ob,
    'Om':Om,
    'p2':p2, 
    'p3':p3, 
    'kappa':kappa, 
    'd_tid':d_tid,
    'p2_arr':p2_arr, 
    'p3_arr':p3_arr, 
    'kappa_arr':kappa_arr,
    'log_mc':log_mc, 
    'run_all':run_all, 
    'log_mh_min':log_mh_min, 
    'log_mh_max':log_mh_max, 
    'N':N, 
    'seed':seed,
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
    'log_Mmin':log_Mmin, 
    'pr':pr, 
    't_lag':t_lag,
    'base':base, 
    'base_tree':base_tree, 
    'base_halo':base_halo,
    'redshift_snap':redshift_snap, 
    'subs':subs,
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
    'UVB_constraint':UVB_constraint, 
    'verbose':verbose,
    }