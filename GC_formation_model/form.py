# Licensed under BSD-3-Clause License - see LICENSE

# important note about randomness:
# Any slight modification will change the random number generation entirely! 
# To keep repeatability of the model, please construct a new random generator
# for the need of new random numbers

import os
import sys
import time

import numpy as np
from scipy.interpolate import interp1d

from . import astro_utils
from . import schechter_interp
from . import loader

# Define globular cluster class to store data about each cluster
class GC(object): 
    def __init__(self, mass, origin_mh, origin_redshift, metallicity, 
        origin_sm, origin_mgas, is_mpb, idform, snapnum) : 

        self.mass = mass
        self.originHaloMass = origin_mh
        self.origin_redshift = origin_redshift
        self.metallicity = metallicity
        self.origin_sm = origin_sm
        self.origin_mgas = origin_mgas
        self.is_mpb = is_mpb
        self.idform = idform
        self.snapnum = snapnum


# Calculate gas mass given stellar mass, halo mass, redshift using scaling 
# relations. Double power law for SM-Mg relation, then scale with redshift. 
# Revise if gas fraction exceeds accreted baryon fraction. As described in 
# Choksi, Gnedin, and Li (2018).
def gasMass(SM, Mh, z, params) : 
    tdep = params['tdep']; sigma_gas = params['sigma_gas']

    slope = 0.33
    if SM < 1e9:
        slope = 0.19
    log_ratio = 0.05 - 0.5 -  slope*(np.log10(SM) - 9.0) # log10(Mg/M*)
    if z < 3: # Gas fraction saturates at z > 3
        if z < 2: 
            # Strong ssfr evolution at z < 2
            log_ratio += (3.0-tdep)*np.log10((1.+z)/3.) + (3.0-tdep)*np.log10(3.) 
        else:
            # Weak ssfr evolution at z > 2 (Lilly+)
            log_ratio += (1.7-tdep)*np.log10((1.+z)/3.) + (3.0-tdep)*np.log10(3.) 
    else:
        # Weak ssfr evolution at z > 2 
        log_ratio += (1.7-tdep)*np.log10((1.+3)/3.) + (3.0-tdep)*np.log10(3.) 
    if sigma_gas > 0:
        log_ratio += params['rng'].normal(0, sigma_gas)
    ratio = 10**log_ratio
    Mg = SM*ratio
    fstar = SM/(params['cosmo'].fb*Mh)
    fgas = Mg/(params['cosmo'].fb*Mh)

    if params['UVB_constraint'] == 'MG10':
        # Model extragalactic UVB as in Muratov & Gnedin 2010
        Mc = 3.6e9*np.exp(-.6*(1+z))/params['h100']
        Mc_min = 1.5e10*(180**(-.5))/(params['cosmo'].E(z)*params['h100'])
        if Mc < Mc_min : 
            Mc = Mc_min
        fin = 1/((1+(Mc/Mh))**3)

    elif params['UVB_constraint'] == 'KM22':
        # Model extragalactic UVB as in Kravtsov & Manwadkar 2022
        s = lambda x,y: (1 + (2**(y/3)-1)*(x**y))**(-3/y)
        z_re = 6
        gam = 15
        beta = z_re / (np.log(1.82e3*np.exp(-0.63*z_re)-1))**(1/gam)
        Mc = 1.69e10*np.exp(-.63*z)/(1+np.exp(
            np.clip((z/beta)**gam, 0, 500))) # too large number results in overflow
        fin = s(Mc/Mh,2)

    else:
        # no constraint
        return Mg
    
    if fstar+fgas > fin: 
        fgas = fin-fstar
        Mg = fgas*params['cosmo'].fb*Mh

    return Mg


# Galaxy stellar mass-metallicity relation, with additional redshift evolution.
def MMR(SM, z, params) : 

    local = params['mmr_slope']*(np.log10(SM) - params['mmr_pivot']) 
    evolution = params['mmr_evolution']*np.log10(1+z)
    fe_h = local - evolution

    # The old implementation. Now we use Gaussian process
    if not params['gaussian_process']:
        if params['sigma_mg'] > 0:
            # Scatter of galactic MMR
            fe_h += params['rng'].normal(0, params['sigma_mg']) 

    # Metallicity saturates at slightly supersolar values
    if fe_h > params['max_feh']:  
        fe_h = params['max_feh']
    return fe_h 

def clusterFormation(Mg, halomass, redshift, metallicity, SM, is_mpb, subid, 
    mgc_to_mmax, Mmin, ug52, snapnum, params) : 

    gc_list = []

    # Total mass of all GCs formed in cluster formation event
    Mgc = (3e-5/params['cosmo'].fb)*Mg*params['p2'] 

    log_Mgc = np.log10(Mgc)
    if Mgc < Mmin: # Not enough mass to form a single cluster of mass Mmin
        if params['no_random_at_formation']:
            if params['sigma_mc'] > 0:
                cluster_metallicity = metallicity + params['rng'].normal(0, params['sigma_mc'])
            else:
                cluster_metallicity = metallicity
            gc_list.append(GC(Mgc, halomass, redshift, cluster_metallicity, SM, Mg, is_mpb, subid, snapnum))
            return gc_list
        if params['log_mc'] < 10: 
            mt = np.logspace(np.log10(Mmin), params['log_mc'], num = 500)
            ntot = ug52 - schechter_interp.upper_gamma2(params['log_mc'])
            cum = np.array([(ug52 - schechter_interp.upper_gamma2(np.log10(mv)))/
                ntot for mv in mt])
            mt = np.squeeze(mt)
            cum = np.squeeze(cum) 
            r_to_m = interp1d(cum, mt)

        for i in range(params['low_mass_attempt_N']):
            r = params['rng'].random()
            if params['log_mc'] < 10: 
                M = r_to_m(r)
            else: 
                M = Mmin/(1-r*(1.-Mmin/Mmax)) 
            if params['rng'].random() < Mgc / M / params['low_mass_attempt_N']:
                if params['sigma_mc'] > 0:
                    cluster_metallicity = metallicity + params['rng'].normal(0, params['sigma_mc'])
                else:
                    cluster_metallicity = metallicity
                gc_list.append(GC(M, halomass, redshift, cluster_metallicity, SM, Mg, is_mpb, subid, snapnum))
        return gc_list

    if params['form_nuclear_cluster']:
        # Optimal sampling: form the nuclear cluster 
        log_Mmax = mgc_to_mmax(log_Mgc)
        if log_Mmax > log_Mgc:
            log_Mmax = log_Mgc

        Mmax = 10**log_Mmax
        if params['sigma_mc'] > 0:
            cluster_metallicity = metallicity + params['rng'].normal(0, params['sigma_mc'])
        else:
            cluster_metallicity = metallicity
        gc_list.append(GC(Mmax, halomass, redshift,cluster_metallicity, SM, Mg, 
            is_mpb, subid, snapnum))
        mass_sum = Mmax

    else:
        log_Mmax = params['log_mc']
        mass_sum = 0

    # Calculate the cumulative distribution r(<M), and invert it numerically
    # If CIMF cutoff mass is below 10, then calcualte r(M) numerically, using 
    # the full expression including gamma functions
    # If CIMF cutoff mass is above 10, then we'll just use the analytic solution 
    # for r(M) in the power-law limit. this speeds up computation 
    # *significantly*! (factor of >5 or so)
    if params['log_mc'] < 10: 
        mt = np.logspace(np.log10(Mmin), log_Mmax, num = 500)
        ntot = ug52 - schechter_interp.upper_gamma2(log_Mmax)
        cum = np.array([(ug52 - schechter_interp.upper_gamma2(np.log10(mv)))/
            ntot for mv in mt])
        mt = np.squeeze(mt)
        cum = np.squeeze(cum) 
        r_to_m = interp1d(cum, mt)

    while mass_sum < Mgc:
        r = params['rng'].random()
        if params['log_mc'] < 10: 
            M = r_to_m(r)
        else: 
            M = Mmin/(1-r*(1.-Mmin/Mmax)) 

        ###### old method: 
        # Make sure the final cluster drawn doesn't exceed the total mass to be 
        # formed. it may produce some clusters below Mmin, but shouldn't really 
        # matter (will disrupt)
        ###### new method: 
        # keep the last gc with prob. = (Mgc - mass_sum) / M
        if mass_sum + M > Mgc: 
            # M = Mgc-mass_sum # old method
            if params['rng'].random() > (Mgc - mass_sum) / M:
                break
        mass_sum += M 
        if params['sigma_mc'] > 0:
            cluster_metallicity = metallicity + params['rng'].normal(0, params['sigma_mc'])
        gc_list.append(GC(M, halomass, redshift, cluster_metallicity, SM, Mg, is_mpb, subid, snapnum))


    # low mass GCs
    if params['low_mass']:
        DMgc = (Mgc - Mmax) * \
            (schechter_interp.upper_gamma1(params['log_Mmin_low_mass']) - 
                schechter_interp.upper_gamma1(params['log_Mmin'])) / \
            (schechter_interp.upper_gamma1(params['log_Mmin']) - 
                schechter_interp.upper_gamma1(log_Mmax))

        mt = np.logspace(params['log_Mmin_low_mass'], np.log10(Mmin), num = 500)
        ntot = (schechter_interp.upper_gamma2(params['log_Mmin_low_mass']) -
            schechter_interp.upper_gamma2(np.log10(Mmin)))
        cum = np.array([(schechter_interp.upper_gamma2(params['log_Mmin_low_mass']) - 
            schechter_interp.upper_gamma2(np.log10(mv)))/ntot for mv in mt])
        mt = np.squeeze(mt)
        cum = np.squeeze(cum) 
        r_to_m = interp1d(cum, mt)


        mass_sum = 0
        while mass_sum < DMgc:
            r = params['rng'].random()
            M = r_to_m(r)

            # Make sure the final cluster drawn doesn't exceed the total mass to be 
            # formed. it may produce some clusters below Mmin, but shouldn't really 
            # matter (will disrupt)
            if mass_sum+M > DMgc: 
                # M = DMgc-mass_sum
                if params['rng'].random() > (DMgc - mass_sum) / M:
                    break
            mass_sum += M 
            if params['sigma_mc'] > 0:
                cluster_metallicity = metallicity + params['rng'].normal(0, params['sigma_mc'])
            gc_list.append(GC(M, halomass, redshift, cluster_metallicity, SM, Mg, is_mpb, subid, snapnum))

    return gc_list

def gaussian_process(rng, thist, sig, l=2):
    
    t1, t2 = np.meshgrid(thist,thist)
    cov = sig**2 * np.exp(-(t1-t2)**2 / 2 / l**2)

    return rng.multivariate_normal(mean=np.zeros(len(thist)), cov=cov)

def gaussian_process_sm(rng, thist, zhist, l=2):
    
    t1, t2 = np.meshgrid(thist,thist)
    z1, z2 = np.meshgrid(zhist,zhist)
    var = (0.218+0.023*z1/(1+z1)) * (0.218+0.023*z2/(1+z2))
    cov = var * np.exp(-(t1-t2)**2 / 2 / l**2)

    return rng.multivariate_normal(mean=np.zeros(len(thist)), cov=cov)

# organize tree in a convenient order
def organize_tree(tree, params):

    # update tree with new index
    def update(idx, old_tree):
        [m, fp, sp, descid, subid, snapnum, mpi, subfindid] = old_tree
        return np.array([m[idx], fp[idx], sp[idx], descid[idx], subid[idx], 
            snapnum[idx], mpi[idx], subfindid[idx]])

    # remove spikes in mass
    def despike(arr, n=1):
        data = np.copy(arr) # in Msun
        for i in range(n, len(arr)-n):
            if np.max(arr[i-n:i+n+1]) / np.min(arr[i-n:i+n+1]) > 10 and np.max(arr[i-n:i+n+1]) > 1e10:
                data[i] = 0
        return data

    m = tree['SubhaloMass'] * 1e10 / params['cosmo'].h # 0
    fp = tree['FirstProgenitorID'] # 1
    sp = tree['NextProgenitorID'] # 2
    descid = tree['DescendantID'] # 3
    subid = tree['SubhaloID'] # 4
    snapnum = tree['SnapNum'] # 5
    mpi = tree['MainLeafProgenitorID'] # 6
    subfindid = tree['SubfindID'] # 7
    old_tree = np.array([m, fp, sp, descid, subid, snapnum, mpi, subfindid], dtype=int)

    idx_valid = (m > 10**params['log_Mhmin']) & (fp != -1)
    new_tree = update(idx_valid, old_tree)

    # new_tree[0] = despike(new_tree[0]) # remove spikes in mass

    snapnum_uique = np.unique(new_tree[5])[::-1]
    halos2 = []
    dfeh2 = []
    dsm2 = []

    for mpiv in np.unique(new_tree[6]):

        # get all the halos along this branch
        idx_mpi = np.where(new_tree[6] == mpiv)[0] 
        mpi_tree = update(idx_mpi, new_tree)

        #sort in order of ascending snapnum
        idx_snapnum = np.argsort(mpi_tree[5])
        mpi_tree = update(idx_snapnum, mpi_tree)
        
        mpi_tree[0] = despike(mpi_tree[0]) # remove spikes in mass

        zlist = [params['redshift_snap'][s] for s in mpi_tree[5]]
        tlist = [params['cosmo'].cosmicTime(z, units = 'Gyr') for z in zlist]
        dfeh = gaussian_process(params['rng'], tlist, params['sigma_mg'], l=params['gauss_l'])
        if params['regen_feh']:
            dfeh = gaussian_process(params['rng_feh'], tlist, params['sigma_mg'], l=params['gauss_l'])
        dsm = gaussian_process_sm(params['rng'], tlist, zlist, l=params['gauss_l_sm'])

        halos2.append(mpi_tree[:,0])
        dfeh2.append(dfeh[0])
        dsm2.append(dsm[0])
        snapnum_uique = np.unique(mpi_tree[5])[::-1]

        i = 0
        while i <= len(snapnum_uique) - 1:
            snapnum_now = snapnum_uique[i]
            halo_now = mpi_tree[:,i]
            j = i+1
            found = False

            # search for the next halo along this branch that has delta M > 0 
            # relative to halo_now (i.e., no mass decrease)
            while not found and j < len(snapnum_uique): 
                snapnum_next = snapnum_uique[j]
                desc = mpi_tree[:,j]
                dfeh_desc = dfeh[j]
                dsm_desc = dsm[j]
                if desc[0] < halo_now[0]:
                    j += 1
                else:
                    halos2.append(np.array([desc[0], halo_now[4], desc[2], desc[3],
                        desc[4], desc[5], desc[6], desc[7]], dtype=int)) #set fp of this halo to h_id of halo_now
                    dfeh2.append(dfeh_desc)
                    dsm2.append(dsm_desc)
                    found = True
            i = j

    return np.array(halos2, dtype=int).T, np.array(dfeh2), np.array(dsm2)

# Main calculation: loop over all halos in the merger tree
def form(params):

    if params['verbose']:
        print('\n########## formation model started ##########')

    Mmin = 10**params['log_Mmin']

    schechter_interp.init(10**params['log_mc'])
    mgc_to_mmax = schechter_interp.generate(10**params['log_mc'], mmin = Mmin)
    ug52 = schechter_interp.upper_gamma2(params['log_Mmin'])

    for num_run, hid_num in enumerate(params['subs']):
        if params['verbose']:
            print(' NO. %d, halo id: %d'%(num_run,hid_num))

        params['rng'] = np.random.default_rng(params['seed']+hid_num) # initialize seed
        if params['regen_feh']:
            params['rng_feh'] = np.random.default_rng(params['seed_feh']+hid_num+1) # initialize seed for feh

        t0 = time.time()

        tree = loader.load_merger_tree(params['base_tree'], hid_num)

        t1 = time.time() # t1 - t0 is time for loading tree
        if params['verbose']:
            print('  - load tree: %.2f s'%(t1-t0))

        organized_tree, dfeh, dsm = organize_tree(tree, params)
        if len(organized_tree) < 1:
            continue
        m, fp, sp, descid, subid, snapnum, mpi, subfindid = organized_tree

        t2 = time.time() # t2 - t1 is time for re-organizing tree
        if params['verbose']:
            print('  - organize tree: %.2f s'%(t2-t1))

        msub = np.max(m)
        mpbi = mpi[m == np.max(m)][0]

        # Go through each halo along the tree and look for events satisfying Rm > p3.
        sm_arr = np.zeros(len(m))
        gal_feh_arr = np.zeros(len(m))
        clusters = []
        exceed_stellar_label = []
        for i in range(len(m)) : # For each halo in the merger tree
            mass = m[i] # Mass of this halo
            fpID = fp[i] # ID of the main progenitor
            znow = params['redshift_snap'][snapnum[i]] # Current redshift
            sm1 = astro_utils.SMHM(mass, znow, scatter = False) # mean sm

            if fpID == -1 or len(subid[subid == fpID]) == 0:
                # Then we've reached the first point along this track of the tree 
                if params['gaussian_process']:
                    if params['sm_scat']:
                        sm_arr[i] = sm1*(10**dsm[i])
                    else:
                        sm_arr[i] = sm1
                else:
                    # Assign a "seed" stellar mass which we will grow self-consistently
                    sm_arr[i] = astro_utils.SMHM(mass, znow, 
                        scatter = params['sm_scat'])
                continue

            # Identify index of first progenitor in data 
            progIdx = np.where(subid == fpID)[0][0]
            progMass = m[progIdx] # Get mass of fprogenitor  

            ratio = mass/progMass - 1 # Calculate merger ratio, Rm = dMh/Mh

            zbefore = params['redshift_snap'][snapnum[progIdx]]

            dt = ( params['cosmo'].cosmicTime(znow, units = 'Gyr') - 
                params['cosmo'].cosmicTime(zbefore, units = 'Gyr') ) 
            ratio = ratio/dt # This is just (dMh/Mh)/dt

            if params['gaussian_process']:
                # Evolve sm using Gaussian process
                sm1 = astro_utils.SMHM(mass, znow, scatter = False)
                if params['sm_scat']:
                    SM = np.max([sm1*(10**dsm[i]), sm_arr[progIdx]])
                else:
                    SM = sm1
                sm_arr[i] = SM

                # Evolve feh using Gaussian process
                gal_feh1 = MMR(SM, znow, params)
                gal_feh = gal_feh1 + dfeh[i]
                if gal_feh > params['max_feh']:  
                    gal_feh = params['max_feh']
                gal_feh_arr[i] = gal_feh

            else:
                # Evolve stellar mass self-consistently as described in CGL18
                sm1 = astro_utils.SMHM(mass, znow, scatter = False)
                sm2 = astro_utils.SMHM(progMass, zbefore, scatter = False)
                dsm = sm1-sm2
                a = 1./(1+znow)

                if params['sm_scat']:
                    dsm *= 10**params['rng'].normal(0,.218 - .023*(a-1)) 
                SM  = sm_arr[progIdx] + dsm

                if SM < 0: 
                    # Only happens in a couple very weird cases at very high redshift
                    SM = sm_arr[progIdx]
                sm_arr[i] = SM

            Mg = gasMass(SM, mass, znow, params) 
            if (ratio > params['p3']) and (Mg > 0):
                if params['gaussian_process']:
                    galaxy_metallicity = gal_feh
                else:
                    galaxy_metallicity = MMR(SM, znow, params) 
                is_mpb = mpi[i] == mpbi
                new_clusters = clusterFormation(Mg, mass, znow, 
                    galaxy_metallicity, SM, is_mpb, subfindid[i], mgc_to_mmax, 
                    Mmin, ug52, snapnum[i], params)
                clusters.extend(new_clusters)
                new_clusters_mass_tot = np.sum([cluster.mass for cluster in new_clusters])
                if new_clusters_mass_tot > sm_arr[i] - sm_arr[progIdx]:
                    # more GC mass than stellar mass, give a label
                    exceed_stellar_label.extend([1] * len(new_clusters))
                else:
                    exceed_stellar_label.extend([0] * len(new_clusters))

        # All GCs that form, regardless of survival -- for use w/ allcat.txt
        GC_mets = np.array([cluster.metallicity for cluster in clusters])
        GC_log_masses = np.log10(np.array(
            [cluster.mass for cluster in clusters]))
        GC_redshifts = np.array(
            [cluster.origin_redshift for cluster in clusters])
        GC_mhost_tform = np.array(
            [cluster.originHaloMass for cluster in clusters])
        GC_log_mhost_tform = np.log10(GC_mhost_tform)
        GC_log_mstar_tform = np.log10(np.array(
            [cluster.origin_sm for cluster in clusters]))
        GC_log_mgas_tform = np.log10(np.array(
            [cluster.origin_mgas for cluster in clusters]))
        GC_ismpb = np.array([cluster.is_mpb for cluster in clusters]).astype(int)
        GC_idform = np.array([cluster.idform for cluster in clusters])
        GC_snapnum = np.array([cluster.snapnum for cluster in clusters])

        logmsub = np.log10(msub)
        # logms = np.log10(astro_utils.SMHM(10**logmsub, 0.0, scatter = False))
        logms = np.log10(np.max(sm_arr))*np.ones(len(clusters))
        # GC_log_mstar_tform = np.log10(sm_arr[m == np.max(m)][0])*np.ones(len(clusters))
        
        #hid_num = int(fname[0:fname.find(".")])
        # hid_num = sub['id']

        output = [hid_num*np.ones(len(clusters)),logmsub*np.ones(len(clusters)),
            logms*np.ones(len(clusters)),GC_log_mhost_tform,GC_log_mstar_tform,
            GC_log_masses,GC_redshifts,GC_mets,GC_ismpb, GC_idform,GC_snapnum]

        output = np.array(output).T

        if num_run == 0:
            save_output = output
        else:
            save_output = np.concatenate((save_output,output))
        
        t3 = time.time() # t3 - t2 is time for modeling gcs 
        if params['verbose']:
            print('  - model gc: %.2f s'%(t3-t2))
            print('  total time: %.2f s'%(t3-t0))
            print('  number of GCs:', len(clusters))

    header = ('subfindID(z=0) | logMh(z=0) | logM*(z=0) | logMh(zform) | logM*(zform)' +
        ' | logM(tform) | zform | feh | isMPB | subfindID(zfrom) | snapnum(zform) \n')

    if params['regen_feh']:
        save_path = params['resultspath']+params['allcat_name'][:-4]+'_regen_feh_s-%d_sigg-%g_l-%g.txt'%(
            params['seed_feh'], params['sigma_mg'], params['gauss_l'])
    else:
        save_path = params['resultspath']+params['allcat_name']

    np.savetxt(save_path, save_output, header=header, 
        fmt='%d %6.3f %6.3f %6.3f %6.3f %6.3f %5.3f %6.3f %d %d %d')

    if 'exceed_stellar' in params and params['exceed_stellar']:
        np.savetxt(save_path[:-4]+'_exceed_stellar.txt', exceed_stellar_label, fmt='%d')

    if params['verbose']:
        print('########## formation model done ##########')
