# Licensed under BSD-3-Clause License - see LICENSE

import time
import warnings

import numpy as np
from scipy.interpolate import interp1d, LinearNDInterpolator
import scipy.spatial as sp

from . import astro_utils
from . import loader

def calc_eig(tree, pos_gc, pot_gc, pos, pot, d_tid):
    # phi in km^2/s^2
    # pos and d_tid in kpc/h

    grid = np.zeros([27,len(pos_gc),3])

    # grid[0]  = pos_gc + d_tid * np.array([-1,-1,-1])
    grid[1]  = pos_gc + d_tid * np.array([-1,-1, 0])
    # grid[2]  = pos_gc + d_tid * np.array([-1,-1, 1])
    grid[3]  = pos_gc + d_tid * np.array([-1, 0,-1])
    grid[4]  = pos_gc + d_tid * np.array([-1, 0, 0])
    grid[5]  = pos_gc + d_tid * np.array([-1, 0, 1])
    # grid[6]  = pos_gc + d_tid * np.array([-1, 1,-1])
    grid[7]  = pos_gc + d_tid * np.array([-1, 1, 0])
    # grid[8]  = pos_gc + d_tid * np.array([-1, 1, 1])
    grid[9]  = pos_gc + d_tid * np.array([ 0,-1,-1])
    grid[10] = pos_gc + d_tid * np.array([ 0,-1, 0])
    grid[11] = pos_gc + d_tid * np.array([ 0,-1, 1])
    grid[12] = pos_gc + d_tid * np.array([ 0, 0,-1])
    grid[13] = pos_gc + d_tid * np.array([ 0, 0, 0])
    grid[14] = pos_gc + d_tid * np.array([ 0, 0, 1])
    grid[15] = pos_gc + d_tid * np.array([ 0, 1,-1])
    grid[16] = pos_gc + d_tid * np.array([ 0, 1, 0])
    grid[17] = pos_gc + d_tid * np.array([ 0, 1, 1])
    # grid[18] = pos_gc + d_tid * np.array([ 1,-1,-1])
    grid[19] = pos_gc + d_tid * np.array([ 1,-1, 0])
    # grid[20] = pos_gc + d_tid * np.array([ 1,-1, 1])
    grid[21] = pos_gc + d_tid * np.array([ 1, 0,-1])
    grid[22] = pos_gc + d_tid * np.array([ 1, 0, 0])
    grid[23] = pos_gc + d_tid * np.array([ 1, 0, 1])
    # grid[24] = pos_gc + d_tid * np.array([ 1, 1,-1])
    grid[25] = pos_gc + d_tid * np.array([ 1, 1, 0])
    # grid[26] = pos_gc + d_tid * np.array([ 1, 1, 1])

    pot_grid = np.zeros([27,len(pos_gc)])
    for i in range(27):
        if i == 13:
            pot_grid[i] = pot_gc
            continue
        if (i == 0) or (i == 2) or (i == 6) or (i == 8) or (i == 18) or (i == 20) or (i == 24) or (i == 26):
            continue

        idxs = tree.query(grid[i], k=8)[1]
        for j in range(len(idxs)):
            idx = idxs[j]
            phi = LinearNDInterpolator(pos[idx], pot[idx])

            pot_grid[i][j] = phi(grid[i][j])

    for i in range(13):
        idx_nan_1 = np.where(np.isnan(pot_grid[i]))[0]
        idx_nan_2 = np.where(np.isnan(pot_grid[26-i]))[0]
        if (not len(idx_nan_1)) and (not len(idx_nan_2)):
            continue
        idx_0 = np.intersect1d(idx_nan_1, idx_nan_2)
        idx_1 = np.setdiff1d(idx_nan_1, idx_0)
        idx_2 = np.setdiff1d(idx_nan_2, idx_0)

        pot_grid[i][idx_0] = pot_gc[idx_0]
        pot_grid[26-i][idx_0] = pot_gc[idx_0]
        pot_grid[i][idx_1] = 2*pot_gc[idx_1] - pot_grid[26-i][idx_1]
        pot_grid[26-i][idx_2] = 2*pot_gc[idx_2] - pot_grid[i][idx_2]

    # tidal tensor
    Txx = (pot_grid[4] + pot_grid[22] - 2*pot_grid[13]) / d_tid**2
    Tyy = (pot_grid[10] + pot_grid[16] - 2*pot_grid[13]) / d_tid**2
    Tzz = (pot_grid[12] + pot_grid[14] - 2*pot_grid[13]) / d_tid**2
    Txy = (pot_grid[1] + pot_grid[25] - pot_grid[7] - pot_grid[19]) / 4 / d_tid**2
    Txz = (pot_grid[3] + pot_grid[23] - pot_grid[5] - pot_grid[21]) / 4 / d_tid**2
    Tyz = (pot_grid[9] + pot_grid[17] - pot_grid[11] - pot_grid[15]) / 4 / d_tid**2

    T = np.array([
        [Txx, Txy, Txz],
        [Txy, Tyy, Tyz],
        [Txz, Tyz, Tzz]]).T

    return 0.48 * np.linalg.eig(T)[0] # in Gyr^-2

def calc_eig_new(tree, pos_gc, pos, mass):
    idxs = tree.query(pos_gc, k=64)[1]
    pos0 = np.transpose(pos[idxs], axes=(1,0,2)) - pos_gc
    r0 = np.sqrt(np.sum(pos0**2, axis=2))
    Txx = 2*pos0[:,0]**2 - pos0[:,1]**2 - pos0[:,2]**2
    Tyy = 2*pos0[:,1]**2 - pos0[:,2]**2 - pos0[:,0]**2
    Tzz = 2*pos0[:,2]**2 - pos0[:,0]**2 - pos0[:,1]**2
    Txy = 3*pos0[:,0]*pos0[:,1]
    Txz = 3*pos0[:,0]*pos0[:,2]
    Tyz = 3*pos0[:,1]*pos0[:,2]

    T = np.sum(np.array([
        [Txx, Txy, Txz],
        [Txy, Tyy, Tyz],
        [Txz, Tyz, Tzz]])*mass[idxs]/r0**5, axis=2)

    return 20642 * np.linalg.eig(T)[0] # in Gyr^-2

# get tidal tensor for one galaxy
def get_tid_unit(i, gcid, hid_root, idx_beg, idx_end, params, k=-1):
    mpb_only = params['mpb_only']
    d_tid = params['d_tid'] * params['h100'] # in kpc/h
    z_list = params['redshift_snap']
    full_snap = params['full_snap']

    t0 = time.time()

    # load merger tree
    fields = ['SnapNum', 'SubfindID', 'SubhaloMass']
    tree = loader.load_merger_tree(params['base_tree'], hid_root[i], fields)

    # existing GCs at this snapshot
    idx_exist_gc = np.arange(idx_beg[i], idx_end[i])

    tag = np.zeros([len(idx_exist_gc), len(full_snap)], dtype=int)
    eig_1 = np.zeros([len(idx_exist_gc), len(full_snap)])
    eig_2 = np.zeros([len(idx_exist_gc), len(full_snap)])
    eig_3 = np.zeros([len(idx_exist_gc), len(full_snap)])

    if k < 0:
        iterlist = range(len(full_snap))
    else:
        iterlist = [k]

    for j in iterlist:

        if not params['skip'] is None:
            if params['skip'][0] == i and params['skip'][1] == j:
                continue

        t1 = time.time()
        t2 = 0 # load halo
        t3 = 0 # build tree
        t4 = 0 # calc eig

        snap = full_snap[j]
        scale_a = 1 / (1 + z_list[snap])

        # all subhalos at this snapshot
        idx_sub = np.where((tree['SnapNum']==snap) & 
            (tree['SubhaloMass'] > 10**(params['log_Mhmin']-10)*params['h100']))[0]
        subfindid = tree['SubfindID'][idx_sub]

        count = 0 # found GCs
        # loop over all subhalos and load densities
        for subid in subfindid:
            t22 = time.time()
            # first, consider all GCs represented by dm
            fields = ['Coordinates', 'ParticleIDs', 'Potential']
            cutout = loader.load_halo(params['base_halo'], hid_root[i], 
                subid, snap, 'dm', fields)
            pos = cutout['Coordinates'] * scale_a # in kpc/h
            pid = cutout['ParticleIDs'].astype(int)
            pot = cutout['Potential'] / scale_a # in km^2/s^2

            # second, consider all GCs represented by stars
            fields = ['Coordinates', 'ParticleIDs', 'Potential']
            cutout = loader.load_halo(params['base_halo'], hid_root[i], 
                subid, snap, 'stars', fields)
            if cutout['count'] > 0:
                pos = np.concatenate((pos, cutout['Coordinates'] * scale_a)) # in kpc/h
                pid = np.concatenate((pid, cutout['ParticleIDs'].astype(int)))
                pot = np.concatenate((pot, cutout['Potential'] / scale_a)) # in km^2/s^2

            # find intersections, xy is useless
            xy, idx_1, idx_2 = np.intersect1d(pid, 
                gcid[idx_exist_gc], return_indices=True)

            pos_gc = pos[idx_1]
            pot_gc = pot[idx_1]

            if len(xy) == 0: # if gc particles not found 
                continue

            count += len(xy)

            fields = ['Coordinates', 'Potential']
            cutout = loader.load_halo(params['base_halo'], hid_root[i], 
                subid, snap, 'gas', fields)
            if cutout['count'] > 0:
                pos = np.concatenate((pos, cutout['Coordinates'] * scale_a)) # in kpc
                pot = np.concatenate((pot, cutout['Potential'] / scale_a)) # in km^2/s^2
            
            t33 = time.time()
            t2 += (t33 - t22) # load halo

            kdtree = sp.KDTree(pos)

            t44 = time.time()
            t3 += (t44 - t33) # build tree

            try:
                eig = calc_eig(kdtree, pos_gc, pot_gc, pos, pot, d_tid) # in Gyr^-2
            except sp.qhull.QhullError as e:
                # QhullError may occur. No perfect solution yet. Here we skip this subhalo
                t4 += time.time() - t44 # calc eig
                warnings.warn('A QhullError happens at NO. %d snap %d)!'%(i,full_snap[j]))
                continue

            # update the tag and eig matrices
            tag[idx_2,j] = np.ones(len(xy), dtype=int)

            eig_sorted = np.sort(eig, axis=1)

            eig_1[idx_2,j] = eig_sorted[:,2]
            eig_2[idx_2,j] = eig_sorted[:,1]
            eig_3[idx_2,j] = eig_sorted[:,0]

            t4 += time.time() - t44 # calc eig

        if params['verbose']:
            print('  * NO. %d, hid: %d, snap: %d, %d/%d found, time: %.1fs'%(
                i, hid_root[i], snap, count, len(idx_exist_gc), time.time()-t1))
            print('   - load halo: %.1fs, build tree: %.1fs, calc eig: %.1fs'%(t2,t3,t4))

    if params['verbose'] and k < 0:
        print(' NO. %d, halo id: %d completed, total time: %.1f s'%(i, hid_root[i], time.time()-t0))

    if k < 0:
        return tag, eig_1, eig_2, eig_3

    return tag[:,j], eig_1[:,j], eig_2[:,j], eig_3[:,j]

def get_tid(params):

    if params['verbose']:
        print('\n########## tidal tensor calculation started ##########')

    time_lag = params['t_lag']
    redshift_snap = params['redshift_snap']
    full_snap = params['full_snap']
    fin_name = params['resultspath'] + params['allcat_name']
    gcid_name = fin_name[:-4] + '_gcid.txt'
    root_name = fin_name[:-4] + '_offset_root.txt'
    offset_name = fin_name[:-4] + '_offset.txt'

    # load GC id
    gcid, quality = np.loadtxt(gcid_name, ndmin=2, unpack=True, dtype='int64')

    # load root offset
    hid_root, idx_beg, idx_end = np.loadtxt(
        root_name, ndmin=2, unpack=True, dtype='int64')[:3]

    tag = np.zeros([len(gcid), len(full_snap)], dtype=int)
    eig_max = np.zeros([len(gcid), len(full_snap)])
    eig_1 = np.zeros([len(gcid), len(full_snap)])
    eig_2 = np.zeros([len(gcid), len(full_snap)])
    eig_3 = np.zeros([len(gcid), len(full_snap)])

    for i in range(len(hid_root)):
        if params['verbose']:
            print(' NO. %d, halo id: %d'%(i,hid_root[i]))

        tag_i, eig_1_i, eig_2_i, eig_3_i = get_tid_unit(i, gcid, hid_root, idx_beg, idx_end, params)

        tag[idx_beg[i]:idx_end[i]] = tag_i
        eig_1[idx_beg[i]:idx_end[i]] = eig_1_i
        eig_2[idx_beg[i]:idx_end[i]] = eig_2_i
        eig_3[idx_beg[i]:idx_end[i]] = eig_3_i

    # save data
    np.savetxt(fin_name[:-4]+'_tidtag.txt', tag, fmt='%d')
    np.savetxt(fin_name[:-4]+'_tideig1.txt', eig_1, fmt='%.3e')
    np.savetxt(fin_name[:-4]+'_tideig2.txt', eig_2, fmt='%.3e')
    np.savetxt(fin_name[:-4]+'_tideig3.txt', eig_3, fmt='%.3e')

    if params['verbose']:
        print('########## tidal tensor calculation done ##########')