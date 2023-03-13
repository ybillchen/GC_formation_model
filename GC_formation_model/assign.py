# Licensed under BSD-3-Clause License - see LICENSE

import time

import numpy as np

from . import loader

def assign(params):

    if params['verbose']:
        print('\n########## assigning model started ##########')

    time_lag = params['t_lag']
    redshift_snap = params['redshift_snap']
    fin_name = params['resultspath'] + params['allcat_name']
    root_name = fin_name[:-4] + '_offset_root.txt'
    offset_name = fin_name[:-4] + '_offset.txt'

    hid, logmh, logms, logmh_form, logms_form, logm_form, z_form, feh, \
        ismpb, hid_form, snapnum_form = np.loadtxt(fin_name, ndmin=2, unpack=True)

    # setting indices to int type
    hid = hid.astype(int)
    ismpb = ismpb.astype(int)
    hid_form = hid_form.astype(int)
    snapnum_form = snapnum_form.astype(int)

    hid_root, idx_root_beg, idx_root_end, idx_beg_in_off, idx_end_in_off = np.loadtxt(
        root_name, ndmin=2, unpack=True, dtype=int)

    snap_form_offset, hid_offset, idx_beg, idx_end = np.loadtxt(
        offset_name, ndmin=2, unpack=True, dtype=int)

    gcid = []
    quality = [] 
    # quality:
    # 2: good (star with t_s > t_form - time_lag):
    # 1: half (star with t_s < t_form - time_lag)
    # 0: bad (dm particle)

    # loop over all subhalos 
    for j in range(len(hid_root)):
        if params['verbose']:
            print(' NO. %d, halo id: %d'%(j,hid_root[j]))

        params['rng'] = np.random.default_rng(params['seed']) # initialize seed

        tree = loader.load_merger_tree(params['base_tree'], hid_root[j])

        # loop over all GC formation events
        for i in range(idx_beg_in_off[j], idx_end_in_off[j]):

            idx_in_tree = np.where( (tree['SnapNum']==snap_form_offset[i]) &
                (tree['SubfindID']==hid_offset[i]) )[0][0]

            num_gc = idx_end[i] - idx_beg[i] # number of GCs formed at the subhalo

            # cosmic time of GC formation
            t_form = params['cosmo'].cosmicTime(
                redshift_snap[snap_form_offset[i]], units = 'Gyr')
            t_before = params['cosmo'].cosmicTime(
                redshift_snap[snap_form_offset[i]-1], units = 'Gyr')
            dt = t_form - t_before
            scale_a = 1 / (1 + redshift_snap[snap_form_offset[i]])

            hpos = tree['SubhaloPos'][idx_in_tree]

            # load cutout snapshot of the subhalo
            fields = ['Coordinates', 'GFM_StellarFormationTime', 'ParticleIDs']
            cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                hid_offset[i], snap_form_offset[i], 'stars', fields)

            if cutout['count'] == 0:
                # if there is no sellar particles in this subhalo, use dm particles
                fields = ['Coordinates', 'ParticleIDs']
                cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                    hid_offset[i], snap_form_offset[i], 'dm', fields)

                pos = cutout['Coordinates']
                dmid = cutout['ParticleIDs'].astype(int)
                x = pos[:,0] - hpos[0]
                y = pos[:,1] - hpos[1]
                z = pos[:,2] - hpos[2]
                r2_dm = x**2 + y**2 + z**2

                # select the closest dm particles to center
                idx_sort_r1_dm = np.argsort(r2_dm)
                dmid = dmid[idx_sort_r1_dm]

                for dm in dmid[:num_gc]:
                    gcid.append(dm)
                    quality.append(0)

                continue

            a_s = cutout['GFM_StellarFormationTime']
            sid = cutout['ParticleIDs'].astype(int)
            pos = cutout['Coordinates']

            x = pos[:,0] - hpos[0]
            y = pos[:,1] - hpos[1]
            z = pos[:,2] - hpos[2]
            r2_s = (x**2 + y**2 + z**2) * (scale_a / params['cosmo'].h)**2 # in kpc

            # select stars with valid scale factor (0-1)
            idx_valid_s = np.where( (a_s>=0) & (a_s<=1) & (r2_s<9))[0] 
            a_s = a_s[idx_valid_s]
            sid = sid[idx_valid_s]

            # cosmic time of stars
            z_s = -1 + 1./a_s # redshift of stars
            t_s = np.array(params['cosmo'].cosmicTime(z_s, units = 'Gyr'))

            # minimum time lag is time_lag
            idx_in_lag_min = np.where((t_s < t_form) & (t_s > t_form - time_lag))[0]

            # maximum time lag is half the time interval between snapshots
            idx_in_lag_max = np.where((t_s < t_form) & 
                (t_s > t_form - max(dt/2, time_lag) ))[0]

            if len(idx_in_lag_max) < num_gc:
                # if there are not enough sellar particles in this subhalo
                # first, fill all stellar particles
                for s, t in zip(sid[idx_in_lag_max], t_s[idx_in_lag_max]):
                    gcid.append(s)
                    if t > t_form - time_lag:
                        quality.append(2)
                    else:
                        quality.append(1)

                # second, use dm particles
                hpos = tree['SubhaloPos'][idx_in_tree]
                fields = ['Coordinates', 'ParticleIDs']
                cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                    hid_offset[i], snap_form_offset[i], 'dm', fields)

                pos = cutout['Coordinates']
                dmid = cutout['ParticleIDs'].astype(int)
                x = pos[:,0] - hpos[0]
                y = pos[:,1] - hpos[1]
                z = pos[:,2] - hpos[2]
                r2_dm = x**2 + y**2 + z**2

                # select the closest dm particles to center
                idx_sort_r1_dm = np.argsort(r2_dm)
                dmid = dmid[idx_sort_r1_dm]

                for dm in dmid[:num_gc-len(idx_in_lag_max)]:
                    gcid.append(dm)
                    quality.append(0)

            elif len(idx_in_lag_min) < num_gc:
                # if there are enough stellar particles in this subhalo, but not 
                # enough in time_lag use the most newly formed stelar particels
                sid = sid[idx_in_lag_max]
                t_s = t_s[idx_in_lag_max]

                # sort by formation time
                idx_sort_t_s = np.argsort(t_s)[::-1]
                sid = sid[idx_sort_t_s]
                t_s = t_s[idx_sort_t_s]

                for s, t in zip(sid[:num_gc], t_s[:num_gc]):
                    gcid.append(s)
                    if t > t_form - time_lag:
                        quality.append(2)
                    else:
                        quality.append(1)

            else:
                # if there are enough stellar particles in time_lag, 
                # just randomly pick
                sid = sid[params['rng'].choice(idx_in_lag_min, num_gc, replace=False)]
                for s in sid:
                    gcid.append(s)
                    quality.append(2)

    output = np.array([gcid, quality], dtype=int).T
    header = 'GC ID'
    np.savetxt(fin_name[:-4]+'_gcid.txt', output, fmt='%d', header=header)

    if params['verbose']:
        print('########## assigning model done ##########')