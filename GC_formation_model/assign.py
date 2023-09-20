# Licensed under BSD-3-Clause License - see LICENSE

import time

import numpy as np
import scipy.spatial as sp

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

    if 'full_snap_only' in params and params['full_snap_only']:
        offset_full_snap_name = offset_name[:-4] + '_full_snap.txt'
        snap_form_offset_full_snap, hid_offset_full_snap, idx_beg, idx_end = np.loadtxt(
            offset_full_snap_name, ndmin=2, unpack=True, dtype=int)

    gcid = []
    quality = [] 
    # quality: 
    # if params['collisionless_only'] == False
    # 2: good (star with t_s > t_form - time_lag):
    # 1: half (star with t_s < t_form - time_lag)
    # 0: bad (dm particle)

    # if params['collisionless_only'] == True and params['assign_at_peaks'] == True
    # 3: formed in peaks
    # 2: formed in peak surroundings
    # 1: formed in scale radius
    # 0: outside scale radius

    # if params['collisionless_only'] == True and params['assign_at_peaks'] == False
    # 0: all cases

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

            if 'collisionless_only' in params and params['collisionless_only']:
                rs = tree['ScaleRad'][idx_in_tree] / params['h100'] # in kpc

                fields = ['Coordinates', 'ParticleIDs', 'Masses']
                cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                    hid_offset[i], snap_form_offset[i], 'dm', fields)

                pos = cutout['Coordinates']
                dmid = cutout['ParticleIDs'].astype(int)
                x = pos[:,0] - hpos[0]
                y = pos[:,1] - hpos[1]
                z = pos[:,2] - hpos[2]
                r2_dm = x**2 + y**2 + z**2

                # select the closest dm particles to center
                idx_sort_r2_dm = np.argsort(r2_dm)
                dmid = dmid[idx_sort_r2_dm]
                r2_dm = r2_dm[idx_sort_r2_dm]
                mass = cutout['Masses'][idx_sort_r2_dm] # in Msun

                # select dm particles in rs*frac_rc
                idx_in = np.where(r2_dm < (rs*params['frac_rs'])**2)[0]

                if params['assign_at_peaks']:
                    # The following are mass densities. In ART simulations,
                    # densities are stored in cells (gas). So, we keep the 
                    # name 'gas' here. If you store densities in particles,
                    # like Gadget, the line below should load the same cutout
                    # as 'dm'.
                    fields = ['Coordinates', 'Density']
                    cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                    hid_offset[i], snap_form_offset[i], 'gas', fields)
                    pos = cutout['Coordinates']
                    xg = pos[:,0] - hpos[0]
                    yg = pos[:,1] - hpos[1]
                    zg = pos[:,2] - hpos[2]
                    dens = cutout['Density']
                    r2_g = xg**2 + yg**2 + zg**2
                    idx_g_in = np.where(r2_g < (rs*params['frac_rs'])**2)[0]

                if params['assign_at_peaks'] and len(idx_in) > 0 and len(idx_g_in) > 0:
                    x = x[idx_sort_r2_dm]
                    y = y[idx_sort_r2_dm]
                    z = z[idx_sort_r2_dm]
                    pos = np.array([x[idx_in], y[idx_in], z[idx_in]]).T
                    kdtree = sp.KDTree(pos)

                    if params['verbose']:
                        print('# of dm parts.: %d, # of cells: %d, rs: %.2f kpc'%
                            (len(idx_in),len(idx_g_in),rs))
                    tot_mass = np.sum(mass[idx_in])
                    tot_vol = (4*np.pi/3) * (rs*params['frac_rs'])**3
                    avg_dens = tot_mass / tot_vol
                    min_dens = 30 * avg_dens # in Msun / kpc^3
                    xg = xg[idx_g_in]
                    yg = yg[idx_g_in]
                    zg = zg[idx_g_in]
                    pos_g = np.array([xg,yg,zg]).T
                    kdtree_g = sp.KDTree(pos_g)
                    dens = dens[idx_g_in]
                    r2_g = r2_g[idx_g_in]
                    idx_sort_dens_g = np.argsort(dens)[::-1]
                    peaks = [0] # idx in idx_sort_dens_g
                    t0 = time.time()
                    tot = 0
                    suc = 0
                    idx_in_peaks = np.zeros(len(idx_sort_dens_g), dtype='i8')
                    idx = idx_sort_dens_g[0]
                    idxs2 = np.array(kdtree_g.query_ball_point(
                        np.array([xg[idx],yg[idx],zg[idx]]).T, 
                        params['rmin_peaks_finding']), dtype='i8')
                    idx_in_peaks[idxs2] = 1
                    g_dens = [dens[idx_sort_dens_g[0]]]
                    for k in range(1,len(idx_sort_dens_g)):
                        idx = idx_sort_dens_g[k]
                        if dens[idx] < min_dens: # TODO! unit conversion!
                            if params['verbose']:
                                print('No. %d of %d below min density'%(k+1, len(idx_sort_dens_g)))
                            break
                        in_peaks = False
                        if idx_in_peaks[idx] > 0:
                            in_peaks = True
                        if not in_peaks:
                            tot += 1
                            idxs = kdtree_g.query(np.array([xg[idx],yg[idx],zg[idx]]).T, 
                                k=np.min([16,len(idx_g_in)]))[1]
                            idxs = idxs[idxs!=idx]
                            if np.max(dens[idxs]) < dens[idx]:
                                suc += 1
                                peaks.append(k)
                                g_dens.append(dens[idx_sort_dens_g[k]])
                                idxs2 = np.array(kdtree_g.query_ball_point(
                                    np.array([xg[idx],yg[idx],zg[idx]]).T, 
                                    params['rmin_peaks_finding']), dtype='i8')
                                idx_in_peaks[idxs2] = 1
                    if params['verbose']:
                        print('# of candidate peaks: %d, # of valid peaks: %d'%(tot, suc))
                    if len(idx_sort_dens_g) == 0:
                        peaks = []
                    t1 = time.time()

                    x_peaks = xg[idx_sort_dens_g[peaks]]
                    y_peaks = yg[idx_sort_dens_g[peaks]]
                    z_peaks = zg[idx_sort_dens_g[peaks]]
                    peaks_dm = kdtree.query(np.array([x_peaks,y_peaks,z_peaks]).T, k=1)[1] # idx in idx_in

                    u, ind = np.unique(peaks_dm, return_index=True) # remove replication
                    peaks_dm = u[np.argsort(ind)].astype('i8')
                    g_dens = np.array(g_dens)[np.sort(ind)]
                    peaks_mass = []
                    peaks_idxs = []
                    peaks_surrounding = []
                    for k in range(len(peaks_dm)):
                        idxs_peak_i = np.array(kdtree.query_ball_point(
                            pos[peaks_dm[k]], params['peak_radius']), dtype='i8')
                        peaks_mass.append(np.sum(mass[idx_in[idxs_peak_i]]))
                        peaks_idxs.append(idxs_peak_i[idxs_peak_i!=peaks_dm[k]])
                        peaks_surrounding.extend(peaks_idxs[k])

                        pos_list.append(pos[peaks_dm[k]])
                        hid_list.append(hid_offset[i])

                    peaks_surrounding = np.unique(peaks_surrounding)
                    len_old = len(peaks_surrounding)

                    if params['evenly_distribute']:
                        peaks_surrounding = np.r_[list(zip_longest(*peaks_idxs))].flatten()
                        peaks_surrounding = peaks_surrounding[peaks_surrounding != np.array(None)]
                        u, ind = np.unique(peaks_surrounding, return_index=True)
                        peaks_surrounding = u[np.argsort(ind)].astype('i8')

                    assert len_old == len(peaks_surrounding)
                    if len(peaks_surrounding) == 0:
                        peaks_surrounding = []

                    t2 = time.time()
                    if params['verbose']:
                        print('time find peaks: %.2f s, time link parts.: %.2f s'%(t1-t0, t2-t1))
                        print('# of dm peaks: %d'%len(peaks_dm))
                        print('# of GCs: %d; # of parts. in peaks (old): %d'%(num_gc, len_old+len(peaks_dm)))
                        print('# of GCs: %d; # of parts. in peaks (new): %d'%(num_gc, len(peaks_surrounding)+len(peaks_dm)))

                    mask_not_peak = np.ones(len(idx_in), dtype=bool)
                    mask_not_peak[peaks_dm] = False
                    mask_not_peak[peaks_surrounding] = False
                    permuted_idx_in = np.r_[idx_in[peaks_dm], np.random.permutation(idx_in[peaks_surrounding]), 
                        np.random.permutation(idx_in[mask_not_peak])] 
                    dmid[idx_in] = dmid[permuted_idx_in]
                else:
                    dmid[idx_in] = np.random.permutation(dmid[idx_in])

                # avoid duplicated id
                k = 0
                num_now = 0
                while num_now < num_gc:
                    k += 1
                    if dmid[k] in gcid:
                        pass
                    else:
                        gcid.append(dmid[k])
                        if params['assign_at_peaks']:
                            if k < len(peaks_dm):
                                quality.append(3)
                            elif k < len(peaks_dm) + len(peaks_surrounding):
                                quality.append(2)
                            elif k < len(idx_in):
                                quality.append(1)
                            else:
                                quality.append(0)
                        else:
                            quality.append(0)
                        num_now += 1
                assert num_now == num_gc

            else: # params['collisionless_only'] == False
                # load cutout snapshot of the subhalo
                fields = ['Coordinates', 'GFM_StellarFormationTime', 'ParticleIDs']
                if 'full_snap_only' in params and params['full_snap_only']:
                    cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                        hid_offset_full_snap[i], snap_form_offset_full_snap[i], 'stars', fields)
                else:
                    cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                        hid_offset[i], snap_form_offset[i], 'stars', fields)

                if cutout['count'] == 0:
                    # if there is no sellar particles in this subhalo, use dm particles
                    fields = ['Coordinates', 'ParticleIDs']

                    if 'full_snap_only' in params and params['full_snap_only']:
                        cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                            hid_offset_full_snap[i], snap_form_offset_full_snap[i], 'dm', fields)
                    else:
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
                idx_valid_s = np.where( (a_s>=0) & (a_s<=1) & (r2_s<params['rmax_form']**2))[0] 
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
                    if 'full_snap_only' in params and params['full_snap_only']:
                        cutout = loader.load_halo(params['base_halo'], hid_root[j], 
                            hid_offset_full_snap[i], snap_form_offset_full_snap[i], 'dm', fields)
                    else:
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