# Licensed under BSD-3-Clause License - see LICENSE

import numpy as np
from scipy.interpolate import interp1d

class MassLoss(object):
    def __init__(self, path):
        self.flost, self.t_solar, self.t_subsolar = np.loadtxt(path, usecols=(1,2,3), unpack=True)
        self.massFraction = np.vectorize(self.massFraction_scalar)

    # Metallicity dependent stellar evolution, as calculated by Prieto & Gnedin 
    def massFraction_scalar(self, fe_h, t_alive):
        ssub = interp1d(self.t_subsolar, self.flost, kind='linear')
        ssolar = interp1d(self.t_solar, self.flost, kind='linear')

        t_alive *= 1e9 # transfer Gyr to yr
        if fe_h >= 0: # supersolar metallicity
            return ssolar(t_alive)
        elif fe_h <= -0.69: # if Z<= 0.2Zsun
            return ssub(t_alive)
        else: # intermediate case: interpolate between .2*Z_sun and Z_sun
            slope = (self.t_solar - self.t_subsolar)/0.69

            # new intermediate MS lifetime table
            t_int= self.t_subsolar + slope*(fe_h - -0.69) 

            sint = interp1d(t_int, self.flost, kind = 'linear') 
            return sint(t_alive)

# fix missing snaps for P
def fix_P(P, tag):
    for i in range(len(P)):
        if tag[i][-1] == 0:
            continue
        for j in range(len(P[0])-1, -1, -1):
            if tag[i][j] == 0:
                P[i][j] = P[i][j+1]
                tag[i][j] = tag[i][j+1]
    return P, tag

def evolve(params, snap_range=None, return_t_disrupt=False):

    if params['verbose']:
        print('\n########## evolution started ##########')

    z_list = params['redshift_snap']
    redshift_snap = params['redshift_snap']
    full_snap = params['full_snap']
    fin_name = params['resultspath'] + params['allcat_name']
    gcid_name = fin_name[:-4] + '_gcid.txt'
    root_name = fin_name[:-4] + '_offset_root.txt'
    offset_name = fin_name[:-4] + '_offset.txt'
    tid1_name = fin_name[:-4] + '_tideig1.txt'
    tid2_name = fin_name[:-4] + '_tideig2.txt'
    tid3_name = fin_name[:-4] + '_tideig3.txt'
    tag_name = fin_name[:-4] + '_tag.txt'

    if snap_range is None:
        snap_range = len(full_snap)

    ml = MassLoss(params['path_massloss'])

    # load GC catalog
    hid, logmh, logms, logmh_form, logms_form, logm_form, z_form, feh, \
        ismpb, hid_form, snapnum_form = np.loadtxt(fin_name, unpack=True)

    # setting indices to int type
    hid = hid.astype(int)
    ismpb = ismpb.astype(int)
    hid_form = hid_form.astype(int)
    snapnum_form = snapnum_form.astype(int)

    z_form = z_list[snapnum_form] 
    t_form = params['cosmo'].cosmicTime(z_form, units='Gyr')

    # load GC id
    gcid, quality = np.loadtxt(gcid_name, ndmin=2, unpack=True, dtype='int64')

    # load root offset
    hid_root, idx_beg, idx_end, idx_beg_in_off, idx_end_in_off = np.loadtxt(
        root_name, ndmin=2, unpack=True, dtype='int64')

    # calculate P_inverse and load tag
    if params['disrupt_mode'] == 'constant':
        P_inverse = np.ones([len(gcid), len(full_snap)]) / params['pr']
        tag = np.ones(P_inverse.shape, dtype=int)

    elif params['disrupt_mode'] == 'tidal':
        tideig1 = np.loadtxt(tid1_name)
        tideig2= np.loadtxt(tid2_name)
        tideig3 = np.loadtxt(tid3_name)
        tideige = tideig1 - tideig3
        P_inverse = params['kappa'] * np.sqrt(tideige) / 150
        tag_name = fin_name[:-4] + '_tidtag.txt'
        tag = np.loadtxt(tag_name, dtype='int64')
        P_inverse, tag = fix_P(P_inverse, tag)

    else:
        if params['verbose']:
            print('Error disrupt_mode!')
        return 0

    m_now = 10**logm_form
    snap_now = np.copy(snapnum_form)
    t_disrupt = -1 * np.ones(len(m_now))

    for i in range(len(hid_root)):
        if params['verbose']:
            print(' NO. %d, subhalo id: %d'%(i,hid_root[i]))
        # t0 = time.time()

        for j in range(snap_range):
            # t1 = time.time()

            snap = full_snap[j]
            z_snap = z_list[snap]
            scale_a = 1 / (1 + z_snap)

            # existing and survived GCs at this snapshot
            idx_exist_gc = np.where((snapnum_form[idx_beg[i]:idx_end[i]] < snap) & 
                (tag[idx_beg[i]:idx_end[i],snap_range-1] > 0) & 
                (m_now[idx_beg[i]:idx_end[i]] > 0))[0]

            # jump to the next if no existing GC
            if len(idx_exist_gc) == 0:
                continue
            idx_exist_gc += idx_beg[i] # idx in the whole catalog

            # redshift and time of GCs before disruption
            z_now = z_list[snap_now[idx_exist_gc]] 
            t_now = params['cosmo'].cosmicTime(z_now, units = 'Gyr')
            # cosmic time of the current snapshot
            t_snap = params['cosmo'].cosmicTime(z_snap, units = 'Gyr')

            mi = 10**logm_form[idx_exist_gc]
            t_tid = 10 * (mi/2e5)**params['disrupt_x'] * (m_now[idx_exist_gc]/mi)**params['disrupt_y'] / \
                np.clip(P_inverse[idx_exist_gc,j], 1e-10, None)

            t_iso = 1e10 # i.e., no iso disruption. 17 * (m_now[idx_exist_gc]/2e5) # in Gyr

            idx_tid = np.where(np.greater(t_iso, t_tid))[0]
            idx_iso = np.where(np.less_equal(t_iso, t_tid))[0]

            if len(idx_tid):
            # tidal dominated
                gamma = params['disrupt_y']
                dt = t_tid[idx_tid] / gamma
                k1 = 1 - (t_snap-t_now[idx_tid]) / dt
                k1 = np.clip(k1, 0, 1)
                m_now[idx_exist_gc[idx_tid]] = m_now[idx_exist_gc[idx_tid]] * k1**(1/gamma)

                idx_d = np.where((dt - (t_snap-t_now[idx_tid])) <= 0)[0] # disrupted
                t_disrupt[idx_exist_gc[idx_tid[idx_d]]] = t_now[idx_tid[idx_d]] + dt[idx_d]

            if len(idx_iso) > 0:
            # iso dominated
                m_now[idx_exist_gc[idx_iso]] = m_now[idx_exist_gc[idx_iso]] - \
                    (t_snap-t_now[idx_iso]) * (2e5/17) * np.ones(len(idx_iso))

            # update snap_now
            snap_now[idx_exist_gc] = snap * np.ones(len(idx_exist_gc), dtype=int)

        if len(idx_exist_gc):
            m_now[idx_exist_gc] = m_now[idx_exist_gc] * (1-ml.massFraction(feh[idx_exist_gc], 
                t_snap - t_form[idx_exist_gc]))

    idx_valid = np.where(
        (snapnum_form < full_snap[snap_range-1]) & (tag[:,snap_range-1]>0) & (m_now > 0))[0]

    m_out = -1 * np.ones(len(m_now))
    m_out[idx_valid] = np.log10(m_now[idx_valid])

    # save data
    if params['disrupt_mode'] == 'constant':
        np.savetxt(fin_name[:-4]+'_p-%g_logm_snap%d.txt'%(params['pr'],full_snap[snap_range-1]), m_out, fmt='%.3f')
    elif params['disrupt_mode'] == 'tidal':
        np.savetxt(fin_name[:-4]+'_k-%g_logm_snap%d.txt'%(params['kappa'],full_snap[snap_range-1]), m_out, fmt='%.3f')

    if return_t_disrupt:
        if params['disrupt_mode'] == 'constant':
            np.savetxt(fin_name[:-4]+'_p-%g_t_disrupt_snap%d.txt'%(params['pr'],full_snap[snap_range-1]), t_disrupt, fmt='%.3f')
        elif params['disrupt_mode'] == 'tidal':
            np.savetxt(fin_name[:-4]+'_k-%g_t_disrupt_snap%d.txt'%(params['kappa'],full_snap[snap_range-1]), t_disrupt, fmt='%.3f')

    if params['verbose']:
        print('########## evolution done ##########')