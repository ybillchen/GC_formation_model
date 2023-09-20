# Licensed under BSD-3-Clause License - see LICENSE

import numpy as np

from . import loader

def get_offset(params):

    if params['verbose']:
        print('\n########## offset calculation started ##########')

    fin_name = params['resultspath'] + params['allcat_name']

    hid, logmh, logms, logmh_form, logms_form, logm_form, z_form, feh, \
        ismpb, hid_form, snapnum_form = np.loadtxt(fin_name, ndmin=2, unpack=True)

    # setting indices to int type
    hid = hid.astype(int)
    ismpb = ismpb.astype(int)
    hid_form = hid_form.astype(int)
    snapnum_form = snapnum_form.astype(int)

    i = 0 # pointer of the begining of a group of GCs formed at a same subhalo
    num_run = 0

    # root offset, for each galaxy
    subfindids_root = []
    begins_root = []
    ends_root = []

    # detailed offset, for each GC formation event
    snapnums = []
    subfindids = []
    begins = []
    ends = []

    # find the group of GCs formed at a same subhalo. note, those GCs must
    # be stored together, so the following O(n) algorithm is good enough
    for j in range(len(hid) + 1):
        # at the end of a subhalo
        if (j >= len(hid)) or (hid[np.clip(j-1,0,len(hid)-1)] != hid[np.clip(j,0,len(hid)-1)]):
            ends_root.append(j)

        # at the begining of a subhalo
        if j == 0 or hid[j-1] != hid[np.clip(j,0,len(hid)-1)]:
            if params['verbose']:
                print(' NO. %d, halo id: %d'%(num_run,params['subs'][num_run]))

            num_run += 1

            subfindids_root.append(hid[j])
            begins_root.append(j)

        if (j < len(hid)) and (snapnum_form[i] == snapnum_form[np.clip(j,0,len(hid)-1)]) and (hid_form[i] == hid_form[np.clip(j,0,len(hid)-1)]):
            continue

        snapnums.append(snapnum_form[i])
        subfindids.append(hid_form[i])
        begins.append(i)
        ends.append(j)

        i = j # update i

    # calculate locations of root offset in detailed offset
    begins_in_detail = [0]
    ends_in_detail = []

    idx_in_root = 0
    for k in range(len(ends)):
        if ends[k] > ends_root[idx_in_root]:
            begins_in_detail.append(k)
            ends_in_detail.append(k)
            idx_in_root += 1
    ends_in_detail.append(len(ends))

    output_root = np.array([subfindids_root, begins_root, ends_root, 
        begins_in_detail, ends_in_detail], dtype=int).T
    header = 'SubfindID(z=0) | BeginIdx | EndIdx+1 | BeginIdx in offset | EndIdx in offset'
    np.savetxt(fin_name[:-4]+'_offset_root.txt', output_root, fmt='%d ', header=header)

    output = np.array([snapnums, subfindids, begins, ends], dtype=int).T
    header = 'Snapnum | SubfindID | BeginIdx | EndIdx+1'
    np.savetxt(fin_name[:-4]+'_offset.txt', output, fmt='%d ', header=header)

    if params['verbose']:
        print('########## offset calculation done ##########')

def find_next_full_snap(params):
    if params['verbose']:
        print('\n########## finding next full snap started ##########')

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

    snap_form_offset_new = np.copy(snap_form_offset)
    hid_offset_new = np.copy(hid_offset)

    for j in range(len(hid_root)):
        if params['verbose']:
            print(' NO. %d, halo id: %d'%(j,hid_root[j]))


        tree = loader.load_merger_tree(params['base_tree'], hid_root[j])
        mpbi = tree['MainLeafProgenitorID'][(tree['SubfindID']==hid_root[j])&\
            (tree['SnapNum']==np.max(params['full_snap']))][0]

        # loop over all GC formation events
        for i in range(idx_beg_in_off[j], idx_end_in_off[j]):
            if snap_form_offset[i] in params['full_snap']:
                continue

            snap = snap_form_offset[i] 

            idx_in_tree = np.where( (tree['SnapNum']==snap) &
                (tree['SubfindID']==hid_offset[i]) )[0][0]

            desc_id = tree['DescendantID'][idx_in_tree]

            idx_desc = np.where( (tree['SnapNum']>snap) &
                (tree['SubhaloID']==desc_id) )[0]
            assert len(idx_desc) == 1
           
            snap = tree['SnapNum'][idx_desc[0]]

            while not snap in params['full_snap']:

                desc_id = tree['DescendantID'][idx_desc[0]]

                idx_desc = np.where( (tree['SnapNum']>snap) &
                    (tree['SubhaloID']==desc_id) )[0]
                assert len(idx_desc) == 1

                snap = tree['SnapNum'][idx_desc[0]]

            snap_form_offset_new[i] = snap
            hid_offset_new[i] = tree['SubfindID'][idx_desc[0]]

    output = np.array([snap_form_offset_new, hid_offset_new, idx_beg, idx_end], dtype=int).T
    header = 'Snapnum | SubfindID | BeginIdx | EndIdx+1'
    np.savetxt(fin_name[:-4]+'_offset_full_snap.txt', output, fmt='%d ', header=header)

    if params['verbose']:
        print('########## finding next full snap done ##########')


def offset(params):
    get_offset(params)
    if 'full_snap_only' in params and params['full_snap_only']:
        find_next_full_snap(params)