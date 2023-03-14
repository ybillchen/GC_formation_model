# Licensed under BSD-3-Clause License - see LICENSE

import numpy as np

def offset(params):

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
                print(' NO. %d, halo id: %d'%(num_run,hid[num_run]))

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