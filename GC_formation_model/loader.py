# Licensed under BSD-3-Clause License - see LICENSE

import h5py

def load_merger_tree(base, hid, fields=None):

    res = {}

    if fields is None:
        fields = ['SubhaloMass', 'FirstProgenitorID', 'SubhaloID', 'SnapNum', 
            'MainLeafProgenitorID', 'NextProgenitorID', 'DescendantID', 
            'SubfindID', 'SubhaloPos']

    filename = base + 'merger_tree_%d.hdf5'%hid

    with h5py.File(filename, 'r') as f:
        for field in fields:
            res[field] = f[field][:]

    return res

def load_halo(base, hid_root, hid, snap, parttype, fields=None):

    res = {}

    if fields is None:
        fields = ['Coordinates', 'ParticleIDs']

    filename = base + 'halo_%d.hdf5'%hid_root

    with h5py.File(filename, 'r') as f:
        res['count'] = d.attrs['count']
        if res['count'] == 0:
            for field in fields:
                res[field] = []
        else:
            d = f['snap_%d_halo_%d'%(snap,hid)][parttype]
            for field in fields:
                res[field] = d[field][:]

    return res