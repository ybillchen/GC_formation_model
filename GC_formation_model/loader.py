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
        d = f['snap_%d_halo_%d'%(snap,hid)][parttype]
        res['count'] = d.attrs['count']
        for field in fields:
            res[field] = d[field][:]

    return res