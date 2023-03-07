import h5py

def load_merger_tree(base, hid):

    res = {}

    fields = ['SubhaloMass', 'FirstProgenitorID', 'SubhaloID', 'SnapNum', 
        'MainLeafProgenitorID', 'NextProgenitorID', 'DescendantID', 
        'SubfindID']

    filename = base + 'merger_tree_%d.hdf5'%hid
    
    with h5py.File(filename, 'r') as f:
        for field in fields:
            res[field] = f[field][:]

    return res