import os
import h5py
import illustris_python_true as il

def handle_merger_tree(sim_base, save_base, hid, fields=None, dtypes=None):

    if fields is None:
        fields = ['SubhaloMass', 'FirstProgenitorID', 'SubhaloID', 'SnapNum', 
            'MainLeafProgenitorID', 'NextProgenitorID', 'DescendantID', 
            'SubfindID', 'SubhaloPos']
    if dtypes is None:
        dtypes = ['f8', 'i8', 'i8', 'i8', 'i8', 'i8', 'i8', 'i8', 'f8']

    tree = il.sublink.loadTree(sim_base, 99, hid, fields, False)
    
    filename = save_base + 'merger_tree_%d.hdf5'%hid
    if os.path.exists(filename):
        os.remove(filename)

    with h5py.File(filename, 'w') as f:
        for field, dtype in zip(fields, dtypes):
            data = tree[field]
            f.create_dataset(field, data=data, dtype=dtype)

def handle_halo(sim_base, save_base, hid, parttypes=None, fields_list=None, dtypes_list=None):

    if parttypes is None:
        parttypes = ['dm', 'stars', 'gas']
    if fields_list is None:
        fields_list = [
            ['Coordinates', 'ParticleIDs'],
            ['Coordinates', 'GFM_StellarFormationTime', 'ParticleIDs'],
            ['Coordinates', 'ParticleIDs']]
    if dtypes_list is None:
        dtypes_list = [['f8', 'i8'], ['f8', 'f8', 'i8'], ['f8', 'i8']]

    tree = il.sublink.loadTree(sim_base, 99, hid, ['SubfindID', 'SnapNum'], False)

    filename = save_base + 'halo_%d.hdf5'%hid
    if os.path.exists(filename):
        os.remove(filename)

    with h5py.File(filename, 'w') as f:
        for h, s in zip(tree['SubfindID'], tree['SnapNum']):
            d = f.create_group('snap_%d_halo_%d'%(s,h))

            for parttype, fields, dtypes in zip(parttypes, fields_list, dtypes_list):
                d2 = d.create_group(parttype)
                cutout = il.snapshot.loadSubhalo(sim_base, s, h, parttype, fields=fields)

                for field, dtype in zip(fields, dtypes):
                    data = cutout[field]
                    d2.create_dataset(field, data=data, dtype=dtype)

if __name__ == '__main__':

    # base path of TNG50-1
    sim_base = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L35n2160TNG/output/'
    save_base = '/n/holyscratch01/vogelsberger/billchen/temp/'
    hid = 523889

    handle_merger_tree(sim_base, save_base, hid) 
    handle_halo(sim_base, save_base, hid)