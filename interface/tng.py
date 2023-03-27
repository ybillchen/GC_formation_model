import os
import numpy as np
import h5py
from multiprocessing import Pool
import illustris_python_true as il

h100 = 0.6774
full_snap = [2,3,4,6,8,11,13,17,21,25,33,40,50,59,67,72,78,84,91,99]

def handle_merger_tree(sim_base, save_base, hid, fields=None, dtypes=None):

    if fields is None:
        fields = ['SubhaloMass', 'FirstProgenitorID', 'SubhaloID', 'SnapNum', 
            'MainLeafProgenitorID', 'NextProgenitorID', 'DescendantID', 
            'SubfindID', 'SubhaloPos', 'SubhaloVel']
    if dtypes is None:
        dtypes = ['f8', 'i8', 'i8', 'i8', 'i8', 'i8', 'i8', 'i8', 'f8', 'f8']

    tree = il.sublink.loadTree(sim_base, 99, hid, fields, False)

    filename = save_base + 'merger_tree_%d.hdf5'%hid
    if os.path.exists(filename):
        os.remove(filename)

    with h5py.File(filename, 'w') as f:
        for field, dtype in zip(fields, dtypes):
            f.create_dataset(field, data=tree[field], dtype=dtype)

    print('handle_merger_tree hid: %d done'%hid)

def handle_halo(sim_base, save_base, hid, mh_min=1e8, parttypes=None, fields_list=None, dtypes_list=None):

    if parttypes is None:
        parttypes = ['dm', 'stars', 'gas']
    if fields_list is None:
        fields_list = [
            ['Coordinates', 'Velocities', 'ParticleIDs', 'Potential'],
            ['Coordinates', 'Velocities', 'ParticleIDs', 'GFM_StellarFormationTime', 'Potential'],
            ['Coordinates', 'Velocities', 'ParticleIDs', 'Potential']]
    if dtypes_list is None:
        dtypes_list = [['f8', 'f8', 'i8', 'f8'], ['f8', 'f8', 'i8', 'f8', 'f8'], ['f8', 'f8', 'i8', 'f8']]

    tree = il.sublink.loadTree(sim_base, 99, hid, ['SubhaloMass', 'SubfindID', 'SnapNum'], False)
    ntot = len(tree['SubhaloMass'])

    filename = save_base + 'halo_%d.hdf5'%hid
    if os.path.exists(filename):
        os.remove(filename)


    with h5py.File(filename, 'w') as f:
        for h, s, m in zip(tree['SubfindID'], tree['SnapNum'], tree['SubhaloMass']):
            
            d = f.create_group('snap_%d_halo_%d'%(s,h))

            for parttype, fields, dtypes in zip(parttypes, fields_list, dtypes_list):
                if not s in full_snap:
                    fs = fields[:-1]
                    ds = dtypes[:-1]
                else:
                    fs = fields
                    ds = dtypes
                d2 = d.create_group(parttype)
                cutout = il.snapshot.loadSubhalo(sim_base, s, h, parttype, fields=fs)

                if cutout['count'] == 0: # no particle
                    d2.attrs['count'] = 0
                    for field, dtype in zip(fs, ds):
                        d2.create_dataset(field, data=[], dtype=dtype)

                else:
                    d2.attrs['count'] = len(cutout[fs[0]])
                    for field, dtype in zip(fs, ds):
                        d2.create_dataset(field, data=cutout[field], dtype=dtype)

    print('handle_halo hid: %d done'%hid)

if __name__ == '__main__':

    # base path of TNG50-1
    sim_base = '$sim_base'
    save_base = '$save_base'
    hid = 523889

    handle_merger_tree(sim_base, save_base, hid) 
    handle_halo(sim_base, save_base, hid)