import h5py
import illustris_python_true as il

def handle_merger_tree(sim_base, save_base, hid):

    fields = ['SubhaloMass', 'FirstProgenitorID', 'SubhaloID', 'SnapNum', 
        'MainLeafProgenitorID', 'NextProgenitorID', 'DescendantID', 'SubfindID']
    dtypes = ['f8', 'i8', 'i8', 'i8', 'i8', 'i8', 'i8', 'i8']

    tree = il.sublink.loadTree(sim_base, 99, hid, fields, False)

    filename = save_base + 'merger_tree_%d.hdf5'%hid

    with h5py.File(filename, 'w') as f:
        for field, dtype in zip(fields, dtypes):
            data = tree[field]
            f.create_dataset(field, data=data, dtype=dtype)

if __name__ == '__main__':

    # base path of TNG50-1
    sim_base = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L35n2160TNG/output/'
    save_base = '/n/holyscratch01/vogelsberger/billchen/temp/'
    hid = 523889

    handle_merger_tree(sim_base, save_base, hid) 