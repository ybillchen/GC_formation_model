# Interface for simulations

Before running the model, you must convert the simulation outputs into a certain data structure. Currently, the model only works for TNG and simulations with ART. You may want to write your own interface to do the conversion.

## Data structure

For computational efficiency, we use the `hdf5` format.

### Merger tree

- File name: `merger_tree_[haloID].hdf5`

### Particle outputs

- File name: `halo_[haloID].hdf5`