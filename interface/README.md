# Interface for simulations

Before running the model, you must convert the simulation outputs into a certain data structure. Currently, the model only works for TNG and simulations with ART. You may want to write your own interface to do the conversion.

## Data structure

For computational efficiency, we use the `hdf5` format. For each halo (galaxy), the model needs two files: merger tree and particle output. 

### Merger tree

- File name: `merger_tree_[haloID].hdf5`
	- `dataset`: `SnapNum`: snapshot number
	- `dataset`: `SubfindID`: halo ID in the halo finder
	- `dataset`: `SubhaloID`: halo ID in the merger tree (may differ from `SubfindID`)
	- `dataset`: `SubhaloMass`: total mass of the halo in 1e10 Msun/h
	- `dataset`: `FirstProgenitorID`: `SubhaloID` of the first progenitor
	- `dataset`: `MainLeafProgenitorID`: `SubhaloID` of the leaf halo along this branch
	- `dataset`: `NextProgenitorID`: `SubhaloID` of the second progenitor
	- `dataset`: `DescendantID`: `SubhaloID` of the first descendant
	- `dataset`: `SubhaloPos`: positions of the halo in comoving kpc/h

### Particle output

- File name: `halo_[haloID].hdf5`
	- `group`: `snap_[snapNum]_halo_[subhaloID]`: all subhalos in `merger_tree_[haloID].hdf5`
		- `group`: `dm` or `stars` or `gas`: particle type
			- `dataset`: `Coordinates`: coordinates of the particle in comoving kpc/h
			- `dataset`: `ParticleIDs`: unique ID of the particle
			- `dataset`: `Potential`: potential of the particle in (km/s)^2/a
			- `dataset`: `GFM_StellarFormationTime`: formation time of `stars` particle in scale factor