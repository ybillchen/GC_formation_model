# Interface for simulations

Before running the model, you must convert the simulation outputs into a certain data structure. Currently, the model only works for TNG and simulations with ART. You may want to write your own interface to do the conversion.

## Data structure

For computational efficiency, we use the `hdf5` format. For each halo (galaxy), the model needs two files: merger tree and particle output. 

### Merger tree

- File name: `merger_tree_[haloID].hdf5`
	- `dataset (int)`: `SnapNum`: snapshot number
	- `dataset (int)`: `SubfindID`: ID of the halo in the halo finder
	- `dataset (int)`: `SubhaloID`: ID of the halo in the merger tree (may differ from `SubfindID`)
	- `dataset (float)`: `SubhaloMass`: total mass of the halo (in 10^10 Msun/h)
	- `dataset (int)`: `FirstProgenitorID`: `SubhaloID` of the first progenitor
	- `dataset (int)`: `MainLeafProgenitorID`: `SubhaloID` of the leaf halo along this branch
	- `dataset (int)`: `NextProgenitorID`: `SubhaloID` of the second progenitor
	- `dataset (int)`: `DescendantID`: `SubhaloID` of the first descendant
	- `dataset (float)`: `SubhaloPos`: positions of the halo (in comoving kpc/h)
	- `dataset (float)`: `SubhaloVel`: velocities of the halo (in km/s)

### Particle output

- File name: `halo_[haloID].hdf5`
	- `group`: `snap_[snapNum]_halo_[subhaloID]`: all subhalos in `merger_tree_[haloID].hdf5`
		- `group`: `dm` or `stars` or `gas`: particle type
			- `attrs (int)`: `count`: number of particles in this halo of a certain particle type
			- `dataset (float)`: `Coordinates`: coordinates of the particle (in comoving kpc/h)
			- `dataset (float)`: `Velocities`: velocities of the particle (in km/s)
			- `dataset (int)`: `ParticleIDs`: unique ID of the particle
			- `dataset (float)`: `Potential`: potential of the particle (in (km/s)^2/a)
			- `dataset (float)`: `GFM_StellarFormationTime`: formation time of star particle (in scale factor)