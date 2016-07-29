# Depth-aware Motion Magnification

(c) Julian Kooij, 2016, Delft University of Technology


## Installation instructions

1. Clone this repository, e.g. to a local `depthaware-momag' directory

2. Download the [EVM Matlab code](http://people.csail.mit.edu/mrub/evm/#code) from the paper

	"Eulerian Video Magnification for Revealing Subtle Changes in the World"
	Hao-Yu Wu, Michael Rubinstein, Eugene Shih, John Guttag, Fredo Durand and William T. Freeman
	ACM Transactions on Graphics (Proc. SIGGRAPH 2012)

3. Unzip the EVM Matlab code to `depthaware-momag/matlab/external/`

3. Unzip the archive (see below) containing the example sequences `depthaware-momag/data/`, you should have directories
```
	depthaware-momag/data/sequence1/
	depthaware-momag/data/sequence2/
	depthaware-momag/data/sequence3/
	depthaware-momag/data/sequence4/
```

4. Compile the mex code. in Matlab go the the `depthaware-momag/matlab/` directory, and run
```
	% add paths
	startup

	% build mex code
	cd mex
	build_bilatspyr_mex_posix   % on Linux
	build_bilatspyr_mex_windows % for MS Windows
```


## Running depth-aware motion magnification
In Matlab go the the `depthaware-momag/matlab/` directory, and run
```
	startup % add paths
	run_all
```
Results should have been written to `depthaware-momag/output/`


## Data
Get the data from **TODO** (+/- 2GB).
Unzip each sequence archive in `depthaware-momag/data/`.

The frames of the sequence are stored as uncompressed MJPEG files (created in Matlab), which you should be able to open in Matlab or python/OpenCV.
Most important files are
*	`kinect.mj2`: the standard HD RGB color video
*	`kinect_depth.mj2`: depth frames in the original low-res depth image space
*	`kinect_map.mj2`: depth frames projected to the HD color space





