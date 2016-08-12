# Depth-aware Motion Magnification

(c) Julian Kooij, 2016, Delft University of Technology


## Installation instructions

1. Clone this repository, e.g. to a local `depthaware-momag' directory

2. Install Eero P. Simoncelli's matlabPyrTools, e.g. in the `external/` subdirectory:
```
	cd depthaware-momag/external/
	git clone https://github.com/LabForComputationalVision/matlabPyrTools
```
	Follow the instructions to build the mex files in the matlabPyrTools/MEX/ directory.

	*Note*: you might get an error in `reconSpyr.m` line 95, and `reconSpyrLevs.m` line 41.
	In that case, replace the start of those lines,
```
	res = upConv(...
```
by
```
	upConv( ...
```
	
	You can then either add `matlabPyrTools/` and `matlabPyrTools/MEX/` to your Matlab path,
	or place the dependency in `depthaware-momag/external/` such that `depthaware-momag/startup.m` will add it to your path when needed.

3. Unzip the data archive(s) (see below) containing the example sequences `depthaware-momag/data/`. You should now have directories
```
	depthaware-momag/data/sequence1/
	depthaware-momag/data/sequence2/
	depthaware-momag/data/sequence3/
	depthaware-momag/data/sequence4/
```

4. Compile the mex code. in Matlab go the the `depthaware-momag/matlab/` directory, and run
```
	% go to the project directory
	cd depthaware-momag/matlab

	% add paths
	startup

	% build mex code
	cd mex
	build_bilatspyr_mex_posix   % for Linux
	build_bilatspyr_mex_windows % for MS Windows
```
*NOTE*: to compile the mex code, you might need to adjust the include/library paths, and/or download the [Eigen template library](http://eigen.tuxfamily.org/index.php?title=Main_Page) first.

## Running depth-aware motion magnification
In Matlab go the the `depthaware-momag/matlab/` directory, and run
```
	startup % add paths
	run_all
```
Results should have been written to `depthaware-momag/output/`

If you open `run_all`, you can see that 


## Data
Get the data from the [Technology in Motion (TIM) website](https://tim.lumc.nl/site/en/research/downloads/).

All sequences together are +/- 2GB.
Unzip each sequence archive in `depthaware-momag/data/`.

The frames of the sequence are stored as uncompressed MJPEG files (created in Matlab), which you should be able to open in Matlab or python/OpenCV.
Most important files are
*	`kinect.mj2`: the standard HD RGB color video
*	`kinect_depth.mj2`: depth frames in the original low-res depth image space
*	`kinect_map.mj2`: depth frames projected to the HD color space





