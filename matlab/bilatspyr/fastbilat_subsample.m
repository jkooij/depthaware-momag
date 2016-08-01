% Subsample the extended representation, keep cells at given step interval
%
% Based upon bilateralFilter.m of Paris and Durand's 
%   "A Fast Approximation of the Bilateral Filter using a Signal Processing Approach", (ECCV 2006).
%   See http://people.csail.mit.edu/jiawen/#code
%
% created by Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%
function x3d = fastbilat_subsample(x3d, step, start)
    x3d.gridData = x3d.gridData(start(1):step(1):end, start(2):step(2):end, :);
    x3d.gridWeights = x3d.gridWeights(start(1):step(1):end, start(2):step(2):end, :);
    size2d = [size(x3d.gridData,1), size(x3d.gridData,2)];
    x3d.di = x3d.di(1:size2d(1), 1:size2d(2));
    x3d.dj = x3d.dj(1:size2d(1), 1:size2d(2));
    x3d.map = x3d.map(start(1):step(1):end, start(2):step(2):end);
    x3d.size2d = size2d;
end
