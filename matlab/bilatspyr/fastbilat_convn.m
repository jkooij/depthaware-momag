% N-d convolution on both data and weights
%
% Based upon bilateralFilter.m of Paris and Durand's 
%   "A Fast Approximation of the Bilateral Filter using a Signal Processing Approach", (ECCV 2006).
%   See http://people.csail.mit.edu/jiawen/#code
%
% created by Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%
function x3d = fastbilat_convn(x3d, kernel3d)
    % convolution
    x3d.gridData = convn( x3d.gridData, kernel3d, 'same' );
    x3d.gridWeights = convn( x3d.gridWeights, kernel3d, 'same' );
end
