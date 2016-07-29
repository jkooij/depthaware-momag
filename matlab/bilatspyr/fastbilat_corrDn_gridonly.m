% Apply EVM_Matlab's corrDn (on data ONLY) at each depth layer, including 'Dn' subsampling
%  Uses MEX implementation corrDnBatch
%
% Based upon bilateralFilter.m of Paris and Durand's 
%   "A Fast Approximation of the Bilateral Filter using a Signal Processing Approach", (ECCV 2006).
%   See http://people.csail.mit.edu/jiawen/#code
%
% created by Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%
% SEE ALSO fastbilat_corr_gridonly
function x3d = fastbilat_corrDn_gridonly(x3d, lofilt, edgemethod, step, start)
    x3d.gridData = corrDnBatch(x3d.gridData, lofilt, edgemethod, step, start);
    x3d.gridWeights = x3d.gridWeights(start(1):step(1):end, start(2):step(2):end, :);
    
    size2d = [size(x3d.gridData,1), size(x3d.gridData,2)];

    x3d.di = x3d.di(1:size2d(1), 1:size2d(2));
    x3d.dj = x3d.dj(1:size2d(1), 1:size2d(2));
    x3d.map = x3d.map(start(1):step(1):end, start(2):step(2):end);
    x3d.size2d = size2d;    
end

