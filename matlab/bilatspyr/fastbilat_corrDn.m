% Apply EVM_Matlab's corrDn on both data and weights each depth layer 
%  Uses MEX implementation corrDnBatch
%
% Based upon bilateralFilter.m of Paris and Durand's 
%   "A Fast Approximation of the Bilateral Filter using a Signal Processing Approach", (ECCV 2006).
%   See http://people.csail.mit.edu/jiawen/#code
%
% created by Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%
function s = fastbilat_corrDn(x3d, lofilt, edgemethod, step, start)
    response = sum(lofilt(:));
    
    s = struct;
    s.gridData = corrDnBatch(x3d.gridData, lofilt, edgemethod, step, start);
    s.gridWeights = corrDnBatch(x3d.gridWeights, lofilt / response, edgemethod, step, start);
    size2d = [size(s.gridData,1), size(s.gridData,2)];

    s.di = x3d.di(1:size2d(1), 1:size2d(2));
    s.dj = x3d.dj(1:size2d(1), 1:size2d(2));
    s.map = x3d.map(start(1):step(1):end, start(2):step(2):end);
    s.size2d = size2d;
    %fprintf('%d x %d --> %d x %d\n', x3d.size2d(1), x3d.size2d(2), s.size2d(1), s.size2d(2));
end

