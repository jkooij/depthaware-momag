% Normalize each cell by dividing it by its weight
%
% Based upon bilateralFilter.m of Paris and Durand's 
%   "A Fast Approximation of the Bilateral Filter using a Signal Processing Approach", (ECCV 2006).
%   See http://people.csail.mit.edu/jiawen/#code
%
% created by Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%
function x3d = fastbilat_normalize(x3d)
    % normalize
    thresh = 0;
    mask = (x3d.gridWeights <= thresh);
    normalizedGrid = x3d.gridData ./ x3d.gridWeights;
    normalizedGrid( mask ) = 0; % put 0s where it's undefined
    
    x3d.gridData = normalizedGrid;
    x3d.gridWeights = cast(~mask, 'like', x3d.gridWeights);
end
