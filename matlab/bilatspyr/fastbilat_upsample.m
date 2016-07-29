% Sample extended representation at the depth coordinates.
%
% Based upon bilateralFilter.m of Paris and Durand's 
%   "A Fast Approximation of the Bilateral Filter using a Signal Processing Approach", (ECCV 2006).
%   See http://people.csail.mit.edu/jiawen/#code
%
% created by Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%
function output = fastbilat_upsample(x3d)
    if 0
        % interpn takes rows, then cols, etc
        % i.e. size(v,1), then size(v,2), ...
        output = interpn( x3d.gridData, x3d.di, x3d.dj, x3d.map, 'linear');
    else
        % faster
        mapl = floor(x3d.map);
        maph = ceil(x3d.map);
        alpha = x3d.map - mapl;

        d = size(x3d.gridData, 3);
        mapl(mapl < 1) = 1; maph(maph < 1) = 1;
        mapl(mapl > d) = d; maph(maph > d) = d;
        
        indl = sub2ind(size(x3d.gridData), x3d.di, x3d.dj, mapl);
        indh = sub2ind(size(x3d.gridData), x3d.di, x3d.dj, maph);
        
        output = x3d.gridData(indl) .* (1-alpha) + x3d.gridData(indh) .* alpha;
    end
end


