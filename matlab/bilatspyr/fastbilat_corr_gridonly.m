% Apply EVM_Matlab's corrDn (on data ONLY) at each depth layer, but WITHOUT 'Dn' subsampling
%  Uses MEX implementation corrDnBatch
%
% Based upon bilateralFilter.m of Paris and Durand's 
%   "A Fast Approximation of the Bilateral Filter using a Signal Processing Approach", (ECCV 2006).
%   See http://people.csail.mit.edu/jiawen/#code
%
% created by Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%
% SEE ALSO fastbilat_corrDn_gridonly
function s = fastbilat_corr_gridonly(x3d, lofilt, edgemethod)
    s = x3d;
    if isreal(lofilt)
        s.gridData = corrDnBatch(x3d.gridData, lofilt, edgemethod);
    else
        lofilt_r = real(lofilt);
        lofilt_i = imag(lofilt);
        gdata_r = corrDnBatch(x3d.gridData, lofilt_r, edgemethod);
        gdata_i = corrDnBatch(x3d.gridData, lofilt_i, edgemethod);
        s.gridData = gdata_r + 1j .* gdata_i;
    end
end
