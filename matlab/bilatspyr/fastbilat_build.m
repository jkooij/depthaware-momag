% Create an extended representation for fast bilateral filtering
%
% Based upon bilateralFilter.m of Paris and Durand's 
%   "A Fast Approximation of the Bilateral Filter using a Signal Processing Approach", (ECCV 2006).
%   See http://people.csail.mit.edu/jiawen/#code
%
% created by Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%
function s = fastbilat_build(im, dMapIndices, downsampledDepth)

    [inputHeight, inputWidth] = size(im);
    
    gridData = zeros( inputHeight, inputWidth, downsampledDepth );
    gridWeights = zeros( inputHeight, inputWidth, downsampledDepth );
    
    % compute downsampled indices
    [ jj, ii ] = meshgrid( 0 : inputWidth - 1, 0 : inputHeight - 1 );

    di = ii + 1;
    dj = jj + 1;
    dz = round( dMapIndices );
    dz(dz < 1) = 1; dz(dz > downsampledDepth) = downsampledDepth;
    
    % perform scatter (there's probably a faster way than this)
    % normally would do downsampledWeights( di, dj, dk ) = 1, but we have to
    % perform a summation to do box downsampling
    for k = 1 : numel( dz ),

        dataZ = im( k ); % traverses the image column wise, same as di( k )
        if isnan( dataZ ) continue; end

        dik = di( k );
        djk = dj( k );
        dzk = dz( k );
        if isnan(dzk); continue; end
        if (dzk < 1 || dzk > downsampledDepth); continue; end

        gridData( dik, djk, dzk ) = gridData( dik, djk, dzk ) + dataZ;
        gridWeights( dik, djk, dzk ) = gridWeights( dik, djk, dzk ) + 1;
    end
    
    if 0
        % DEBUG?
        % Improve visual quality of motion magnification by
        % keeping background values to ALL occlusion layers (which are
        % closer to the camera)
        gridData = cumsum(gridData, 3);
        gridWeights = cumsum(gridWeights, 3);
    end
    
    s = struct;
    s.gridData = gridData;
    s.gridWeights = gridWeights;
    s.di = di;
    s.dj = dj;
    s.map = dMapIndices;
    s.size2d = [size(s.gridData,1), size(s.gridData,2)];
end