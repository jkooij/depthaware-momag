% [PYR, INDICES, STEERMTX, HARMONICS] = BUILD_BILATSPYR(IM, HEIGHT, FILTFILE, EDGES)
%
% NOTE: Adapted from buildSpyr, originally by
%   Eero Simoncelli, 6/96.
%   See http://www.cis.upenn.edu/~eero/steerpyr.html for more
%   information about the Steerable Pyramid image decomposition.
%
% Construct a steerable pyramid on matrix IM.  Convolutions are
% done with spatial filters.
%
% HEIGHT (optional) specifies the number of pyramid levels to build. Default
% is maxPyrHt(size(IM),size(FILT)). 
% You can also specify 'auto' to use this value.
%
% FILTFILE (optional) should be a string referring to an m-file that
% returns the rfilters.  (examples: 'sp0Filters', 'sp1Filters',
% 'sp3Filters','sp5Filters'.  default = 'sp1Filters'). EDGES specifies
% edge-handling, and defaults to 'reflect1' (see corrDn).
%
% PYR is a vector containing the N pyramid subbands, ordered from fine
% to coarse.  INDICES is an Nx2 matrix containing the sizes of
% each subband.  This is compatible with the MatLab Wavelet toolbox.
% See the function STEER for a description of STEERMTX and HARMONICS.
%

function [pyr,pind,precomp] = build_bilatspyr(im, ht, filtfile, edgemethod, dMap, dMin, dMax, dSigma)
    %-----------------------------------------------------------------
    % DEFAULTS:

    if (exist('filtfile') ~= 1) || isempty(filtfile)
      filtfile = 'sp1Filters';
    end

    if (exist('edges') ~= 1) || isempty(edgemethod)
      edgemethod= 'reflect1';
    end

    if (isstr(filtfile) & (exist(filtfile) == 2))
       [lo0filt,hi0filt,lofilt,bfilts,steermtx,harmonics] = eval(filtfile);
       precomp = [];
    elseif isstruct(filtfile)
       precomp = filtfile;
    else
      fprintf(1,'\nUse buildSFpyr for pyramids with arbitrary numbers of orientation bands.\n');
      error('FILTFILE argument must be the name of an M-file containing SPYR filters.');
    end
       
    if isempty(precomp)
        precomp = struct;
        precomp.steermtx = steermtx;
        precomp.harmonics = harmonics;

        precomp.lofilt = lofilt;
        precomp.lo0filt = lo0filt;
        precomp.hi0filt = hi0filt;
        precomp.bfilts = bfilts;
        precomp.hi0filt = hi0filt;
    else
        lofilt = precomp.lofilt;
        lo0filt = precomp.lo0filt;
        hi0filt = precomp.hi0filt;
        bfilts = precomp.bfilts;        
    end
    
    max_ht = maxPyrHt(size(im), size(lofilt,1));
    if ( (exist('ht') ~= 1) | (ht == 'auto') )
      ht = max_ht;
    else
      if (ht > max_ht)
        error(sprintf('Cannot build pyramid higher than %d levels.',max_ht));
      end
    end
    
    dDelta = dMax - dMin;
    dRange = dSigma;

    derivedDRange = dSigma / dRange;
    kernelDepth = 4 * derivedDRange + 1;
    halfKernelDepth = floor( kernelDepth / 2 );
    
    %paddingZ = floor( kernelDepth-1 ) + 1;
    paddingZ = 0; % FIXME: is paddingZ needed?
    
    downsampledDepth = floor( dDelta / dRange ) + 1 + 2 * paddingZ;

    % DEBUG: set kernelDepth to size of downsampled depth
    %kernelDepth = 2 * floor(downsampledDepth/2) + 1;
    %halfKernelDepth = floor( kernelDepth / 2 );

    % create depth kernel
    gridZ = reshape([0 : kernelDepth - 1] - halfKernelDepth, 1, 1, []);
    kernelZ = exp( -.5 * (gridZ .* gridZ) / derivedDRange ); 

    dMapNorm = ( dMap - dMin ) / dRange + paddingZ + 1;    
    
    % compute for each layer the depth at the center of the bin
    dcenters = ([1:downsampledDepth] - 1 - paddingZ) * dRange + dMin;
    
    %% --- create bilateral pyramid
    % -- step 0 ---
    % build 3D extended representation of image
    im3d = fastbilat_build(im, dMapNorm, downsampledDepth);
    
    % keep track of valid input data
    origD = im3d.gridData;
    origW = im3d.gridWeights;
    mask = (im3d.gridWeights ~= 0);

    % -- step 1 ---
    % apply 2D spatial smoothing kernel
    im3d = fastbilat_tdistconv(im3d, .1);
    %im3d = fastbilat_gaussconv(im3d, 1);
    
    % perform 1D depth smoothing kernel
    im3d = fastbilat_convn(im3d, kernelZ);
    
    % -- step 2 ---
    % element-wise division by weight total
    im3d = fastbilat_normalize(im3d);

    % -- step 3 ---
    % replace data again
    im3d.gridData(mask) = origD(mask);
    im3d.gridWeights(mask) = origW(mask);
    %% ---
    
    % intial lowpass/highpass filters
    lo03d = fastbilat_corrDn(im3d, lo0filt, edgemethod, [1 1], [1 1]);
    hi03d = fastbilat_corrDn(im3d, hi0filt, edgemethod, [1 1], [1 1]);
    
    %% build level recursively
    outputs = buildSpyrLevs(lo03d, ht, lofilt, bfilts, edgemethod, kernelZ);
    outputs = {hi03d, outputs{:}};
    
    %[pyr, pind] = out_to_pyr_pind(outputs);
    [pyr, pind] = out_to_pyr_pind_map(outputs);
    pyr.dcenters = dcenters;
end

function [pyr, pind] = out_to_pyr_pind_map(outputs)
    % number of sample points per output band
    nlevels = numel(outputs);
    pind = [];
    for lvl = 1:nlevels
        output = outputs{lvl};   
        for band = 1:numel(output)
            bsize = size(output(band).gridData);
            pind = [pind; bsize];
        end
    end
    if size(pind,2) < 3; pind(:,3) = 1; end;
    
    blen = pind(:,1).*pind(:,2); % length (# of elements) of each band
    brange = [cumsum([1; blen(1:end-1)]), cumsum(blen)]; % (start,stop) index per band
    plen = brange(end); % total length of pyramid
    
    depth = pind(1,3);
    
    pyr = zeros(plen, depth, 'single');
    pmap = {};
    
    j = 0;
    for lvl = 1:numel(outputs)
        output = outputs{lvl};
       
        for b = 1:numel(output)
            j = j + 1;
            
            band = output(b);
            
            bpyr = band.gridData;
            bpyr = reshape(bpyr, [], depth);
            bmap = band.map;

            pyr(brange(j,1):brange(j,2),:) = bpyr;
        end
        
        bmap = single(bmap);
        pmap{end+1} = bmap;
    end
    
    out = struct;
    out.pyr = pyr;
    out.pmap = pmap;
    
    pyr = out;
end

function x3d = fastbilat_tdistconv(x3d, sSigma)
    % convolution with t-distribution kernel
    x3d.gridData = tdistconv( x3d.gridData, sSigma );
    x3d.gridWeights = tdistconv( x3d.gridWeights, sSigma );
end

function x3d = fastbilat_gaussconv(x3d, sSigma)
    % convolution with gauss kernel
    x3d.gridData = gaussconv( x3d.gridData, sSigma );
    x3d.gridWeights = gaussconv( x3d.gridWeights, sSigma );
end

function outputs = buildSpyrLevs(lo03d,ht,lofilt,bfilts,edgemethod,kernelZ)
    % [PYR, INDICES] = buildSpyrLevs(LOIM, HEIGHT, LOFILT, BFILTS, EDGES)
    %
    % Recursive function for constructing levels of a steerable pyramid.  This
    % is called by buildSpyr, and is not usually called directly.

    % Eero Simoncelli, 6/96.

    cur3d = lo03d;
    %cur3d = fastbilat_corrDn(cur3d, lofilt_G, edgemethod, [1 1], [1 1]);      
    %cur3d = fastbilat_convn(cur3d, kernelZ);
    %cur3d = fastbilat_normalize(cur3d);
    %cur3d = fastbilat_corr_gridonly(cur3d, lofilt_R, edgemethod);

    
    if (ht <= 0)
        %cur3d = fastbilat_normalize(cur3d);
        
        if ~isreal(bfilts)
            % FIXME
            % we are creating a complex representation,
            % ensure that the lowest low-level pass is also complex,
            % giving the same result as running a real-pyramid and
            % imag-pyramid and joining these.
            cur3d.gridData = cur3d.gridData + 1j .* cur3d.gridData;
        end
        
        outputs = {cur3d};
      %lo0 = fastbilat_upsample(cur3d);
      %pyr = lo0(:);
      %pind = size(lo0);

    else
        % Assume square filters:
        bfiltsz =  round(sqrt(size(bfilts,1)));

        %bands = zeros(prod(cur3d.size2d),size(bfilts,2));
        %bind = zeros(size(bfilts,2),2);

        boutput = [];
        for b = 1:size(bfilts,2)
            filt = reshape(bfilts(:,b),bfiltsz,bfiltsz);

            band3d = cur3d;
            band3d = fastbilat_corr_gridonly(band3d, filt, edgemethod);
            %band3d = fastbilat_normalize(band3d);
                
            boutput = [boutput, band3d];
            
            %band = fastbilat_upsample(band3d);
            %bands(:,b) = band(:);
            %bind(b,:)  = size(band);
        end

        % prepare for next layer
        lo3d = lo03d;
        lo3d = fastbilat_corrDn_gridonly(lo3d, lofilt, edgemethod, [2 2], [1 1]);

        noutputs = buildSpyrLevs(lo3d, ht-1, lofilt, bfilts, edgemethod, kernelZ);

        outputs = {boutput, noutputs{:}};
    end

end

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

function output = fastbilat_upsample(x3d)
    % no rounding!

    if 0
        % interpn takes rows, then cols, etc
        % i.e. size(v,1), then size(v,2), ...
        output = interpn( x3d.gridData, x3d.di, x3d.dj, x3d.map, 'linear');
    else
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

function x3d = fastbilat_convn(x3d, kernelZ)
    % convolution
    x3d.gridData = convn( x3d.gridData, kernelZ, 'same' );
    x3d.gridWeights = convn( x3d.gridWeights, kernelZ, 'same' );
end

function x3d = fastbilat_normalize(x3d)
    % normalize
    thresh = 0;
    mask = (x3d.gridWeights <= thresh);
    normalizedGrid = x3d.gridData ./ x3d.gridWeights;
    normalizedGrid( mask ) = 0; % put 0s where it's undefined
    
    x3d.gridData = normalizedGrid;
    x3d.gridWeights = cast(~mask, 'like', x3d.gridWeights);
end

function x3d = fastbilat_normalize_OLD(x3d)
    % normalize
    gridWeights = x3d.gridWeights; 
    gridWeights( gridWeights == 0 ) = -2; % avoid divide by 0, won't read there anyway
    normalizedGrid = x3d.gridData ./ gridWeights;
    normalizedGrid( gridWeights < -1 ) = 0; % put 0s where it's undefined
    
    x3d.gridData = normalizedGrid;
    x3d.gridWeights(x3d.gridWeights ~= 0) = 1;
end


function x3d = fastbilat_subsample(x3d, step, start)
    x3d.gridData = x3d.gridData(start(1):step(1):end, start(2):step(2):end, :);
    x3d.gridWeights = x3d.gridWeights(start(1):step(1):end, start(2):step(2):end, :);
    size2d = [size(x3d.gridData,1), size(x3d.gridData,2)];
    x3d.di = x3d.di(1:size2d(1), 1:size2d(2));
    x3d.dj = x3d.dj(1:size2d(1), 1:size2d(2));
    x3d.map = x3d.map(start(1):step(1):end, start(2):step(2):end);
    x3d.size2d = size2d;
end

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

function x3d = fastbilat_corrDn_gridonly(x3d, lofilt, edgemethod, step, start)
    x3d.gridData = corrDnBatch(x3d.gridData, lofilt, edgemethod, step, start);
    x3d.gridWeights = x3d.gridWeights(start(1):step(1):end, start(2):step(2):end, :);
    
    size2d = [size(x3d.gridData,1), size(x3d.gridData,2)];

    x3d.di = x3d.di(1:size2d(1), 1:size2d(2));
    x3d.dj = x3d.dj(1:size2d(1), 1:size2d(2));
    x3d.map = x3d.map(start(1):step(1):end, start(2):step(2):end);
    x3d.size2d = size2d;    
end

function DEBUG_SCAT_VIEW(x3d, layer)
    whas x3d.gridData(:,:,l) x3d.gridWeights
    sfigure(2);
    clf
    subplot(2,2,1)
    imagesc(x3d.gridData(:,:,layer))
    subplot(2,2,2)
    imagesc(x3d.gridWeights(:,:,layer))
    subplot(2,2,3)
    imagesc(x3d.gridData(:,:,layer) ./ x3d.gridWeights(:,:,layer), [0 1])
    subplot(2,2,4)
    imagesc(x3d.map);
    %imagesc(abs(x3d.map - layer) < .5)
    imagesc(fastbilat_upsample(x3d), [0 1])
    set(findobj(gcf, 'type', 'axes'), 'tag', 'x3d')
    linkaxes(findobj('tag', 'x3d'))
    drawnow;
end
