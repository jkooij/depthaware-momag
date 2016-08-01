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
    im3d = fastbilat_tdistconv(im3d); % DEBUG: t-distribution instead of Gaussian for more stable tales
    %im3d = fastbilat_gaussconv(im3d);
    
    % perform 1D depth smoothing kernel
    im3d = fastbilat_convn(im3d, kernelZ);
    
    % -- step 2 ---
    % element-wise division by weight total
    im3d = fastbilat_normalize(im3d);

    % -- step 3 ---
    % replace data again
    im3d.gridData(mask) = origD(mask);
    im3d.gridWeights(mask) = origW(mask);
    %% -- step 4 --
    
    % intial lowpass/highpass filters
    lo03d = fastbilat_corrDn(im3d, lo0filt, edgemethod, [1 1], [1 1]);
    hi03d = fastbilat_corrDn(im3d, hi0filt, edgemethod, [1 1], [1 1]);
    
    %% build level recursively
    outputs = build_bilatspyr_levs(lo03d, ht, lofilt, bfilts, edgemethod, kernelZ);
    outputs = {hi03d, outputs{:}};
    
    % remove all the nested cell arrays, simplify the datastructure
    [pyr, pind] = out_to_pyr_pind_map(outputs);
    pyr.dcenters = dcenters;
end

function outputs = build_bilatspyr_levs(lo03d, ht, lofilt, bfilts, edgemethod, kernelZ)
    % OUTPUTS = build_bilatspyr_levs(LOIM, HEIGHT, LOFILT, BFILTS, EDGES, KERNELZ)
    % Recursive function for constructing levels of a bilateral steerable pyramid. 
    %
    % Adapted from buildSpyrLevs by Eero Simoncelli, 6/96.
    %

    cur3d = lo03d;
    
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

    else
        % Assume square filters:
        bfiltsz =  round(sqrt(size(bfilts,1)));

        boutput = [];
        for b = 1:size(bfilts,2)
            filt = reshape(bfilts(:,b),bfiltsz,bfiltsz);

            band3d = cur3d;
            band3d = fastbilat_corr_gridonly(band3d, filt, edgemethod);
                
            boutput = [boutput, band3d];
        end

        % prepare for next layer
        lo3d = lo03d;
        lo3d = fastbilat_corrDn_gridonly(lo3d, lofilt, edgemethod, [2 2], [1 1]);

        noutputs = build_bilatspyr_levs(lo3d, ht-1, lofilt, bfilts, edgemethod, kernelZ);

        outputs = {boutput, noutputs{:}};
    end

end

function [pyr, pind] = out_to_pyr_pind_map(outputs)
    % Remove all the nested cell arrays, simplify the datastructure

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
