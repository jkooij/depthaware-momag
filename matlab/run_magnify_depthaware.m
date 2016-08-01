function run_magnify_depthaware(seq_num)
% Run Depth-Aware Motion Maginification method
% seq_num
%   selects input sequence, see script_setup_sequence
%
% Author: Julian Kooij, Delft University of Technology, 2015
%   "Depth-Aware Motion Magnification", (ECCV 2016)
%

% load sequence data
script_setup_sequence

% create output directory, if it does not yet exist
mkdir_check(seq_outdir)

% Turn on to see if image and depth are aligned
%  see script_setup_sequence
DEBUG_SHOW_FOREGROUND_MASK = 1; 

cache_filepath = fullfile(seq_outdir, 'preprocessed.mat');
if ~exist(cache_filepath, 'file')

    % load the color/depth images, perform pre-processing
    script_preprocess_sequence;

    % save the pre-processed data
    save(cache_filepath, 'Ir_ybrs', 'DMrs');
    
else
    % load stored pre-processed result
    load(cache_filepath);
end

method_name = 'bilatspyr'; % name of our proposed method (used in output filenames only)

%% select the filters
filt_r = sprintf('sp%dFilters', order); % real
filt_i = sprintf('sp%dFilters_imag', order); % imag
filt_c = sprintf('sp%dFilters_complex', order); % complex
% precomp_X will cache the filters once they are loaded
precomp_r = filt_r; precomp_i = filt_i; precomp_c = filt_c;

%% compute bilateral pyramids
T = size(Ir_ybrs, 4);
xPyrs = struct('pyr', cell(1,T), 'pmap', cell(1,T), 'dcenters', cell(1,T));

for t = 1:T
    fprintf('building bilateral pyramid t = %d / %d\n', t, T);
    
    Ir_ybr = Ir_ybrs(:,:,:,t);
    DMr = DMrs(:,:,t);
    
    Ir = Ir_ybr(:,:,1);

    [xpyr, xpind, precomp_c] = build_bilatspyr( ...
        Ir, ...
        lvl, ...
        precomp_c, ...
        [], ...
        DMr, ...
        edgeMin, ...
        edgeMax, ...
        sigmaRange ...
    );   

    % use single instead of double precision to reduce memory usage
    xpyr.pyr = single(xpyr.pyr);
    
    xPyrs(t) = xpyr;
end

%% perform temporal filtering

switch TEMPORAL_LOWPASS_FILTER_METHOD
    case 'center-time',
        % take an intermediate time-point as reference
        t = ceil(numel(xPyrs)/2);
        m_phi = angle(xPyrs(t).pyr);
        m_mag = abs(xPyrs(t).pyr);

        [pyr, pind] = recon_bilatspyr(xPyrs(t), xpind);

    case 'mean',
        % --------------------------------------------------------------
        % NOTE: mean filter can take A LOT of memory, time to compute
        % --------------------------------------------------------------
        clear pyrs3d phis3d mags3d
        pyrs3d = cat(3, xPyrs.pyr);

        % compute temporal filters
        phis3d = angle(pyrs3d);
        %phis3d = norm_angle_series(phis3d, 3);
        m_phi = circ_mean(phis3d,[],3);
        clear phis3d

        mags3d = abs(pyrs3d);
        clear pyrs3d
        %m_phi = mean(phis3d,3);
        m_mag = mean(mags3d,3);
        clear mags3d

    otherwise,
        error('unknown TEMPORAL_LOWPASS_FILTER_METHOD value');
end

%% reconstruct magnified

depth_to_dfactor = @(DM) exp(-((DM - mag_depth_mean)/mag_depth_stddev).^2 / 2);

%setminmax = @(x, minx, maxx) min(max(mod(x+pi/2,2*pi)-pi/2, minx), maxx);
setminmax = @(x, minx, maxx) x; % NO scale limiting

D = size(xPyrs(1).pyr,2); % pyramid depth
dlayers = 1:D; % depth layers to magnify
dfactor = depth_to_dfactor(xPyrs(1).dcenters(dlayers));

clear M_frames
for t = 1:T
    %%
    fprintf('magnifying frame t = %d / %d\n', t, T);
    
    xpyr = xPyrs(t);
    DMr = DMrs(:,:,t);

    % magnify given bands
    pyr = xpyr.pyr(:,dlayers);
    pind = xpind(:,1:2);

    phi = angle(pyr);
    mag = abs(pyr);

    phi = bsxfun(@times, setminmax(phi - m_phi(:,dlayers),-pi/2, pi/2), (dfactor .* (phase_factor-1))+1) + m_phi(:,dlayers);
    pyr2 = mag .* exp(1j * phi);

    bDMr = [];
    for b = 1:size(pind,1);
        %% select depth map of correct size
        bsize = pind(b,:);
        if any(bsize ~= size(bDMr))
            bDMr = imresize(DMr, bsize);
        end

        % create magnification mask
        %bmask = exp(-((bDMr - mag_depth_mean)/mag_depth_stddev).^2 / 2);
        bmask = depth_to_dfactor(bDMr);
        %bmask(:) = 1;

        bmask = bmask(:);
        bmask = repmat(bmask, [1 D]);

        % compose band with magnified result
        bidxs = pyrBandIndices(pind, b);
        pyr(bidxs,:) = bmask .* pyr2(bidxs,:) + (1 - bmask) .* pyr(bidxs,:);
    end

    % DEBUG
    %pyr(:) = 0;
    %pyr(320097) = 1000 * exp(1j * 2*pi*14e-1);

    xpyr.pyr(:,dlayers) = pyr;

    %% reconstruct each layer
    I_rs = [];
    for d2 = 1:D
        pyr = xpyr.pyr(:,d2);

        if 0 && ismember(d2, dlayers)
            pyr = double(pyr);
            I_r = reconSpyr(real(pyr), pind, filt_r);
            I_i = reconSpyr(imag(pyr), pind, filt_i);
            I_r = (I_r + I_i)/2;
        else
            pyr = real(pyr);
            pyr = double(pyr);
            I_r = reconSpyr(pyr, pind, filt_r);
        end

        I_rs = cat(3, I_rs, I_r);
        %imagesc(I_r, [0 1]);
        %drawnow
    end

    %% make image

    % resample
    pmap = xpyr.pmap{1};
    %pmap(:) = 6; % DEBUG: visualize specific layer
    I_r = map_upsample(pmap, I_rs);

    % slightly more clear output
    I_r = I_r .* OUTPUT_MULTIPLY_CONTRAST;
    I_r = I_r + OUTPUT_INCREASE_CONTRAST;

    Ir_ybr = Ir_ybrs(:,:,:,t);
    Ir_ybr(:,:,1) = I_r;
    Ir_rgb = ycbcr2rgb(double(Ir_ybr));

    sfigure(3);
    cla
    if 1
        % just show/grab output
        imagesc(Ir_rgb, [0 1]);
        axis image

        M_frames(t) = im2frame(Ir_rgb);
        %M_frames(t) = getframe(gcf);
    else
        % compare output to original
        subplot(1,2,1)
        imagesc(ycbcr2rgb(double(Ir_ybrs(:,:,:,t))), [0 1]);
        axis image

        subplot(1,2,2);
        imagesc(Ir_rgb, [0 1]);
        axis image

        M_frames(t) = getframe(gcf);
    end

    drawnow
    
end

%% save output

% save the pre-processed data
output_filename = sprintf('%s_%s_f%.1f.mat', seq_name, method_name, phase_factor);
output_filepath = fullfile(seq_outdir, output_filename);
fprintf('saving %s ...\n', output_filename);
save(output_filepath, 'M_frames', 'method_name', 'phase_factor', 'seq_num', 'seq_name');

%%
output_filename = sprintf('%s_%s_f%.1f.avi', seq_name, method_name, phase_factor);
output_filepath = fullfile(seq_outdir, output_filename);
fprintf('saving %s ...\n', output_filename);
movie2avi(M_frames, output_filepath, 'Fps', 30, 'Compression', 'none');


%% replay

sfigure(3);
clf;
movie(M_frames, 1, 30)
