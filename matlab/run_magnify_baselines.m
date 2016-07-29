close all; clear all;

%% select the sequence
seq_num = 3;

% ** Baseline composition method **
% unmagnified : reconstruct original
% standard    : standard magnification approach                [original]
% blend       : compose original+magnification at image level  [Elgharib,CVPR'15]

%method_name = 'unmagnified'; % [No motion magnification]
%method_name = 'standard';    % [Wadhwa'13]
method_name = 'blend';       % [Elgharib,CVPR'15]

%% load sequence data
script_setup_sequence

% create output directory, if it does not yet exist
mkdir_check(seq_outdir)

% Turn on to see if image and depth are aligned
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

%% select the filters
filt_r = sprintf('sp%dFilters', order); % real
filt_i = sprintf('sp%dFilters_imag', order); % imag

% pyramid settings
USE_SCF = 0; % use Steerable-Complex pyramid?

%% compute pyramids
T = size(Ir_ybrs, 4);
Pyrs = [];

for t = 1:T
    fprintf('building pyramid t = %d / %d\n', t, T);
    
    Ir_ybr = Ir_ybrs(:,:,:,t);
    Ir = Ir_ybr(:,:,1);

    if USE_SCF
        % computes complex pyramid line
        [pyr, pind] = buildSCFpyr(Ir, lvl, order);
    else
        % apply real and imaginary filters
        [pyr, pind] = buildSpyr(double(Ir), lvl, filt_r);
        [pyr_i, ~] = buildSpyr(double(Ir), lvl, filt_i);
        pyr = single(pyr + 1j * pyr_i);
    end
        
    Pyrs = [Pyrs, pyr];
end

%% perform temporal filtering

switch TEMPORAL_LOWPASS_FILTER_METHOD
    case 'center-time',
        % take an intermediate time-point as reference
        t = ceil(size(Pyrs,2)/2);
        m_phi = angle(Pyrs(:,t));
        m_mag = abs(Pyrs(:,t));

    case 'mean',
        %% compute temporal filters
        phis = angle(Pyrs);
        mags = abs(Pyrs);

        phis = norm_angle_series(phis, 2);
        m_phi = circ_mean(phis,[],2);
        %m_phi = mean(phis,2);
        m_mag = mean(mags,2);

    otherwise,
        error('unknown TEMPORAL_LOWPASS_FILTER_METHOD value');
end

%% reconstruct magnified
depth_to_dfactor = @(DM) exp(-((DM - mag_depth_mean)/mag_depth_stddev).^2 / 2);

%setminmax = @(x, minx, maxx) min(max(mod(x+pi/2,2*pi)-pi/2, minx), maxx);
setminmax = @(x, minx, maxx) x; % NO scale limiting

clear M_frames
for t = 1:T
    %%
    fprintf('magnifying frame t = %d / %d\n', t, T);
    
    pyr = Pyrs(:,t);
    DMr = DMrs(:,:,t);
    
    % select baseline magnification method
    switch method_name
        case 'unmagnified'
            % reconstruct original            
            if USE_SCF
                I_r = reconSCFpyr(pyr, pind);
            else
                I_r = reconSpyr(double(real(pyr)), pind, filt_r);
            end
            
        case 'standard'
            % standard magnification approach

            % create fully magnified pyramid
            phi = angle(pyr); mag = abs(pyr);
            phi = setminmax(phi - m_phi,-pi/2, pi/2)*phase_factor + m_phi;
            mag = (mag - m_mag)*mag_factor + m_mag;
            pyr = mag .* exp(1j * phi);
            
            % construct image
            if USE_SCF
                I_r = reconSCFpyr(pyr, pind);
            else
                I_r = reconSpyr(double(real(pyr)), pind, filt_r);
            end

        case 'blend',
            % compose original+magnification at image level
            
            % create fully magnified pyramid
            phi = angle(pyr); mag = abs(pyr);
            phi = setminmax(phi - m_phi,-pi/2, pi/2)*phase_factor + m_phi;
            mag = (mag - m_mag)*mag_factor + m_mag;
            pyr2 = mag .* exp(1j * phi);
            
            I_r2 = reconSpyr(double(real(pyr2)), pind, filt_r); % magnified image
            I_r = reconSpyr(double(real(pyr)), pind, filt_r); % original image
            
            bmask = depth_to_dfactor(DMr);
            I_r = bmask .* I_r2 + (1 - bmask) .* I_r;
    end

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
