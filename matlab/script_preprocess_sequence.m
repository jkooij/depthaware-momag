% ------------------------------------
% run SCRIPT_SETUP_SEQUENCE first
% ------------------------------------

if ~exist('DEBUG_SHOW_FOREGROUND_MASK', 'var');
    DEBUG_SHOW_FOREGROUND_MASK = 1;
end

% target size
selected_rect = rectregion(rect(1), rect(2), rect(3), rect(4));
bsize = ([selected_rect.height, selected_rect.width] + [1 1]) / scale;

Ir_ybrs = [];
DMrs = [];
t = 0;

% open video file readers
vC = VideoReader(color_filepath); % color video
vD = VideoReader(depth_filepath); % depth video

%% load all frames, compute bilateral pyramids

while hasFrame(vC)
    %%
    t = t + 1;
    fprintf('.');

    % load color and depth frame
    I = readFrame(vC);
    DM = readFrame(vD);
    
    % color image
    Ir_rgb = imcrop(I, selected_rect.x0y0wh + [0 0 -1 -1]);
    Ir_rgb = imresize(Ir_rgb, bsize, 'bilinear');
    Ir_rgb = im2double(Ir_rgb);
    Ir_ybr = rgb2ycbcr(Ir_rgb);
    Ir_ybr = single(Ir_ybr);    
            
    % depth image
    DMr = batch_imcrop(DM, selected_rect.x0y0wh + [0 0 -1 -1] + [-dm_offset 0 0]);
    DMr = DMr(:,:,1); % keep only z-world coordinate (actual depth), express it in mm
    DMr = im2double(DMr) * 4000; % 0-1 range = [0mm 4000mm]

    % [Preprocessing] median filter
    DMr(:,:,1) = medfilt2(DMr(:,:,1), [1 1]*depthmap_med_fsize);
    % [Preprocessing] orginal (max) filter
    DMr = ordfilt2(DMr, depthmap_ord_fsize(1), ones(1,depthmap_ord_fsize(1))); % max filter
    DMr = ordfilt2(DMr, depthmap_ord_fsize(end), ones(depthmap_ord_fsize(end),1)); % max filter
    
    DMr(DMr < edgeMin) = edgeMin; 
    DMr = batch_imresize(DMr, bsize, 'nearest'); % nearest = don't interpolate depth values
    
    if DEBUG_SHOW_FOREGROUND_MASK
        %% DEBUG: visualize image within depth region
        mask = (DMr > (mag_depth_mean-mag_depth_stddev)) & ...
               (DMr < (mag_depth_mean+mag_depth_stddev));
        sfigure(3);
        clf;
        imagesc(bsxfun(@times, Ir_rgb, mask));
        axis image
        drawnow
    end

    %% store data
    Ir_ybrs = cat(4, Ir_ybrs, Ir_ybr);    
    DMrs = cat(3, DMrs, DMr);
    
end

T = t;
fprintf('\n', T);
fprintf('found %d frames\n', T);