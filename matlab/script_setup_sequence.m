% select the sequence
if ~exist('seq_num', 'var') || isempty(seq_num);
    seq_num = 1;
end
seq_name = sprintf('sequence%d', seq_num);

fprintf('--- Processing "%s" ---\n', seq_name);

% location where data can be found
base_datadir = fullfile(depthaware_momag_basedir, 'data');

% location to store (intermediate) results
base_outdir = fullfile(depthaware_momag_basedir, 'output');

color_filename = 'kinect.mj2';
depth_filename = 'kinect_map.mj2';

seq_datadir = fullfile(base_datadir, seq_name);
seq_outdir = fullfile(base_outdir, seq_name);
color_filepath = fullfile(seq_datadir, color_filename);
depth_filepath = fullfile(seq_datadir, depth_filename);

%% -- settings --

% body-part computation?
boundary_thresh = 2900;
body_dist_thresh = 220; % mm        
%body_dist_thresh = 180; % mm        

% pyramid settings
lvl = 3;
order = 1;

TEMPORAL_LOWPASS_FILTER_METHOD = 'mean';

OUTPUT_INCREASE_CONTRAST = 0;
OUTPUT_MULTIPLY_CONTRAST = 1;        %%http://nl.mathworks.com/help/matlab/ref/ispc.html


switch seq_name
    case 'sequence1',
        % pyramid settings
        lvl = 3;
        order = 3; % <-- !
        
        % [Preprocessing] input-image downscaling
        scale = 2;
        rect = [934 185 359 346];
        
        % bilateral depth settings
        edgeMin = 1900;
        edgeMax = 2900;
        sigmaRange = 100e0;

        % Magnification settings
        mag_depth_mean = 2300;
        mag_depth_stddev = 200;        
        mag_factor = 1; % magnitude
        phase_factor = 10; % phase
        
        % [Preprocessing] depth-map to image offset
        dm_offset = [4 1];        
        depthmap_med_fsize = 15;
        depthmap_ord_fsize = 13;
        
        
        OUTPUT_INCREASE_CONTRAST = 0;
        TEMPORAL_LOWPASS_FILTER_METHOD = 'center-time';
        
    case 'sequence2',  
        % pyramid settings
        lvl = 3;        
        order = 1;
        
        % [Preprocessing] input-image downscaling
        scale = 2;
        rect = [879 199 332 253];
        
        % Magnification settings
        mag_depth_mean = 2440;
        mag_depth_stddev = 200;
        mag_factor = 1;
        phase_factor = 3;
                
        % bilateral depth settings
        edgeMin = 2300; % vid has lots of depth noise
        edgeMax = 3000;
        sigmaRange = 100e0;
        
        % [Preprocessing] depth-map to image offset
        dm_offset = [3 1];        
        depthmap_med_fsize = 15;
        depthmap_ord_fsize = [13 10];
        
        OUTPUT_INCREASE_CONTRAST = 0;
        
    case 'sequence3',
        % pyramid settings
        lvl = 3;
        order = 1;
        
        % [Preprocessing] input-image downscaling
        scale = 2;
        rect = [884 352 211 201];
        
        % Magnification settings
        mag_depth_mean = 2200;
        mag_depth_stddev = 100;
        mag_factor = 1;
        phase_factor = 5;
        
        % bilateral depth settings
        edgeMin = 1700;
        edgeMax = 2700;
        sigmaRange = 100e0;
        
        % [Preprocessing] depth-map to image offset
        dm_offset = [0 2];        
        depthmap_med_fsize = 20;
        depthmap_ord_fsize = 10;
        
        OUTPUT_INCREASE_CONTRAST = 0;
        OUTPUT_MULTIPLY_CONTRAST = 2;
        
       
    case 'sequence4',
        % pyramid settings
        lvl = 3;
        order = 1;

        % [Preprocessing] input-image downscaling
        scale = 2;
        rect = [1103 269 199 170];
        
        % Magnification settings
        mag_depth_mean = 2500;
        mag_depth_stddev = 100;
        mag_factor = 1;
        phase_factor = 10;
        
        % bilateral depth settings
        edgeMin = 2300;
        edgeMax = 3000;
        sigmaRange = 100e0;
        
        % [pre-processing]
        dm_offset = [3 4];        
        depthmap_med_fsize = 15;
        depthmap_ord_fsize = [15 15];
        
        OUTPUT_INCREASE_CONTRAST = 0;
        OUTPUT_MULTIPLY_CONTRAST = 2;
        
    otherwise,
        warning('unidentified selection');
end
