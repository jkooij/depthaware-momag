close all; clear all;

%% select the sequence
seq_nums = 1:4;

% ** Baseline composition method **
% unmagnified : reconstruct original
% standard    : standard magnification approach                [original]
% blend       : compose original+magnification at image level  [Elgharib,CVPR'15]

% run ALL magnification methods, and save videos
for seq_num = seq_nums

    % NOTE: first run magnification method will pre-process data and take a bit longer
    
    fprintf('--------------\n');
    fprintf('  SEQUENCE %d\n', seq_num);
    fprintf('--------------\n');
    
    %% motion magnification methods
    run_magnify_baselines(seq_num, 'unmagnified')      % [No motion magnification]
    run_magnify_baselines(seq_num, 'standard')         % [Wadhwa'13]
    run_magnify_baselines(seq_num, 'blend')            % [Elgharib,CVPR'15]
    run_magnify_depthaware(seq_num) % 'bilatspyr' method [Kooij,ECCV'16]

    %% create space-time plots
    run_make_spacetime_slices(seq_num, 'unmagnified')
    run_make_spacetime_slices(seq_num, 'standard')
    run_make_spacetime_slices(seq_num, 'blend')
    run_make_spacetime_slices(seq_num, 'bilatspyr')
    
end
