function run_make_spacetime_slices(seq_num, method_name)

% load sequence data
script_setup_sequence

switch seq_name
    case 'sequence1',
        show_t = 50;
        slice_ys = 50:150; slice_xs = 70;
        
    case 'sequence2',
        show_t = 112;
        slice_xs = 105; slice_ys = 0 + [20:100];
        
    case 'sequence3',
        show_t = 23;
        slice_xs = 1:100; slice_ys = 25;

    case 'sequence4',
        show_t = 60;
        slice_xs = 20; slice_ys = 0 + [5:80];

    otherwise,
        error('please define slice for %s\n', seq_name);
end

% load the stored results
output_filename = sprintf('%s_%s_f%.1f.mat', seq_name, method_name, phase_factor);
output_filepath = fullfile(seq_outdir, output_filename);
fprintf('load %s ...\n', output_filename);
load(output_filepath);

%% collect temporal slice
slices = [];
T = numel(M_frames);
for t = 1:T
    Ir = M_frames(t).cdata;
    slice = Ir(slice_ys, slice_xs, :);
    slice = reshape(slice, 1, [], 3);
    slices = [slices; slice];
end

% create annotated image
dx = 0; dy = 0;
if numel(slice_xs) == 1;
    dx = [-1:1];
    slices = permute(slices, [2 1 3]);
end
if numel(slice_ys) == 1; dy = [-1:1]; end

Ir_annot = M_frames(show_t).cdata;
Ir_annot(slice_ys + dy, slice_xs + dx, 1) = 255;
Ir_annot(slice_ys + dy, slice_xs + dx, 2) = 0;
Ir_annot(slice_ys + dy, slice_xs + dx, 3) = 0;

%%
sfigure(2);
clf;
subplot(2,1,1);
imagesc(Ir_annot);
axis image

subplot(2,1,2);
imagesc(slices);

%%
if 1
    %% save images
    mkdir_check(seq_outdir);

    idstring = sprintf('%s_f%.1f', method_name, phase_factor);
    
    out_filename = sprintf('timeslice_%s.png', idstring);
    out_filepath = fullfile(seq_outdir, out_filename);
    fprintf('saving %s ...\n', out_filename);
    imwrite(slices, out_filepath);

    out_filename = sprintf('timeslice_annot_%s.png', idstring);
    out_filepath = fullfile(seq_outdir, out_filename);
    fprintf('saving %s ...\n', out_filename);
    imwrite(Ir_annot, out_filepath);
end
