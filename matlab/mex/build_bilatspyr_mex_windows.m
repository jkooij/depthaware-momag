%% Windows

pyrtoolsmexpath = fullfile(depthaware_momag_basedir, 'matlab/external/EVM_Matlab/matlabPyrTools/MEX/');
fpath = @(name) fullfile(pyrtoolsmexpath, name);

EigenDir = 'D:\libs\Eigen';

% cd to directory of this script
curdir = pwd;
compdir = fileparts(mfilename('fullpath'));
cd(compdir)

fprintf('compile directory: %s\n', compdir);
fprintf('matlabPyrTools/mex directory: %s\n', pyrtoolsmexpath);

% corrDnBatch.c
mex_cmd = ['mex -I' pyrtoolsmexpath ' -L' pyrtoolsmexpath ' corrDnBatch.c ' fpath('convolve.c') ' ' fpath('wrap.c') ' ' fpath('edges.c')];
mex_cmd = [mex_cmd ' COMPFLAGS="$COMPFLAGS /openmp"']; % Visual C++ OpenMP
fprintf('MEX command: %s\n', mex_cmd);
eval(mex_cmd);

% upConvBatch.c
mex_cmd = ['mex -I' EigenDir ' min_pdist_thresh.c'];
fprintf('MEX command: %s\n', mex_cmd);
eval(mex_cmd);

% restore directory
cd(curdir);