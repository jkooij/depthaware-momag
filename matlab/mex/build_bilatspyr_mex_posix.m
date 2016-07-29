%% Posix

pyrtoolsmexpath = fullfile(depthaware_momag_basedir, 'matlab/external/EVM_Matlab/matlabPyrTools/MEX/');
fpath = @(n) fullfile(pyrtoolsmexpath, n);

%build corrDnBatch
cmd = 'mex';
cmd = [cmd ' -v CC="gcc" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"'];
cmd = [cmd ' -I.'];
cmd = [cmd ' -I' pyrtoolsmexpath];
cmd = [cmd ' -L' pyrtoolsmexpath];
cmd = [cmd ' corrDnBatch.c'];
cmd = [cmd ' ' fpath('convolve.c') ' ' fpath('wrap.c') ' ' fpath('edges.c')];
eval(cmd)

% build upConvBatch    
cmd = 'mex';
cmd = [cmd ' -v CC="gcc" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"'];
cmd = [cmd ' -I.'];
cmd = [cmd ' -I' pyrtoolsmexpath];
cmd = [cmd ' -L' pyrtoolsmexpath];
cmd = [cmd ' upConvBatch.c'];
cmd = [cmd ' ' fpath('convolve.c') ' ' fpath('wrap.c') ' ' fpath('edges.c')];
eval(cmd)

% build min_pdist_thresh
%  REQUIRES THE Eigen LIBRARY
cmd = 'mex';
%cmd = [cmd ' -v CC="gcc" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"'];
cmd = [cmd ' -I.'];
cmd = [cmd ' -I' pyrtoolsmexpath];
cmd = [cmd ' -L' pyrtoolsmexpath];
cmd = [cmd ' min_pdist_thresh.cpp'];
eval(cmd)

