basedir = fileparts(mfilename('fullpath'));

% add files to path
addpath(basedir);
addpath(fullfile(basedir, 'bilatspyr'));
addpath(fullfile(basedir, 'filters'));
addpath(fullfile(basedir, 'mex'));
addpath(fullfile(basedir, 'utils'));