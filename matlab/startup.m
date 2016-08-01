basedir = fileparts(mfilename('fullpath'));

% add files to path
addpath(basedir);
addpath(fullfile(basedir, 'bilatspyr'));
addpath(fullfile(basedir, 'filters'));
addpath(fullfile(basedir, 'mex'));
addpath(fullfile(basedir, 'utils'));

% Ensure that Eero P. Simoncelli's matlabPyrTools are added to the path
%   See https://github.com/LabForComputationalVision/matlabPyrTools
%   Of course, you can also add them to your global path yourself
addpath(fullfile(basedir, 'external', 'matlabPyrTools'));
addpath(fullfile(basedir, 'external', 'matlabPyrTools', 'MEX'));