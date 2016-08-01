function basedir = depthaware_momag_basedir
    basedir = fileparts(mfilename('fullpath'));
    basedir = fileparts(basedir); % strip matlab/ part
end