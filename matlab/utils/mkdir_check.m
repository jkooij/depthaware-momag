function mkdir_check(dirpath)
    if ~isdir(dirpath)
        warning('creating directory %s', dirpath);
        mkdir(dirpath);
    end
end