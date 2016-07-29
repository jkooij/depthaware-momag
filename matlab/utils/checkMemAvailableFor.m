function checkMemAvailableFor(size, type)
    if nargin < 2; type = 'double'; end
    
    [numBytesNeeded, type] = bytesizeof(type);
    numBytesNeeded = numBytesNeeded * prod(size(:));

    % See http://stackoverflow.com/a/14144392/218682
    maxMemFrac = 0.8; % use at most 80% of the available memory

    %% read available memory
    try
        [~,memStats] = memory;
    catch ME
        if strcmp(ME.identifier, 'MATLAB:memory:unsupported')
            % on Linux, the MEMORY function is not supported,
            % just continue
            return
        else
            rethrow(ME)           
        end
    end
     
    %%
    numBytesAvailable = memStats.PhysicalMemory.Available;

    if numBytesNeeded > numBytesAvailable * maxMemFrac
        sizestr = [sprintf('%d', size(1)) sprintf(' x %d', size(2:end))];
        error('MYSIM:OUTOFMEMORY','too much memory needed for %s array of size [%s] (%.2f MB > %.2f MB x %d%% available)', ...
            type, sizestr, numBytesNeeded / 1024.^2, numBytesAvailable / 1024.^2, maxMemFrac*100);
    end    

end

function [bytes, classname] = bytesizeof(type)
    if ischar(type)
        % convert string into a dummy variable to measure with WHOS
        type = lower(type);
        if any(strfind(type, 'complex ') == 1)
            % special case
            type = type(9:end);
            type = cast(1, type);
            type = complex(type);
        else
            % type is a string like 'double', 'uint8'
            type = cast(1, type);
        end
    end
    
    dinfo = whos('type');
    bytes = dinfo.bytes / prod(dinfo.size);
    classname = dinfo.class;
    if dinfo.complex
        classname = ['complex ' classname];
    end
end