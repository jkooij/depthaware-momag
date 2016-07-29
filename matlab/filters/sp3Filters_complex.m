function [lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp3Filters_complex
    [lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp3Filters;
    [~, ~, ~, bfilts_i] = sp3Filters_imag;
    
    bfilts = bfilts + 1j .* bfilts_i;
end