function [lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp1Filters_complex
    [lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp1Filters;
    [~, ~, ~, bfilts_i] = sp1Filters_imag;
    
    bfilts = bfilts + 1j .* bfilts_i;
end