function [lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp1Filters_imag
    [lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp1Filters;
    
    n = size(bfilts,1);
    sz = floor(sqrt(n));
    sz = [sz n/sz];
    for b = 1:size(bfilts,2)
        B = reshape(bfilts(:,b), sz);
        if b <= size(bfilts,2)/2
            % hilbert2 performed over vertical direction
            HB = fftshift(hilbert2(ifftshift(B)));
        else
            % hilbert2 performed over horizontal direction
            HB = fftshift(hilbert2(ifftshift(B')))';
        end
        bfilts(:,b) = HB(:);
    end
end
