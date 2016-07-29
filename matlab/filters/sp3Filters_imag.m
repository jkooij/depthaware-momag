% Steerable pyramid filters.  Transform described  in:
%
% @INPROCEEDINGS{Simoncelli95b,
%	TITLE = "The Steerable Pyramid: A Flexible Architecture for
%		 Multi-Scale Derivative Computation",
%	AUTHOR = "E P Simoncelli and W T Freeman",
%	BOOKTITLE = "Second Int'l Conf on Image Processing",
%	ADDRESS = "Washington, DC", MONTH = "October", YEAR = 1995 }
%
% Filter kernel design described in:
%
%@INPROCEEDINGS{Karasaridis96,
%	TITLE = "A Filter Design Technique for 
%		Steerable Pyramid Image Transforms",
%	AUTHOR = "A Karasaridis and E P Simoncelli",
%	BOOKTITLE = "ICASSP",	ADDRESS = "Atlanta, GA",
%	MONTH = "May",	YEAR = 1996 }

% Eero Simoncelli, 6/96.

function [lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp3Filters_imag
    [lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp3Filters;
    
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
