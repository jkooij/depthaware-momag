buildExtraPyrTools

%% load test image
I = imread('cameraman.tif');
I = im2double(I);
I = imresize(I, size(I)/2);

sfigure(10);
imagesc(I)
axis image

%% multiple images (i.e. a batch)
Is = cat(3, I, I', I, I');

%% get a nice filter
%[lo0filt,hi0filt,lofilt,bfilts,mtx,harmonics] = sp3Filters();
%B = reshape(bfilts(:,1), [9 9]);
B = [
   -0.0008    0.0039    0.0013    0.0007         0   -0.0007   -0.0013   -0.0039    0.0008
    0.0044    0.0045   -0.0038   -0.0004         0    0.0004    0.0038   -0.0045   -0.0044
    0.0123   -0.0059    0.0083   -0.0225         0    0.0225   -0.0083    0.0059   -0.0123
    0.0140   -0.0029    0.0394   -0.1106         0    0.1106   -0.0394    0.0029   -0.0140
    0.0142    0.0085    0.0536   -0.1768         0    0.1768   -0.0536   -0.0085   -0.0142
    0.0140   -0.0029    0.0394   -0.1106         0    0.1106   -0.0394    0.0029   -0.0140
    0.0123   -0.0059    0.0083   -0.0225         0    0.0225   -0.0083    0.0059   -0.0123
    0.0044    0.0045   -0.0038   -0.0004         0    0.0004    0.0038   -0.0045   -0.0044
   -0.0008    0.0039    0.0013    0.0007         0   -0.0007   -0.0013   -0.0039    0.0008
];


%%
tic
%O = corrDn(I, B', 'reflect1', [2 2]);
%O = corrDnBatch(Is, B, 'reflect1', [2 2]);
O = upConvBatch(Is, B, 'reflect1');
%O = convn( Is, B, 'same' );
%O = upConvBatch(Is, B);
toc

whas O

sfigure(10);
imagesc(O(:,:,1));
axis image