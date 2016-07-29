function Outs = batch_imresize(Imgs, varargin)
    if numel(size(Imgs)) > 3
        N = size(Imgs, 4);
        Outs = cell(1,N);
        for j = 1:N
            Outs{j} = imresize(Imgs(:,:,:,j), varargin{:});
        end
        Outs = cat(4, Outs{:});
    else
        N = size(Imgs, 3);
        Outs = cell(1,N);
        for j = 1:N
            Outs{j} = imresize(Imgs(:,:,j), varargin{:});
        end
        Outs = cat(3, Outs{:});
    end
end