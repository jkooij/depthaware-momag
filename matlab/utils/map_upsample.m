function [output, di, dj] = map_upsample(map, gridData, di, dj)
    % no rounding!

    [height, width, depth] = size(gridData);
    if nargin < 3 || isempty(di)
        [ dj, di ] = meshgrid( 1 : width, 1 : height);
    end
    
    mapl = floor(map);
    maph = ceil(map);
    alpha = map - mapl;

    mapl(mapl < 1) = 1; maph(maph < 1) = 1;
    mapl(mapl > depth) = depth; maph(maph > depth) = depth;

    indl = sub2ind(size(gridData), di, dj, mapl);
    indh = sub2ind(size(gridData), di, dj, maph);

    output = gridData(indl) .* (1-alpha) + gridData(indh) .* alpha;
end
