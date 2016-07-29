%% Back-project the extended pyramid to a normal pyramid 
function [pyr, pind] = recon_bilatspyr(xpyr, pind)
    nbands_in_level = spyrNumBands(pind(:,1:2));
    nbands = size(pind,1);
    
    pyr = [];
    for band = 1:nbands
        midx = ceil((band-1) / nbands_in_level)+1;
        
        map = xpyr.pmap{midx};
        
        bpyr = convert_band(xpyr, pind, map, band);
        pyr = [pyr; bpyr(:)];
    end

    pind = pind(:,1:2);
end

function bpyr = convert_band(xpyr, pind, map, b)
    idxs = pyrBandIndices(pind(:,1:2), b);

    gridData = xpyr.pyr(idxs,:);
    gridData = reshape(gridData, pind(b,:));

    bpyr = map_upsample(map, gridData);
end
