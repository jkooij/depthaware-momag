function scat_debug_view(x3d, layer)
    sfigure(2);
    clf
    subplot(2,2,1)
    imagesc(x3d.gridData(:,:,layer))
    subplot(2,2,2)
    imagesc(x3d.gridWeights(:,:,layer))
    subplot(2,2,3)
    imagesc(x3d.gridData(:,:,layer) ./ x3d.gridWeights(:,:,layer), [0 1])
    subplot(2,2,4)
    imagesc(x3d.map);
    %imagesc(abs(x3d.map - layer) < .5)
    imagesc(fastbilat_upsample(x3d), [0 1])
    set(findobj(gcf, 'type', 'axes'), 'tag', 'x3d')
    linkaxes(findobj('tag', 'x3d'))
    drawnow;
end
