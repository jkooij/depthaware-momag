function alphas = norm_angle_series(alphas, dim)
% Similar (but faster) than UNWRAP, fix jumps in angles in range [-pi, pi]
%    ps = unwrap(ps, [], dim);

    checkMemAvailableFor(size(alphas), alphas)

    switch dim
        case 1,
            alphas = cumsum([alphas(1,:); circ_dist(alphas(2:end,:), alphas(1:end-1,:))], 1);
        
        case 2,
            alphas = cumsum([alphas(:,1), circ_dist(alphas(:,2:end), alphas(:,1:end-1))], 2);
    
        case 3,
            alphas = cumsum(cat(3,alphas(:,:,1), circ_dist(alphas(:,:,2:end), alphas(:,:,1:end-1))), 3);
    
        otherwise,
            % NOTE: this will use a lot more memory!
            alphas = shiftdim(alphas, dim-1);
            alphas = cumsum([alphas(1,:); circ_dist(alphas(2:end,:), alphas(1:end-1,:))], 1);

            ndims = numel(size(alphas));
            alphas = shiftdim(alphas, ndims-(dim-1));
    end
end