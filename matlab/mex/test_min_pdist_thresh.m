build_mex

%%
X = rand(4,3);
Y = rand(5, 3);
%Y(2:2:end,:) = Y(1:2:end,:)+.1;

X(:,3) = 0; Y(:,3) = 0;

%% first time is slow
[knn_classes, knn_dist] = min_pdist_thresh(X, Y, 0.1);

%%
widths = ones(1, size(X,1)) * 10e-1;
widths(1) = 4e-1;
widths(2) = 7e-1;

tic
for r = 1:10
    [knn_classes, knn_dist] = min_pdist_thresh(X, Y, widths);
end
toc
knn_classes'
knn_dist'

if 0
    %% comparison
    tic
    for r = 1:10
        [knn_classes, knn_dist] = knnsearch(X, Y);
    end
    toc
    knn_classes'
    knn_dist'
end

% visualize
yj = 3;
xj = knn_classes(yj);
if isnan(xj);
    xj = [];
    dj = NaN;
    wj = NaN;
else
    dj = knn_dist(yj);
    wj = widths(xj);
end

sfigure(4);
cla;
plot(X(:,1), X(:,2), '.');
hold all;
plot(Y(yj,1), Y(yj,2), 'x');
plot_circle(Y(yj,:), dj, 'k--');
plot(X(xj,1), X(xj,2), 'd'); % found closest point
plot_circle(Y(yj,:), wj, 'r--');

axis equal
xlabel('x'); ylabel('y'); zlabel('z');
grid on
axis([0 1 0 1])