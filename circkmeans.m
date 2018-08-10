function [ idx, centroids, residual ] = circkmeans( X, k, max_iter )
%CIRCKMEANS K-means adapted for circular data
%   X is n x m, where n is the number of trajectories and m is the number of features
%   idx contains the cluster assignments: idx(1) = 2 means trajectory 1 belongs to cluster 2
%   centroids is self explainatory
%   residual is the residual error for the k clusters
X = X';
centroids = init_centroids(X, k);
old_idx = zeros(size(X, 1), 1);
conv_iter = 0;

for i=1:max_iter
    idx = findclosest(X, centroids);
    centroids = comp_centroids(X, idx, k);
    
    if sum(idx == old_idx) == size(X, 1) %Solution converged
        conv_iter = conv_iter + 1;
        if conv_iter == 3 %3 consecutive iterations with the same assignments
            break;
        end
    else
        conv_iter = 0;
    end
    old_idx = idx;
end

centroids = centroids';
idx = idx';
residual = sum(sum(1-cos(X'-centroids(:, idx))));

end

function [ dist, min_id ] = comp_dist(centroids, x)
%COMP_DIST Computes the distance from x to the centroids
%   Detailed explanation goes here
% min_norm = inf;
% for add=0:0.1:2*pi
%     d = sin(bsxfun(@minus, centroids, x+add)/2);
%     [d, id] = min(sum(d.^2, 2));
%     if d < min_norm
%         dist = d;
%         min_id = id;
%         min_norm = d;
%     end
% end
dist = sin(bsxfun(@minus, centroids, x)/2);
%dist = min(bsxfun(@minus, centroids, x), 2*pi-bsxfun(@minus, centroids, x));
[dist, min_id] = min(sum(dist.^2, 2));

end

function [ centroids ] = init_centroids( X, k )
%INIT_CENTROIDS Initializes the centroids that are to be used by kmeans
%   Detailed explanation goes here

m = size(X, 1);
centroid_points = [randi(m)]; %1st centroid chosen uniformely at random
D = zeros(m, 1);

% Carefully chosen centroids, give the best results.
% centroids = X([1, 13, 24, 31, 44, 76, 86, 93, 100, 108, 120, 128], :);
% return;

for c=2:k
    for i=1:m
        D(i) = comp_dist(X(centroid_points, :), X(i, :));
    end
    centroid_points = [centroid_points, randsample(m, 1, true, D)];
end
centroids = X(centroid_points, :);

end

function [ idx ] = findclosest( X, centroids )
%FINDCLOSEST Computes centroid membership for every example
%   Detailed explanation goes here
idx = zeros(size(X, 1), 1);
for i=1:size(X, 1)
    [~, idx(i)] = comp_dist(centroids, X(i, :));
end

end

function [ centroids ] = comp_centroids( X, idx, k )
%COMP_CENTROIDS Computes the new centroids
%   Detailed explanation goes here
centroids = zeros(k, size(X, 2));
for i=1:k
    id = find(idx == i);
    centroids(i, :) = angle(sum(exp(1j*X(id, :)), 1));
end

end
