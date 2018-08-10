% Plots the k clusters using k-means
load % INSERT DATASET HERE
% Parameters for dataset
N = 50;
p = 1;  % Smoothing (max=0, min=1)
invariant2rot = 0;  % 1 for rotation invariance, 0 otherwise
k = 4;  % Number of clusters to consider


%Pre-processing: interpolation + conversion to angles
Omega = convert_dataset(X, N, p, invariant2rot);

%kmeans
max_iter = 100;
[id, centroids, residual] = circkmeans(Omega, k, max_iter);

%Plot clusters
for i=1:k
    set(gca, 'ColorOrderIndex', i)
    for j=find(id == i)
        plot(X{j}(1,:), X{j}(2,:));
        hold on;
        set(gca, 'ColorOrderIndex', i);
    end
    axis tight;
    pause;
end
