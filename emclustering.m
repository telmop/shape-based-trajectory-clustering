% Von Mises Mixtures
load % INSERT DATASET HERE
% Parameters for dataset
N = 50;
p = 1;  % Smoothing (max=0, min=1)
invariant2rot = 0;  % 1 for rotation invariance, 0 otherwise
k = 4;  % Number of clusters to consider

%Pre-processing: interpolation + conversion to angles
Omega = convert_dataset(X, N, p, invariant2rot);

%FMM
max_iter = 100;
multiple_ks = 1;
th = 1e-4;
R0 = -5e-5;
c = 5e-5;
[ w, alpha, mu, kappa, loglike ] = em_vonmises(Omega, k, max_iter, multiple_ks, th, R0, c);

%Compute cluster assignments
[~,id] = max(w,[],1);

%Plot clusters
for i=1:k
    set(gca, 'ColorOrderIndex', i)
    for j=find(id == i)
        %scatter(Alpha(1, j), Alpha(end, j));
        plot(X{j}(1,:), X{j}(2,:));
        hold on;
        set(gca, 'ColorOrderIndex', i);
    end
    axis tight;
    pause;
end
