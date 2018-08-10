% Plots the distortion function. Can be used with k-means and SSNMF
load % INSERT DATASET HERE
% Parameters for dataset
N = 50;
p = 1;  % Smoothing (max=0, min=1)
invariant2rot = 0;  % 1 for rotation invariance, 0 otherwise
kmin = 1;  % Minimum number of clusters to consider
kmax = 10;  % Maximum number of clusters to consider


%Pre-processing: interpolation + conversion to angles
Omega = convert_dataset(X, N, p, invariant2rot);

residuals = zeros(kmax-kmin+1, 1);
repeat = 20;
max_iter = 100;

for k=kmin:kmax
    for i=1:repeat
%         %K-MEANS
%         [id, centroids, res] = circkmeans(Omega, k, max_iter);
        %NMF
        [~, ~, ~, res] = angular_ssnmf(Omega, k, k, 1, 0.1);
        residuals(k) = residuals(k) + res;
    end
end
residuals = residuals/repeat;

plot(kmin:kmax, residuals)
