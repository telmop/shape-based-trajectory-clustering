% Plots the consistency function
load % INSERT DATASET HERE
% Parameters for dataset
N = 50;
p = 1;  % Smoothing (max=0, min=1)
invariant2rot = 0;  % 1 for rotation invariance, 0 otherwise
kmin = 1;  % Minimum number of clusters to consider
kmax = 10;  % Maximum number of clusters to consider

%Pre-processing: interpolation + conversion to angles
Omega = convert_dataset(X, N, p, invariant2rot);
repeat = 20;

%NMF
[bestk, bestW, bestH, consistency] = angular_ssnmf(Omega, kmin, kmax, repeat, 0.1);

%Compute cluster assignments
[~, id] = max(bestH, [], 1);

%Plot clusters
plot(2:10, consistency(2:end));

% for i=1:bestk
%     set(gca, 'ColorOrderIndex', i)
%     for j=find(id == i)
%         plot(X{j}(1,:), X{j}(2,:));
%         hold on;
%         set(gca, 'ColorOrderIndex', i);
%     end
%     axis tight;
%     pause;
% end
