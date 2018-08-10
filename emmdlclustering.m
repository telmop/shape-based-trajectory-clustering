% Plots the best clusters found using VMMs
load % INSERT DATASET HERE
% Parameters for dataset
N = 50;
p = 1;  % Smoothing (max=0, min=1)
invariant2rot = 0;  % 1 for rotation invariance, 0 otherwise
kmin = 1;  % Minimum number of clusters to consider
kmax = 10;  % Maximum number of clusters to consider


%Pre-processing: interpolation + conversion to angles
Omega = convert_dataset(X, N, p, invariant2rot);

%FMM + MDL
max_iter = 100;
multiple_ks = 0;
th = 1e-4;
repeat = 20;
R0 = -5e-5;
c = 5e-5;
[ bestk, bestalpha, bestmu, bestkappa, mdl_costs ] = em_vmmdl(Omega, kmin, kmax, max_iter, multiple_ks, th, R0, c, repeat);

%Plot results
disp('Plotting results...');

%Compute cluster assignments
w = [];
for j=1:bestk
    w(j, :) = bestalpha(j)*multimises(Omega, bestmu(:, j), bestkappa(:, j));
end
w = bsxfun(@rdivide, w, sum(w, 1));
[~,id] = max(w,[],1);

%Plot clusters
for i=1:bestk
    %set(gca, 'ColorOrderIndex', i)
    for j=find(id == i)
        set(gca, 'ColorOrderIndex', i);
        plot(X{j}(1,:), X{j}(2,:));
        hold on;
    end
    axis tight;
    %pause;
end
