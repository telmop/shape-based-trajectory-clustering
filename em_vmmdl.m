function [ bestk, bestalpha, bestmu, bestkappa, mdlcosts ] = em_vmmdl( data, kmin, kmax, max_iter, multiple_ks, th, R0, c, repeat )
%EM_VMMDL Applies MDL to the mixture of von Mises algorithm to determine
%the number of clusters
%   Detailed explanation goes here

if nargin < 7
    repeat = 1;
end

if nargin < 6
    th = 1e-4;
end

min_mdl = inf;
[d, m] = size(data);

if multiple_ks
    N = 2*d;
else
    N = d+1;
end

mdlcosts = zeros(length(kmin:kmax), 1);
for k=kmax:-1:kmin
    fprintf('Number of clusters: %d\n', k);
    best_loglike = -inf;
    for i=1:repeat
        [~, alpha, mu, kappa, loglike] = em_vonmises(data, k, max_iter, multiple_ks, th, R0, c);
        if loglike > best_loglike
            best_loglike = loglike;
            mdlcosts(k) = -loglike + (N+1)*0.5*k*log(m);
        end
        mdl_cost = -loglike + (N+1)*0.5*k*log(m);
        if mdl_cost <= min_mdl && ~isnan(mdl_cost)
            min_mdl = mdl_cost;
            bestalpha = alpha;
            bestmu = mu;
            bestkappa = kappa;
            bestk = k;
        end
    end
    
    fprintf('Best MDL_cost: %f\n', min_mdl);
%     if -best_loglike + (N+1)*0.5*log(m) > min_mdl %No longer worth it to keep searching
%         disp('Break condition: loglikelihood is too small (no better solutions for smaller ks).');
%         break;
%     end
end

end

function [ id1, id2 ] = closest_mus(mu)
%CLOSEST_MUS Determines the two closest centroids
%   Detailed explanation goes here

min_dist = inf;

for j=1:size(mu, 2)
    dist = sum(1-cos(bsxfun(@minus, mu, mu(:, j))), 1);
    dist(j) = max(dist)+1; %so we don't select the same mu
    [mind, id] = min(dist);
    
    if mind < min_dist
        min_dist = mind;
        id1 = j;
        id2 = id;
    end
end

end

