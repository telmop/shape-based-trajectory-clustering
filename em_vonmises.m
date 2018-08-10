function [ w, alpha, mu, kappa, loglike ] = em_vonmises( data, k, max_iter, multiple_ks, th, R0, c, init_mu )
%EM_VONMISES Applies the EM algorithm to estimate a Von Mises mixture
%model
%   Detailed explanation goes here
if nargin < 8
    init_mu = [];
end
if nargin < 7
    R0 = 0;
    c = 0;
end

if nargin < 5
    th = 1e-4;
end

if nargin < 4
    multiple_ks = 1;
end

old_loglike = 0;
[alpha, mu, kappa ] = initialize_params(data, k, init_mu);

for a=1:max_iter
    %E-step
    w = compute_w(data, k, alpha, mu, kappa, 0); %Unnormalized
    
    %Compute likelihood
    loglike = sum(log(realmin+sum(w, 1)));
    
    %Termination condition
    if abs(loglike - old_loglike) <= th*abs(old_loglike) || isnan(loglike)
        break;
    end
    old_loglike = loglike;
    
    %Normalize w
    w = bsxfun(@rdivide, w, sum(w, 1));
    
    %M-step
    mu = update_mu(w, data);
    kappa = update_kappa(w, mu, data, multiple_ks, R0, c);
    alpha = update_alpha(w);
end

%fprintf('#iterations: %d, final loglikelihood: %f\n', a, loglike);

end

function [ alpha, mu, kappa ] = initialize_params( data, k, init_mu )
%INITIALIZE_PARAMS Computes the first estimates of the model's parameters
%   Detailed explanation goes here

d = size(data, 1);

alpha = (1/k)*ones(k, 1);

global_mu = angle(sum(exp(1j*data), 2));
glob_kappa = inverse_A1(mean2(cos(bsxfun(@minus, data, global_mu))), 0, 0);

kappa = ones(d, k);

if isempty(init_mu)
    [~, mu] = circkmeans(data, k, 100);
else
    mu = init_mu;
end

end

function [ w ] = compute_w( data, k, alpha, mu, kappa, normalize )
%COMPUTE_W E-Step: computes the w_j^i's
%   Detailed explanation goes here

w = zeros(k, size(data, 2));

for j=1:k
    w(j, :) = alpha(j)*multimises(data, mu(:, j), kappa(:, j));
end

if normalize
    w = bsxfun(@rdivide, w, sum(w, 1)); %Normalize w
end

end

function [ alpha ] = update_alpha( w )
%UPDATE_ALPHA M-Step: updates mixture probabilities
%   Detailed explanation goes here
alpha = mean(w, 2);
end

function [ mu ] = update_mu( w, data )
%UPDATE_MU M-Step: updates means
%   Detailed explanation goes here
mu = angle(exp(1j*data)*w');

end

function [ kappa ] = update_kappa( w, mu, data, multiple_ks, R0, c )
%UPDATE_KAPPA M-Step: updates the dispersion factors
%   Detailed explanation goes here

d = size(data, 1);

k = size(w, 1);
kappa = zeros(d, k);
sw = sum(w, 2);
for j=1:k
    A_k = cos(bsxfun(@minus, data, mu(:, j)))*w(j, :)'/sw(j);
    
    if multiple_ks
        kappa(:, j) = inverse_A1(A_k, R0, c);
    else
        kappa(:, j) = inverse_A1(mean(A_k), R0, c)*ones(d, 1);
    end
end

end