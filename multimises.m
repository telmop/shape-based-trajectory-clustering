function [ p ] = multimises( theta, mu, kappa )
%VONMISES_PDF Multivariate Von Mises (with \Lambda = 0^{dxd})
%   Detailed explanation goes here

p = ones(1, size(theta, 2));

for m=1:size(theta, 2)
    for i=1:length(mu)
        p(m) = p(m)*exp(kappa(i)*(cos(theta(i, m)-mu(i))-1))/(2*pi*besseli(0, kappa(i), 1));
    end
end

end
