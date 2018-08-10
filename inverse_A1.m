function [ y ] = inverse_A1(x, R0, c)
%INVERSE_A1 Approximation of the inverse of A1(x) = I1(x)/I0(x)
%   Detailed explanation goes here

y = zeros(length(x), 1);
x = (x+R0)/(1+c); %Prior on kappa, avoids numerical problems without adulterating likelihood
for i=1:length(x)
    if x(i) < 0.53
        y(i) = 2*x(i) + x(i)^3 + 5*x(i)^5/6;
    elseif x(i) >= 0.53 && x(i) < 0.85
        y(i) = -0.4 + 1.39*x(i) + 0.43/(1-x(i));
    else
        y(i) = 1/(x(i)^3 - 4*x(i)^2 + 3*x(i));
    end
end
end
