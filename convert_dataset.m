function [ Omega ] = convert_dataset( X, N, p, invariant2rot )
%CONVERT_DATASET Converts a 2D trajectory dataset to a dataset of tangent
%angles
%   Detailed explanation goes here

if nargin < 4
    invariant2rot = 0;
end
if invariant2rot
    Omega = zeros(N-1, length(X));
else
    Omega = zeros(N, length(X));
end

for i=1:length(X);
    Omega(:, i) = xy2theta(X{i}, N, p, invariant2rot)';
end

end

function [ theta ] = xy2theta( xy, N, lambda, inv_rot )
%XY2THETA Converts a 2D trajectory to a uniformly spaced sequence of the
%tangent's angles.
%   Detailed explanation goes here
if nargin < 4
    inv_rot = 0;
end

df = diff(xy,1,2);
t = cumsum([0,sqrt([1 1]*(df.*df))]); %Estimates time axis
[cs1, ~] = csaps(t, xy, lambda); %Smothing spline
deriv = fnval(fnder(cs1), t(1):(t(end)-t(1))/(N-1):t(end)); %Estimates 1st derivative at equally space points
theta = atan2(deriv(2, :), deriv(1, :)); %Converts gradients to angles

if inv_rot
    theta = diff(theta);
end
end

