function [ maxpart ] = compute_partitions( X, smooth_p )
%COMPUTE_PARTITIONS Finds the maximum number of characteristic points of a
%trajectory dataset
%   Detailed explanation goes here

maxpart = 0;
for i=1:length(X)
    %Fit smoothing  splines
    df = diff(X{i},1,2);
    t = cumsum([0,sqrt([1 1]*(df.*df))]);
    cs = csaps(t,X{i},smooth_p);
    N = 2*length(X{i});
    xyn = fnval(cs, t(1):(t(end)-t(1))/N:t(end)); %Smoothed
    p = characteristic_points(xyn); %Count characteristic points
    if p > maxpart
        maxpart = p;
    end
end

end

