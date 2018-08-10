function [ partp ] = characteristic_points( track )
%CHARACTERISTIC_POINTS Computes the number of points necessary to describe
% the trajectory using MDL
n = size(track, 2);
start = 1;
len = 1;
partp = 1; %Start point only
%tolerance = 0;
while start + len <= n
    current = start + len;
    cost_par = mdl(track, start, current, 1);
    cost_nopar = mdl(track, start, current, 0);
    if cost_par > cost_nopar
        %if tolerance == 2
        %    current = current - tolerance;
        %    tolerance = 0;
            start = current - 1;
            len = 1;
            partp = partp + 1;
        %else
        %    len = len + 1;
        %    tolerance = tolerance + 1;
        %end
    else
        len = len + 1;
        %tolerance = 0;
    end
end
partp = partp + 1; %Add end point
end

function [ mdl_score ] = mdl( track, start, current, part )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if part
    partition = [start, current];
else
    partition = start:current;
end

L_H = 0;
L_DH = 0;

for j=1:length(partition)-1
    pc1 = track(:, partition(j));
    pc2 = track(:, partition(j+1));
    
    L_H = L_H + ceil(log2(norm(pc2-pc1)));
    
    for k=partition(j):partition(j+1)-1
        p1 = track(:, k);
        p2 = track(:, k+1);
        dp = dperp(pc1, pc2, p1, p2);
        dt = dtheta(pc1, pc2, p1, p2);
        %fprintf('k: %d, prt(j): %d, prt(j+1): %d\n', k, partition(j), partition(j+1));
        
        %To take logarithms
        if dp < 1
            dp = 1;
        end
        if dt < 1
            dt = 1;
        end
        
        L_DH = L_DH + ceil(log2(dp)) + ceil(log2(dt));
        %L_DH_perp = L_DH_perp + dp;
        %L_DH_theta = L_DH_theta + dt;
    end
    %L_DH = log2(L_DH_perp) + log2(L_DH_theta);
end

mdl_score = L_H + L_DH;

end

function [ dp ] = dperp( pc1, pc2, p1, p2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if norm(pc2-pc1) >= norm(p2-p1)
    s = pc1;
    e = pc2;
    p = [p1, p2];
else
    s = p1;
    e = p2;
    p = [pc1, pc2];
end

l1 = distance2segment(s, e, p(:, 1));
l2 = distance2segment(s, e, p(:, 2));
if l1 == 0 && l2 == 0
    dp = 0;
else
    dp = (l1^2+l2^2)/(l1+l2);
end
end

function [ d ] = distance2segment( s, e, p )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if norm(e-s) == 0 || norm(p-s) == 0
    d = 0;
    return;
end

proj_p = ((p-s)'*(e-s)/(norm(e-s)^2))*(e-s);
d = norm(p-s-proj_p);
end

function [ dt ] = dtheta( pc1, pc2, p1, p2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Lc = norm(pc2-pc1);
Lp = norm(p2-p1);

vc = (pc2-pc1)/norm(pc2-pc1);
vp = (p2-p1)/norm(p2-p1);
ctheta = vc'*vp;

if ctheta > 1 %Numerical problems, which lead to a cos(theta) slightly higher than 1 (~1+1e-7)
    ctheta = 1;
end

if ctheta < 0 %Angle > 90º, d_theta = min(Lc, Lp);
    ctheta = 0;
end

dt = min(Lc, Lp)*sin(acos(ctheta));
% theta = atan2(norm(cross([pc2-pc1;0],[p2-p1;0])),dot(pc2-pc1,p2-p1))
% if -pi/2 < theta < pi/2
%     dt = min(Lc, Lp)*sin(theta);
% else
%     dt = min(Lc, Lp);
% end
end
