function [ bestk, bestW, bestH, consst ] = angular_ssnmf( data, kmin, kmax, repeat, beta )
%ANGULAR_SSNMF Applies SSNMF to clustering angular data
%    bestk, bestW, and bestH are the parameters of the best solution
%    consst is the consistency of the best solution
if nargin < 4
    repeat = 1;
end

[d, m] = size(data);
tmp = exp(1j*data);

D = [real(tmp); imag(tmp)];

max_consistency = 0;
consst = [];
option.beta = beta;
for k=kmin:kmax
    fprintf('Number of clusters: %d\n', k);
    C = zeros(m, m);
    for i=1:repeat
        %fprintf('Iteration %d of %d\n', i, repeat);
        [W, H, ~, ~, res] = sparseseminmfnnls(D, k, option);
        C = C + connectivity_matrix(H);
    end
%     
    C = C/repeat; %Mean connectivity matrix
    
    consistency = sum(sum(4*(C-1/2).^2))/m^2; %Computes consistency
    consst(k) = consistency;
    
    if consistency > max_consistency
        bestk = k;
        bestW = angle(W(1:d, :) + 1j*W(d+1:end, :)); %Computes the centroid matrix (in radians)
        bestH = H;
        max_consistency = consistency;
    end
%     fprintf('Best consistency: %f\n', max_consistency); 
end
bestk = kmin;
bestW = W;
bestH = H;
end

function [ C ] = connectivity_matrix( H )
%CONNECTIVITY_MATRIX Computes the connectivity matrix associated with the
%results of SSNMF
%   Detailed explanation goes here

m = size(H, 2);
[~, id] = max(H, [], 1); %Computes cluster assignments

C = zeros(m, m); %Connectivity matrix

for i=1:m
    for j=i:m
        C(i, j) = (id(i) == id(j)); %1 if in the same cluster, 0 otherwise
        C(j, i) = C(i, j); %Matrix is symmetric
    end
end

end

function [ ss ] = sure_score( D, W, H )
%SURE_SCORE Computes the SURE score associated with factorization D \approx
%WH
%   Detailed explanation goes here

[r, m] = size(H);
n = size(W, 1);

E = D-W*H;
sqe = norm(E, 'fro')^2;
sigma2 = sqe/m^2; %Variance estimate

lambda = eig(E'*E);
rho = eig(H*H'*W'*W);

term = 0;
for i=1:r
    for j=1:(m-r)
        term = term + 2*lambda(j)/(lambda(j)-rho(i));
    end
end

ss = sqe + 2*sigma2*((m+n)*r - term);

end
