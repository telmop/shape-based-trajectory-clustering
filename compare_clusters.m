function [ equal ] = compare_clusters( computed, true )
%COMPARE_CLUSTERS Compares the clusters found with the true clusters (only for artificial data)
%   Detailed explanation goes here

kmax = max(true);
adapted = computed;
for k=1:kmax
    ind = find(true == k);
    adapted(computed == computed(ind(1))) = k;
end

equal = sum(adapted==true) == length(true);

end

