function [ SI ] = SilhouetteIndex( X,gnd,metric )
% calcluate the Silhouette index of a clustering
% default metric is 'Euclidean',
if ~exist('metric','var')
    metric = 'Euclidean';
end
s = silhouette(X,gnd,metric);
cL = unique(gnd);
nc = length(cL);
SI = 0;
for i = 1:nc
    idx = find(gnd == cL(i));
    SI = SI + mean(s(idx));
end
SI = SI/nc;


end

