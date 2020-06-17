function DI=dunns(data,ind,distM)   
%%% Dunn's index for clustering compactness and separation measurement
%    dunns(data,ind)
%    data: the data matrix,each row is a point
%    ind   = Indexes for each data point aka cluster to which each data point
%             belongs, label indicator vector or label indicator matrix
%    written by Hong Tao, 2018/6/12.

if isvector(ind)
    ind_name = unique(ind);
    clusters_number = length(ind_name);
    nSmp = length(ind);
    Y = zeros(nSmp,clusters_number);
    for ii = 1:clusters_number
        Y(ind == ind_name(ii),ii) = 1;
    end
elseif ismatrix(ind)
    if size(ind,1) == size(data,1)
        Y = ind;
    else
        Y = ind';
    end
else
    error('unknown label indicator type!');
end
if nargin < 3
    distM=squareform(pdist(data));    %%%% calculate the distance between pairwise data points
end
S = Y*Y';

betwClustDis = min(min(distM(S == 0)));  %%%%% minimum distance between points belonging to different clusters
maxClustDia = max(max(distM(S == 1)));   %%%%% maximum cluster diameter


DI=betwClustDis/maxClustDia;
end