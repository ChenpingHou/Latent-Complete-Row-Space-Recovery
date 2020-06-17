%%%%%%%% run experiments of clustering
clear;
addpath(genpath('clustering measure'))
%% results on UCIdermatology
load('data\UCIdermatology');
viewNum = length(M);
nClass = length(unique(gnd));
nSmp = length(gnd);
%%%% data preprocessing
ixnorm = 0;   % parameter for choosing normalization ways
for v = 1:viewNum
    M{v} =DataNormalization(M{v},ixnorm);
end

r = 9; lambda = 1e-6;
tic;
[ V, S,p, obj_ ] = LCRSR( M,r, lambda);
LCRSR_time = toc;

reptimes = 50;
gnd_pre = zeros(nSmp,reptimes);
results = zeros(reptimes,8);
for rr = 1:reptimes
    gnd_pre(:,rr) = kmeans(V,nClass ,'Replicates',10,'EmptyAction','singleton');
    results(rr,:) = ClusteringMeasure(gnd, gnd_pre(:,rr));
end
meanRes(1,:) = mean(results,1);
meanRes(2,:) = std(results,1);

%% results on Wiki
load('data\Wiki');
viewNum = length(M);
nClass = length(unique(gnd));
nSmp = length(gnd);
%%%% data preprocessing
ixnorm = 0;   % parameter for choosing normalization ways
for v = 1:viewNum
    M{v} =DataNormalization(M{v},ixnorm);
end

r = 5; lambda = 1e3;
tic;
[ V, S,p, obj_ ] = LCRSR( M,r, lambda);
LCRSR_time = toc;

reptimes = 50;
gnd_pre = zeros(nSmp,reptimes);
results = zeros(reptimes,8);
for rr = 1:reptimes
    gnd_pre(:,rr) = kmeans(V,nClass ,'Replicates',10,'EmptyAction','singleton');
    results(rr,:) = ClusteringMeasure(gnd, gnd_pre(:,rr));
end
meanRes(1,:) = mean(results,1);
meanRes(2,:) = std(results,1);