%%%%% determining the low-rank parameter using internal clustering
%%%%% measurement

clear;
addpath(genpath('clustering measure'))
load('data\UCIdermatology');
viewNum = length(M);
nClass = length(unique(gnd));
nSmp = length(gnd);
%%%% data preprocessing
ixnorm = 0;   % parameter for choosing normalization ways
for v = 1:viewNum
    M{v} =DataNormalization(M{v},ixnorm);
end
% 
lambda_set = 10.^[-6:6];
r_set = 4:10;
reptimes = 50;

SI_mat = zeros(length(r_set),length(lambda_set));
for rrk = 1:length(r_set)
    r = r_set(rrk);
    for ll = 1:length(lambda_set)
        lambda = lambda_set(ll);
        [ V, S,p, obj_ ] = LCRSR( M, r, lambda);
        restemp = zeros(reptimes,viewNum);
        gnd_pre_temp = zeros(nSmp,reptimes);
        for rr = 1:reptimes
            gnd_pre_temp(:,rr) = kmeans(V,nClass,'Replicates',10,'EmptyAction','singleton');
            for v = 1:viewNum
                restemp(rr,v) = SilhouetteIndex( M{v},gnd_pre_temp(:,rr) );
            end
        end
        mean_restemp = mean(restemp,2);
        SI_mat(rrk,ll) = mean(mean_restemp);
    end
end
[maxSI, maxSIIdx] = maxMatVal( SI_mat );
opt_r = r_set(maxSIIdx.rowIdx);
opt_lambda = lambda_set(maxSIIdx.colIdx);

[ V_SI, S_SI] = LCRSR( M, opt_r , opt_lambda);
gnd_pre_SI = zeros(nSmp,reptimes);
results_SI = zeros(reptimes,8);
for rr = 1:reptimes
     gnd_pre_SI(:,rr) = kmeans(V_SI,nClass,'Replicates',10,'EmptyAction','singleton');
     results_SI(rr,:) =  ClusteringMeasure(gnd, gnd_pre_SI(:,rr));
end
finalRes_SI(1,:) = mean(results_SI,1);
finalRes_SI(2,:) = std(results_SI);





