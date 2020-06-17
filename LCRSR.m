function [ V, S,p,obj ] = LCRSR( X, r, lambda,maxIter )
% X: cell, input data 
% r: VµÄÎ¬¶È £¨VÎªn*r);
% lambda: the trade-off parameter,  lambda = rho*rou, (rou > 1);
% 
if ~exist('maxIter','var')
    maxIter = 50;
end

viewNum = length(X);
nSmp = size(X{1},1);
I = eye(nSmp);

S0 = cell(viewNum,1);
S = cell(viewNum,1);
for v = 1:viewNum
    S0{v} = zeros(size(X{v}))';
end


p = ones(viewNum,1);
obj = zeros(maxIter,1);
thresh = 1e-4;
for iter = 1:maxIter
    %%%% update V
    sumG = zeros(nSmp);
    X_S = cell(viewNum,1);
    for v = 1:viewNum
        X_S{v} = X{v}' - S0{v};
        sumG = sumG + X_S{v}'*X_S{v}/p(v);
    end
    sumG = max(sumG,sumG');
    [V, eigvalue] = eigs(sumG,r);
    
    
    %%% update S(v)
    rou = zeros(viewNum,1);
    IVV = (I - V*V');
%     rou0 = eigs(IVV,1);  %%%% due to V'*V = I, eigs(IVV,1) = 1
    rou0 = 1;
    for v = 1:viewNum
        rou(v) = rou0/p(v);
        grad_S = (S0{v} - X{v}')*IVV/p(v); 
        temp = S0{v} - grad_S/rou(v);
        S{v} = max(abs(temp)- lambda/rou(v),0).*sign(temp);
    end
    
    %%%% update p
    res = zeros(viewNum,1);
    for v = 1:viewNum
        X_S{v} = X{v}' - S{v};
        res(v) = norm(X_S{v}*IVV,'fro');
    end
    p = res;
%     p = res./sum(res);
    
    obj(iter) = 0;
    for v = 1:viewNum
        obj(iter) = obj(iter) + sum(sum(abs(S{v})));
    end
    obj(iter) = lambda*obj(iter)+ sum(res);
    
    if max(res) < thresh && iter >= 10
        break;
    end
    
    S0 = S;
end
obj(iter+1:end) = [];


