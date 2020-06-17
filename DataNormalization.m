
function [XTrain,xmean, xstd,opar]=DataNormalization(XTrain,ixnorm,ipar)
% function to normalize the input and the output data
% input
%       ixnorm      type of nomrlaization
%                   -1 no normalization
%                   0 vector-wise normalizatin by L2 norm: each sample's
%                     L2 norm equal to 1
%                   1 variable-wise normalizatin by mean and standard
%                     deviation: the result training data has zero mean and
%                     std 1, similar with zscore
%                   2 normalization by projection onto a ball
%                   3 vector-wise normalization by L1 norm: each sample's
%                     L1 norm equal to 1
%                   4 vector-wise normalization by L_infty norm
%                   5 variable-wise normalization by L1 norm; median +
%                     MAD,median of absolute deviation
%                   11 centralization by mean of training
%      XTrain       Data matrix which will be normalized. It assumed the
%                   rows are the sample vetors and the columns are variables
%      XTest        Data matrix which will be normalized. It assumed the
%                   rows are the sample vetors and the columns are
%                   variables.
%                   In variable-wise normalization it herites the means
%                   and standard deviation from the XTrain, otherwise it
%                   is normalized independently
%      ipar         it is considered when ixnorm=2, projection onto ball,
%                   where it is the radius of the ball
%  output
%      XTrain       Data matrix which is the result of the normalization
%                   of input XTrain. It assumed the rows are the sample
%                   vetors and the columns are variables
%      XTest        Data matrix which is the result of the normalization
%                   of input XTest. It assumed the rows are the sample
%                   vetors and the columns are variables.
%      opar         the radius in case of ixnorm=2.
%
if ~exist('ixnorm','var')
    ixnorm = 1;
end
if ~exist('ipar','var')
    ipar = 1;
end
opar=0;
[mtrain,n]=size(XTrain);
xmean = 0;
xstd = 0;

% normaliztion
%          ixnorm=-1;
switch ixnorm
    % element-wise normalization by L2 norm
    case 0
        n=size(XTrain,2);
        xsum1=sum(XTrain.^2,2);
        xsum1=sqrt(xsum1);
        xsum1=xsum1+(xsum1==0);
        
        XTrain=XTrain./(xsum1*ones(1,n)+eps);
        
    case 1
        XTrain = zscore(XTrain);
        XTrain(~isfinite(XTrain)) = 0;
    case 11
        % variable-wise normalizatin by mean and standard deviation
        xmean=mean(XTrain,1);
        xstd=std(XTrain,0,1);
        xstd=xstd+(xstd==0);
        XTrain=(XTrain-ones(mtrain,1)*xmean)./(ones(mtrain,1)*xstd);
    case 12
        xmean=mean(XTrain,1);
        XTrain=(XTrain-ones(mtrain,1)*xmean);
        % normalization by projection onto a ball
    case 2 % stereographic projection
        xmean=mean(XTrain,1);
        XTrain=XTrain-ones(mtrain,1)*xmean;
        
        xsum1=sum(XTrain.^2,2);
        %            R=max(sqrt(xsum1));
        %            R=1;
        R=ipar;
        %              R=mean(sqrt(xsum1));
        opar=R;
        xhom=ones(mtrain,1)./(xsum1+R^2);
        xhom2=xsum1-R^2;
        
        XTrain=[2*R^2*XTrain.*(xhom*ones(1,n)),R*xhom2.*xhom];
        
       
        % element-wise normalization by L1 norm
    case 3
        n=size(XTrain,2);
        xsum1=sum(abs(XTrain),2);
        xsum1=xsum1+(xsum1==0);
        
        XTrain=XTrain./(xsum1*ones(1,n));
        
                % element-wise normalization by L_infty norm
    case 4
        n=size(XTrain,2);
        xsum1=max(abs(XTrain),[],2);
        xsum1=xsum1+(xsum1==0);
        
        XTrain=XTrain./(xsum1*ones(1,n));
        
        % variable-wise normalization based on  L1 norm
    case 5
        xmedian=median(XTrain,1);
        
        xmad=median(abs(XTrain-ones(mtrain,1)*xmedian),1);
        
        xmad=xmad+(xmad==0);
        XTrain=(XTrain-ones(mtrain,1)*xmedian)./(ones(mtrain,1)*xmad);
        xmean = xmedian; xstd = xmad;
        
end

