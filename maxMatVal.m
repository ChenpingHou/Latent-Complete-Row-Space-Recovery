function [ maxVal,maxIdx ] = maxMatVal( M )
% 计算一个矩阵里的最大值，及相应的idx
% maxIdx.rowIdx 对应第一维,maxIdx.colIdx对应第二维
[val,idx] = max(M);
[maxVal,idxx] = max(val);
maxIdx.colIdx = idxx;
maxIdx.rowIdx = idx(idxx);


end

