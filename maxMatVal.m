function [ maxVal,maxIdx ] = maxMatVal( M )
% ����һ������������ֵ������Ӧ��idx
% maxIdx.rowIdx ��Ӧ��һά,maxIdx.colIdx��Ӧ�ڶ�ά
[val,idx] = max(M);
[maxVal,idxx] = max(val);
maxIdx.colIdx = idxx;
maxIdx.rowIdx = idx(idxx);


end

