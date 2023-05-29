function [vecSorted, valSorted] = sortEigen(vec, val, direction)
%SORTEIGEN 此处显示有关此函数的摘要
%   此处显示详细说明
[~, index] = sort(diag(val), direction);
valSorted = val(index, index);
vecSorted = vec(:, index);
end
