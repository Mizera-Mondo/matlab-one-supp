function M = vec2mat(v, r)
%VEC2MAT 此处显示有关此函数的摘要
%   此处显示详细说明
[a,b] = size(v);
s = a/r;
if b == 1 && s == round(s)
    M = reshape(v,r,s);
else
    M = 0;
end
end
