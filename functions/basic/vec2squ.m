function M = vec2squ(v)
%VEC2SQU �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[a,b] = size(v);
sa = floor(sqrt(a));
if b == 1 && a == sa^2
    M = reshape(v,sa,sa);
else
    M = 0;
end
end
