function M = vec2mat(v, r)
%VEC2MAT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[a,b] = size(v);
s = a/r;
if b == 1 && s == round(s)
    M = reshape(v,r,s);
else
    M = 0;
end
end
