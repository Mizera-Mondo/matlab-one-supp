function v = mat2vec(M)
%MAT2VEC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[a,b] = size(M);
v = reshape(M,a*b,1);
end
