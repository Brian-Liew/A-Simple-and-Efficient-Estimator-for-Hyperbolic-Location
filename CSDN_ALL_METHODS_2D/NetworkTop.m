function X = NetworkTop(BSN)
%  �����������������ˣ���������뾶Ϊ1
%  NetworkTop
%    ����˵����
%        BSN:   ���������л�վ��Ŀ 3<=BSN<=7
%  Also see: NetworkTop.


%  �������:
if nargout>1,
    error('Too many output arguments!');
end 
if nargin ~= 1
    error('input arguments error!');
end
if BSN > 7 | BSN < 3,
    error('Overflow!');
end

% 7С���������ˣ�
BS = [0, sqrt(3), 0.5*sqrt(3), -0.5*sqrt(3), -sqrt(3), -0.5*sqrt(3), 0.5*sqrt(3);
      0,        0,         1.5,          1.5,        0,         -1.5,        -1.5];

% BSN��С���������ˣ�
for i = 1 : BSN,
    X(1,i) = BS(1,i);
    X(2,i) = BS(2,i);
end

% ��������
if nargout == 1,
    X;
else
    disp(X);
end
