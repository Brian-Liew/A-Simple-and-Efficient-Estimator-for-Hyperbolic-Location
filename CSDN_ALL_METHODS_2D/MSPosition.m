function X = MSPosition()
%  ��������1/12С�����������MS��λ�ã�����С���뾶Ϊ1
%  MSPosition
%    ����˵����
%       �޲�����
%  Also see: MSPosition.


%  �������:
 if  nargout>1,
        error('Too many output arguments.');
 end

% ��������ƶ�̨λ�ã�
x = sqrt(3)*rand(1)/2;
y = sqrt(3)*x*rand(1)/3;

% ��������
if nargout == 1,
    X = [x, y];
else
    disp([x, y]);
end

