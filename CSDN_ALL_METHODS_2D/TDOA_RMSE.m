function X = TDOA_RMSE(MS, EMS)
%
%TDOA_RMSE ����������ʵ�����߶�λ����RMSE�ļ���
%               - MS  Ϊ�ƶ�̨����ʵλ�ã�
%               - EMS Ϊ�ƶ�̨�Ĺ���λ�á�
%See also: TDOA_RMSE.m


% ������飺
if  nargout ~= 1 & nargout ~= 0,
    error('Too many output arguments.');
end
if nargin ~= 2,
    error('Wrong number of input arguments.');
end

% �㷨��ʼ��
[n,m] = size(MS);

% �����������ܺͣ�
sum = 0;
for i = 1:n,
    sum = (MS(i,1) - EMS(i,1))^2 + (MS(i,2) - EMS(i,2))^2 + sum;
end

% RMSE:
RMSE = sqrt(sum/n);

% ��������
if nargout == 1,
    X = RMSE;
else
    disp(RMSE);
end
    