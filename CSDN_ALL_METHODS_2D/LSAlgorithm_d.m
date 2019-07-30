function X = LSAlgorithm_d(BSN, MSP, R, MeaDist)
%CHANALGORITHM_D ����������ʵ�����߶�λ�е�CHAN�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - MeaDist ��������
%               - R    ΪС���뾶����λ(meter)��
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: ChanAlgorithm_d.m


%   ������飺
if  nargout>1,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end


% �㷨��ʼ��
BS = R*NetworkTop(BSN);
MS = R*MSP;

% �������ʣ�
Q = 0.5*eye(BSN-1); % TDOA��������Э�������

% LS��
% Ri,Ki
K1 = 0;
for i = 1: BSN-1,
    R(i) = MeaDist(i+1) - MeaDist(1);
    K(i) = BS(1,i+1)^2 + BS(2,i+1)^2;
end

% Ga
for i = 1: BSN-1,
    Ga(i,1) = -BS(1, i+1);
    Ga(i,2) = -BS(2, i+1);
    Ga(i,3) = -R(i);
end

% h
for i = 1: BSN-1,
    h(i) = 0.5*(R(i)^2 - K(i) + K1);
end

% �ɣ�14b������B�Ĺ���ֵ��
Za0 = pinv(Ga'*pinv(Q)*Ga)*Ga'*pinv(Q)*h';

% ���:
out = [Za0(1),Za0(2)];

if nargout == 1,
    X = out;
elseif nargout == 0,
    disp(out);
end