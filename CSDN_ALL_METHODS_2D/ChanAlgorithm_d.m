function X = ChanAlgorithm_d(BSN, MSP, R, Noise, MeaDist)
%CHANALGORITHM_D ����������ʵ�����߶�λ�е�CHAN�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - Noise ������
%               - MeaDist ��������
%               - R    ΪС���뾶����λ(meter)��
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: ChanAlgorithm_d.m


%   ������飺
if  nargout>1,
    error('Too many output arguments.');
end
if nargin ~= 5,
    error('Wrong number of input arguments.');
end


% �㷨��ʼ��
BS = R*NetworkTop(BSN);
MS = R*MSP;

% �������ʣ�
c = 3*10^8; % ���ߵ粨�����ٶ�
Q = 0.5*eye(BSN-1); % TDOA��������Э�������

% ��һ��LS��
% Ri
R1 = sqrt(MS(1)^2 + MS(2)^2);
K1 = 0;
for i = 1: BSN-1,
    R0(i) = sqrt((BS(1,i+1) - MS(1))^2 + (BS(2,i+1) - MS(2))^2);
end

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

% ����������Թ���ֵ����B��
B = eye(BSN-1);
for i = 1: BSN-1,
    B(i,i) = sqrt((BS(1,i+1) - Za0(1))^2 + (BS(2,i+1) - Za0(2))^2);
end

% mFI:
mFI = B*Q*B;

% ��һ��LS�����
Za1 = pinv(Ga'*pinv(mFI)*Ga)*Ga'*pinv(mFI)*h';

if Za1(3) < 0,
    Za1(3) = abs(Za1(3));
end

%***************************************************************

% �ڶ���LS��
% ��һ��LS�����Э���
CovZa = pinv(Ga'*pinv(mFI)*Ga);

% sB��
sB = eye(3);
for i = 1: 3,
    sB(i,i) = Za1(i);
end;

% sFI��
sFI = 4*sB*CovZa*sB;

% sGa��
sGa = [1, 0; 0, 1; 1, 1];

% sh
sh  = [Za1(1)^2; Za1(2)^2; Za1(3)^2];

% �ڶ���LS�����
Za2 = pinv(sGa'*pinv(sFI)*sGa)*sGa'*pinv(sFI)*sh;

% ���:
Za = sqrt(abs(Za2));

out = Za;

if nargout == 1,
    X = out;
elseif nargout == 0,
    disp(out);
end