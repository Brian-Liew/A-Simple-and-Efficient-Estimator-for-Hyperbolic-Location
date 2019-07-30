function X = ChanAlgorithm_B(BSN, MSP, Radius, Noise)
%CHANALGORITHM ����������ʵ�����߶�λ�е�CHAN�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - Noise ������
%               - R    ΪС���뾶����λ(meter)��
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: ChanAlgorithm.m


%   ������飺
if  nargout>1,
    error('Too many output arguments.');
end
if nargin<2 | nargin>4,
    error('Wrong number of input arguments.');
end


% �㷨��ʼ��
BS = Radius*NetworkTop(BSN);
MS = Radius*MSP;

% �������ʣ�
Q = eye(BSN-1);

% ��һ��LS��
% Ri
K1 = 0;
for i = 1: BSN,
    R0(i) = sqrt((BS(1,i) - MS(1))^2 + (BS(2,i) - MS(2))^2);
end

for i = 1: BSN-1,
    R(i) = R0(i+1) - R0(1) + Noise*randn(1);
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
Za0 = inv(Ga'*inv(Q)*Ga)*Ga'*inv(Q)*h';

% ����������Թ���ֵ����B��
B = eye(BSN-1);
for i = 1: BSN-2,
    B(i,i) = sqrt((BS(1,i+1) - Za0(1))^2 + (BS(2,i+1) - Za0(2))^2);
end
B(BSN-1,BSN-1) = 0.1;

% FI:
FI = B*Q*B;

% ��һ��LS�����
Za1 = inv(Ga'*inv(FI)*Ga)*Ga'*inv(FI)*h';

if Za1(3) < 0,
    Za1(3) = abs(Za1(3));
%     Za1(3) = 0;
end
%***************************************************************

% �ڶ���LS��
% ��һ��LS�����Э���
CovZa = inv(Ga'*inv(FI)*Ga);

% sB��
sB = eye(3);
for i = 1: 3,
    sB(i,i) = Za1(i);
end

% sFI��
sFI = 4*sB*CovZa*sB;

% sGa��
sGa = [1, 0; 0, 1; 1, 1];

% sh
sh  = [Za1(1)^2; Za1(2)^2; Za1(3)^2];

% �ڶ���LS�����
Za2 = inv(sGa'*inv(sFI)*sGa)*sGa'*inv(sFI)*sh;

Za = sqrt(abs(Za2));

% Za = sqrt(Za2);
% 
% % ���:
% if Za1(1) < 0,
%     out1 = -Za(1);
% else
%     out1 = Za(1);
% end
% if Za2(1) < 0,
%     out2 = -Za(2);
% else
%     out2 = Za(2);
% end
% 
% out = [out1;out2];

out = Za;

if nargout == 1,
    X = out;
elseif nargout == 0,
    disp(out);
end