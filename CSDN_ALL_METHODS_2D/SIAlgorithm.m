function X = SIAlgorithm(BSN, MSP, Radius, Noise)
%SIALGORITHM ����������ʵ�����߶�λ�е�SI�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - Noise ������
%               - R    ΪС���뾶����λ(meter)��
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: SIAlgorithm.m


%   ������飺
if  nargout>1,
    error('Too many output arguments.');
end
if nargin<2 | nargin>4,
    error('Wrong number of input arguments.');
end
if BSN < 3,
    error('The number of BSs must be larger than 3 for this program.');
end
flag = size(MSP);
if flag(1)~=1 | flag(2)~=2,
    error('Wrong position vector!');
end

% ��ʼ������
BS = Radius*NetworkTop(BSN);
MS = Radius*MSP;

% % TDOAЭ�������Q��
% c = 3*10^8; % ���ߵ粨�����ٶ�
% Dev = Noise/(c*c); % TDOA��������

% Ri1
R1 = sqrt(MS(1)^2 + MS(2)^2);
for i = 1: BSN-1,
    R(i) = sqrt((BS(1, i+1) - MS(1))^2 + (BS(2, i+1) - MS(2))^2);
end
for i = 1: BSN-1,
    Ri1(i) = R(i) - R1 + Noise*randn(1);
end
    
% W
W = eye(BSN-1);

% delt
for i = 1: BSN-1,
    K(i) = BS(1, i+1)^2 + BS(2, i+1)^2;
end
for i = 1: BSN-1,
    delt(i) = K(i) - Ri1(i)^2;
end

% Pd orthognol
I = eye(BSN-1);
coef = Ri1*Ri1';
Pd_o = I - (Ri1'*Ri1/coef);
    
% S
for i = 1: BSN-1,
    S(i, 1) = BS(1, i+1);
    S(i, 2) = BS(2, i+1);
end

% �����
    Za = 0.5*inv(S'*Pd_o*W*Pd_o*S)*S'*Pd_o*W*Pd_o*delt';
if nargout == 1,
    X = Za;
elseif nargout == 0,
    disp(Za);
end