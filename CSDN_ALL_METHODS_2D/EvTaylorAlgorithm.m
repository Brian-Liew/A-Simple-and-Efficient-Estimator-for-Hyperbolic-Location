function X = EvTaylorAlgorithm(BSN, MSP, Radius, Noise)
%
%TAYLORALGORITHM_D ����������ʵ�����߶�λ�е�TAYLOR�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - Noise ������
%               - R    ΪС���뾶����λ(meter)��
%               - MeaDist �������롣
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: TaylorAlgorithm_d.m


% ������飺
if  nargout ~= 1 & nargout ~= 0,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end

% ��ʼ������
BS = Radius*NetworkTop(BSN);
MSi = Radius*MSP;

% MeaDist
for i = 1: BSN,
    MeaDist(i) = sqrt((MSi(1) - BS(1,i))^2 + (MSi(2) - BS(2,i))^2) + MeaNoise(Noise);
end

MS = LSAlgorithm_d(BSN, MSP, Radius, MeaDist);

% TDOAЭ�������Q��
c = 3*10^8; % ���ߵ粨�����ٶ�
Dev = Noise/(c*c); % TDOA��������
Q = 0.5*eye(BSN-1); % TDOA��������Э�������
    
% ��ʼ����λ�ã�
iEP = MS;

% h0:
for i = 1: BSN-1,
    h0(i) = MeaDist(i+1) - MeaDist(1);
end

alpha = 100;

% �㷨��ʼ��
%while alpha > 10,
for n = 1: 50,
    % Rn:
    R1 = sqrt(iEP(1)*iEP(1) + iEP(2)*iEP(2));
    for i =1: BSN-1,
        R(i) = sqrt((iEP(1) - BS(1,i+1))*(iEP(1) - BS(1,i+1)) + (iEP(2) - BS(2,i+1))*(iEP(2) - BS(2,i+1)));        
    end
    
    % ht:
    for i = 1: BSN-1,
        h(i) = h0(i) - (R(i) - R1);
    end
    ht = h';
    
    % Gt:
    for i = 1: BSN-1,
        Gt(i, 1) = -iEP(1)/R1 - (BS(1, i+1) - iEP(2))/R(i);
        Gt(i, 2) = -iEP(2)/R1 - (BS(2, i+1) - iEP(2))/R(i);
    end
    
    % delt:
    delt = pinv(Gt'*pinv(Q)*Gt)*Gt'*pinv(Q)*ht;
    
    EP = iEP + delt';

    iEP = EP;
end

% ������:
if nargout == 1,
    X = abs(EP);
elseif nargout == 0,
    disp(abs(EP));
end