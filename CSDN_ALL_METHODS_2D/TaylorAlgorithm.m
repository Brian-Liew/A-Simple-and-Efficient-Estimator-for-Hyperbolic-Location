function X = TaylorAlgorithm(BSN, MSP, Radius, Noise)
%
%TAYLORALGORITHM ����������ʵ�����߶�λ�е�TAYLOR�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - Noise ������
%               - R    ΪС���뾶����λ(meter)��
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: TaylorAlgorithm.m


% ������飺
if  nargout ~= 1 & nargout ~= 0,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end

% ��ʼ������
BS = Radius*NetworkTop(BSN);
MS = Radius*MSP;

% TDOAЭ�������Q��
Q = eye(BSN-1);
    
% ��ʼ����λ�ã�
iEP = MS;

% h0:
for i = 1: BSN,
    MeaDist(i) = sqrt((MS(1) - BS(1,i))^2 + (MS(2) - BS(2,i))^2);
end
for i = 1: BSN-1,
    h0(i) = MeaDist(i+1) - MeaDist(1) + Noise*randn(1);   %TDOA����ֵ
end

% �㷨��ʼ��
for n = 1: 10,
    % Rn:
    R1 = sqrt(iEP(1)^2 + iEP(2)^2);
    for i =1: BSN-1,
        R(i) = sqrt((iEP(1) - BS(1,i+1))^2 + (iEP(2) - BS(2,i+1))^2);        
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
    delt = inv(Gt'*inv(Q)*Gt)*Gt'*inv(Q)*ht;
    
    EP = iEP + delt';

    iEP = EP;
end

% ������:
if nargout == 1,
    X = EP;
elseif nargout == 0,
    disp(EP);
end