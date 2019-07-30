function X = EvChanAlgorithm(BSN, MSP, R, Noise)
%EvChanAlgorithm ʵ�ָĽ�Chan�㷨
%EVCHANALGORITHM ����������ʵ�����߶�λ�е�CHAN�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - Noise ������
%               - R    ΪС���뾶����λ(meter)��
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: EvChanAlgorithm.m


%   ������飺
if  nargout>1,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end

% �㷨��ʼ
MS = R*MSP;
BS = R*NetworkTop(BSN);

% MeaDist
for i = 1: BSN,
    MeaDist(i) = sqrt((MS(1) - BS(1,i))^2 + (MS(2) - BS(2,i))^2) + MeaNoise(Noise);
end

% Chan��Taylor����λ��
EMSC = ChanAlgorithm_d(BSN, MSP, R, Noise,MeaDist);
EMSTC = TaylorAlgorithm_d(BSN, EMSC'/R, R, Noise, MeaDist);
EMSTC = EMSTC';

% �в�
ResC = Residual(MeaDist,EMSC,BSN,R);
ResTC = Residual(MeaDist,EMSTC,BSN,R);

% ����λ��
EMS = (EMSC/ResC + EMSTC/ResTC)/(1/ResC + 1/ResTC);

% ������
if nargout == 1,
    X = EMS;
elseif nargout == 0,
    disp(EMS);
end