function X = FangAlgorithm(MSP, R, Noise)
%  ������ʵ�����߶�λ�е�FANG�㷨
%  FANGALGORITHM
%    ����˵����
%       MSP���ƶ�̨���λ�ã�
%       R��  С���뾶��
%       Noise: ������
%       X��  ����ƶ�̨�Ĺ���λ�á�
%  Also see: FangAlgorithm.


%  ����������:
    if  nargout>1,
        error('Too many output arguments.');
    end
    if nargin~=3,
        error('Wrong number of input arguments.');
    end
%  ��ʼ������
    MS = R*MSP;
    BS = R*NetworkTop(3);
    %
    R1 = sqrt(MS(1)^2 + MS(2)^2);
    R2 = sqrt((BS(1,2) - MS(1))^2 + (BS(2,2) - MS(2))^2);
    R3 = sqrt((BS(1,3) - MS(1))^2 + (BS(2,3) - MS(2))^2);
    %
    R21 = R2 - R1 + MeaNoise(Noise);
    R31 = R3 - R1 + MeaNoise(Noise);
    %
    g = ((R31*BS(1,2))/R21 - BS(1,3))/BS(2,3);
    h = (BS(1,3)^2 + BS(2,3)^2 - R31^2 + R31*R21*(1 - (BS(1,2)/R21)^2))/(2*BS(2,3));
    d = -((1 - (BS(1,2)/R21)*(BS(1,2)/R21)) + g^2);
    e = BS(1,2)*(1 - (BS(1,2)/R21)^2) - 2*g*h;
    f = (R21^2/4)*(1-(BS(1,2)/R21)^2)^2 - h^2;
    
%  �����
    root1 = (-e - sqrt(e^2 - 4*d*f))/(2*d);
    root2 = (-e + sqrt(e^2 - 4*d*f))/(2*d);
    
    if root1 > 0,
        EMSX = root1;
    else
        EMSX = root2;
    end
    
    EMSY = g*EMSX + h;
    EMS = [EMSX, EMSY];
    if nargout == 1,
        X = EMS;
    else
        disp(EMS);
    end