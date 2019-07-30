function X = FangAlgorithm(MSP, R, Noise)
%  本程序实现无线定位中的FANG算法
%  FANGALGORITHM
%    参数说明：
%       MSP：移动台相对位置；
%       R：  小区半径；
%       Noise: 测距误差方差。
%       X：  输出移动台的估计位置。
%  Also see: FangAlgorithm.


%  输入参数检测:
    if  nargout>1,
        error('Too many output arguments.');
    end
    if nargin~=3,
        error('Wrong number of input arguments.');
    end
%  初始参数：
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
    
%  输出：
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