function X = EvChanAlgorithm(BSN, MSP, R, Noise)
%EvChanAlgorithm 实现改进Chan算法
%EVCHANALGORITHM 本函数用于实现无线定位中的CHAN算法
%               - BSN  为基站个数，3 < BSN <= 7；
%               - MSP  为移动台的初始位置, MSx, MSy均为[0,1]之间的数；
%                      特别要注意服务小区与MS之间的关系，MS的位置不能越界。
%               - Noise 测距误差方差。
%               - R    为小区半径，单位(meter)；
%               - X    为移动台经算法处理后的位置.
%See also: EvChanAlgorithm.m


%   参数检查：
if  nargout>1,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end

% 算法开始
MS = R*MSP;
BS = R*NetworkTop(BSN);

% MeaDist
for i = 1: BSN,
    MeaDist(i) = sqrt((MS(1) - BS(1,i))^2 + (MS(2) - BS(2,i))^2) + MeaNoise(Noise);
end

% Chan、Taylor估计位置
EMSC = ChanAlgorithm_d(BSN, MSP, R, Noise,MeaDist);
EMSTC = TaylorAlgorithm_d(BSN, EMSC'/R, R, Noise, MeaDist);
EMSTC = EMSTC';

% 残差
ResC = Residual(MeaDist,EMSC,BSN,R);
ResTC = Residual(MeaDist,EMSTC,BSN,R);

% 估计位置
EMS = (EMSC/ResC + EMSTC/ResTC)/(1/ResC + 1/ResTC);

% 结果输出
if nargout == 1,
    X = EMS;
elseif nargout == 0,
    disp(EMS);
end