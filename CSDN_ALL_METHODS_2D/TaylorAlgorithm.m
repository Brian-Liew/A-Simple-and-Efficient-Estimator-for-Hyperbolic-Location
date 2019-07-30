function X = TaylorAlgorithm(BSN, MSP, Radius, Noise)
%
%TAYLORALGORITHM 本函数用于实现无线定位中的TAYLOR算法
%               - BSN  为基站个数，3 < BSN <= 7；
%               - MSP  为移动台的初始位置, MSx, MSy均为[0,1]之间的数；
%                      特别要注意服务小区与MS之间的关系，MS的位置不能越界。
%               - Noise 测距误差方差．
%               - R    为小区半径，单位(meter)；
%               - X    为移动台经算法处理后的位置.
%See also: TaylorAlgorithm.m


% 参数检查：
if  nargout ~= 1 & nargout ~= 0,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end

% 初始参数：
BS = Radius*NetworkTop(BSN);
MS = Radius*MSP;

% TDOA协方差矩阵Q：
Q = eye(BSN-1);
    
% 初始估计位置：
iEP = MS;

% h0:
for i = 1: BSN,
    MeaDist(i) = sqrt((MS(1) - BS(1,i))^2 + (MS(2) - BS(2,i))^2);
end
for i = 1: BSN-1,
    h0(i) = MeaDist(i+1) - MeaDist(1) + Noise*randn(1);   %TDOA测量值
end

% 算法开始：
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

% 结果输出:
if nargout == 1,
    X = EP;
elseif nargout == 0,
    disp(EP);
end