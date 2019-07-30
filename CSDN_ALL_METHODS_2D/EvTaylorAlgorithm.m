function X = EvTaylorAlgorithm(BSN, MSP, Radius, Noise)
%
%TAYLORALGORITHM_D 本函数用于实现无线定位中的TAYLOR算法
%               - BSN  为基站个数，3 < BSN <= 7；
%               - MSP  为移动台的初始位置, MSx, MSy均为[0,1]之间的数；
%                      特别要注意服务小区与MS之间的关系，MS的位置不能越界。
%               - Noise 测距误差方差．
%               - R    为小区半径，单位(meter)；
%               - MeaDist 测量距离。
%               - X    为移动台经算法处理后的位置.
%See also: TaylorAlgorithm_d.m


% 参数检查：
if  nargout ~= 1 & nargout ~= 0,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end

% 初始参数：
BS = Radius*NetworkTop(BSN);
MSi = Radius*MSP;

% MeaDist
for i = 1: BSN,
    MeaDist(i) = sqrt((MSi(1) - BS(1,i))^2 + (MSi(2) - BS(2,i))^2) + MeaNoise(Noise);
end

MS = LSAlgorithm_d(BSN, MSP, Radius, MeaDist);

% TDOA协方差矩阵Q：
c = 3*10^8; % 无线电波传播速度
Dev = Noise/(c*c); % TDOA测量误差方差
Q = 0.5*eye(BSN-1); % TDOA测量误差的协方差矩阵
    
% 初始估计位置：
iEP = MS;

% h0:
for i = 1: BSN-1,
    h0(i) = MeaDist(i+1) - MeaDist(1);
end

alpha = 100;

% 算法开始：
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

% 结果输出:
if nargout == 1,
    X = abs(EP);
elseif nargout == 0,
    disp(abs(EP));
end