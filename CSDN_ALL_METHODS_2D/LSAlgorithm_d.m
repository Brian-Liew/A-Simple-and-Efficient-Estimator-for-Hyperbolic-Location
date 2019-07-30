function X = LSAlgorithm_d(BSN, MSP, R, MeaDist)
%CHANALGORITHM_D 本函数用于实现无线定位中的CHAN算法
%               - BSN  为基站个数，3 < BSN <= 7；
%               - MSP  为移动台的初始位置, MSx, MSy均为[0,1]之间的数；
%                      特别要注意服务小区与MS之间的关系，MS的位置不能越界。
%               - MeaDist 测量距离差。
%               - R    为小区半径，单位(meter)；
%               - X    为移动台经算法处理后的位置.
%See also: ChanAlgorithm_d.m


%   参数检查：
if  nargout>1,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end


% 算法开始：
BS = R*NetworkTop(BSN);
MS = R*MSP;

% 噪声功率：
Q = 0.5*eye(BSN-1); % TDOA测量误差的协方差矩阵

% LS：
% Ri,Ki
K1 = 0;
for i = 1: BSN-1,
    R(i) = MeaDist(i+1) - MeaDist(1);
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

% 由（14b）给出B的估计值：
Za0 = pinv(Ga'*pinv(Q)*Ga)*Ga'*pinv(Q)*h';

% 输出:
out = [Za0(1),Za0(2)];

if nargout == 1,
    X = out;
elseif nargout == 0,
    disp(out);
end