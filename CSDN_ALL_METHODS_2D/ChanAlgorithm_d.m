function X = ChanAlgorithm_d(BSN, MSP, R, Noise, MeaDist)
%CHANALGORITHM_D 本函数用于实现无线定位中的CHAN算法
%               - BSN  为基站个数，3 < BSN <= 7；
%               - MSP  为移动台的初始位置, MSx, MSy均为[0,1]之间的数；
%                      特别要注意服务小区与MS之间的关系，MS的位置不能越界。
%               - Noise 测距误差方差。
%               - MeaDist 测量距离差。
%               - R    为小区半径，单位(meter)；
%               - X    为移动台经算法处理后的位置.
%See also: ChanAlgorithm_d.m


%   参数检查：
if  nargout>1,
    error('Too many output arguments.');
end
if nargin ~= 5,
    error('Wrong number of input arguments.');
end


% 算法开始：
BS = R*NetworkTop(BSN);
MS = R*MSP;

% 噪声功率：
c = 3*10^8; % 无线电波传播速度
Q = 0.5*eye(BSN-1); % TDOA测量误差的协方差矩阵

% 第一次LS：
% Ri
R1 = sqrt(MS(1)^2 + MS(2)^2);
K1 = 0;
for i = 1: BSN-1,
    R0(i) = sqrt((BS(1,i+1) - MS(1))^2 + (BS(2,i+1) - MS(2))^2);
end

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

% 利用这个粗略估计值计算B：
B = eye(BSN-1);
for i = 1: BSN-1,
    B(i,i) = sqrt((BS(1,i+1) - Za0(1))^2 + (BS(2,i+1) - Za0(2))^2);
end

% mFI:
mFI = B*Q*B;

% 第一次LS结果：
Za1 = pinv(Ga'*pinv(mFI)*Ga)*Ga'*pinv(mFI)*h';

if Za1(3) < 0,
    Za1(3) = abs(Za1(3));
end

%***************************************************************

% 第二次LS：
% 第一次LS结果的协方差：
CovZa = pinv(Ga'*pinv(mFI)*Ga);

% sB：
sB = eye(3);
for i = 1: 3,
    sB(i,i) = Za1(i);
end;

% sFI：
sFI = 4*sB*CovZa*sB;

% sGa：
sGa = [1, 0; 0, 1; 1, 1];

% sh
sh  = [Za1(1)^2; Za1(2)^2; Za1(3)^2];

% 第二次LS结果：
Za2 = pinv(sGa'*pinv(sFI)*sGa)*sGa'*pinv(sFI)*sh;

% 输出:
Za = sqrt(abs(Za2));

out = Za;

if nargout == 1,
    X = out;
elseif nargout == 0,
    disp(out);
end