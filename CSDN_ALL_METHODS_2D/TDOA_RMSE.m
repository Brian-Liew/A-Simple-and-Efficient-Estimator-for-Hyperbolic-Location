function X = TDOA_RMSE(MS, EMS)
%
%TDOA_RMSE 本函数用于实现无线定位精度RMSE的计算
%               - MS  为移动台的真实位置；
%               - EMS 为移动台的估计位置。
%See also: TDOA_RMSE.m


% 参数检查：
if  nargout ~= 1 & nargout ~= 0,
    error('Too many output arguments.');
end
if nargin ~= 2,
    error('Wrong number of input arguments.');
end

% 算法开始：
[n,m] = size(MS);

% 计算均方误差总和：
sum = 0;
for i = 1:n,
    sum = (MS(i,1) - EMS(i,1))^2 + (MS(i,2) - EMS(i,2))^2 + sum;
end

% RMSE:
RMSE = sqrt(sum/n);

% 结果输出：
if nargout == 1,
    X = RMSE;
else
    disp(RMSE);
end
    