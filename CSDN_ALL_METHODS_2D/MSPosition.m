function X = MSPosition()
%  本程序在1/12小区内随机产生MS的位置，假设小区半径为1
%  MSPosition
%    参数说明：
%       无参数。
%  Also see: MSPosition.


%  参数检测:
 if  nargout>1,
        error('Too many output arguments.');
 end

% 随机产生移动台位置：
x = sqrt(3)*rand(1)/2;
y = sqrt(3)*x*rand(1)/3;

% 结果输出：
if nargout == 1,
    X = [x, y];
else
    disp([x, y]);
end

