function X = MeaNoise(Noise)
%  本程序生成测距噪声，服从高斯分布
%  MeaNoise
%    参数说明：
%        Noise:   高斯分布方差乘以光速平方的结果
%  Also see: MeaNoise.


%  参数检测:
if nargout>1,
    error('Too many output arguments!');
end 
if nargin ~= 1
    error('input arguments error!');
end

% 测距误差方差：
Dev = Noise;

% 测距误差：
X = sqrt(Dev)*randn(1);

% 结果输出：
if nargout == 1,
    X;
else
    disp(X);
end
