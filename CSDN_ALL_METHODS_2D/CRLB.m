function X = CRLB( BSN, MSP, R, noise )
%   Cramer-Rao Lower Bound 是无偏估计的理论下界
%   CRML参数说明:
%       BSN: 基站的个数；
%       MS:  移动台的位置，其中MSx/MSy均在[0，1]；
%       R:   小区半径。
%
%       X:   CRLB是一矩阵的形式，主对角线上第一个元素表示横坐标的CRLB；
%            第二个元素表示定位精度的CRLB。
%   See also CRLB.m


%   参数检查：
if  nargout>1,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end
% if BSN <= 3,
%     error('The number of BSs must be larger than 3 for this program.');
% end
flag = size(MSP);
if flag(1)~=1 | flag(2)~=2,
    error('Wrong position vector!');
end

%   初始参数设置:
%     BSN = 4;
%     MS  = [0.5, 0.8];
%     R   = 5000;
    BS  = [0, sqrt(3)*R, 0.5*sqrt(3)*R, -0.5*sqrt(3)*R, -sqrt(3)*R, -0.5*sqrt(3)*R, 0.5*sqrt(3)*R;
           0, 0,         1.5*R,         1.5*R,          0,          -1.5*R,         1.5*R];
    c   = 3*10^8;
    Dev = noise^2/(c*c);
    Q   = 0.5*Dev*(eye(BSN-1) + ones(BSN-1));
    %Q = Dev*eye(BSN-1);
%   算法流程：
    MS = R*MSP;
    
    BSR_1  = sqrt(MS(1)*MS(1) + MS(2)*MS(2)); % R(1)
    
    % R(2) --- R(BSN)
    for i = 1: BSN-1,
        BSR(i) = sqrt((BS(1, i+1) - MS(1))*(BS(1, i+1) - MS(1)) + (BS(2, i+1) - MS(2))*(BS(2, i+1) - MS(2)));
    end
    % Ga0
    for i = 1:BSN -1,
        Ga(i,1) = -BS(1, i+1);
        Ga(i,2) = -BS(2, i+1);
        Ga(i,3) = -BSR(i);
    end
    % Ga'
    mGa = [1, 0; 0, 1; 1, 1];
    % B
    B = zeros(BSN-1, BSN-1);
    for i = 1: BSN-1,
        B(i, i) = BSR(i);
    end
    % B'
    mB = [MS(1), 0, 0; 0, MS(2), 0; 0, 0, BSR_1];
    % B''
    mmB = [MS(1), 0; 0, MS(2);];
    
%   输出
    Crlb = c*c*inv(mmB*mGa'*inv(mB)*Ga'*inv(B)*inv(Q)*inv(B)*Ga*inv(mB)*mGa*mmB);
    
    out = sqrt(Crlb(1,1) + Crlb(2,2));
    
    if nargout == 1, 
        X = out;
    elseif nargout == 0,
        disp( out );
    end