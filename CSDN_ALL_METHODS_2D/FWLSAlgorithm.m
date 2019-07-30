function X = FWLSAlgorithm(BSN, MSP, R)
%
%FWLSALGORITHM 本函数用于实现无线定位中的Friedlander算法
%               - BSN  为基站个数，3 < BSN <= 7；
%               - MSP  为移动台的初始位置, MSx, MSy均为[0,1]之间的数；
%                      特别要注意服务小区与MS之间的关系，MS的位置不能越界。
%               - R    为小区半径，单位(meter)；
%               - X    为移动台经算法处理后的位置.
%See also: FWLSAlgorithm.m


%   参数检查：
if  nargout>1,
    error('Too many output arguments.');
end
if nargin<2 | nargin>3,
    error('Wrong number of input arguments.');
end
if BSN <= 3,
    error('The number of BSs must be larger than 3 for this program.');
end
flag = size(MSP);
if flag(1)~=1 | flag(2)~=2,
    error('Wrong position vector!');
end

%   初始参数：
    BSN = 4;
    MSP = [0.5, 0.5];
    R = 5000; % 
    BS = [0, sqrt(3)*R, 0.5*sqrt(3)*R, -0.5*sqrt(3)*R, -sqrt(3)*R, -0.5*sqrt(3)*R, 0.5*sqrt(3)*R;
          0, 0,         1.5*R,         1.5*R,          0,          -1.5*R,         -1.5*R];
    MS = R*MSP;
    c = 3*10^8; % 无线电波传播速度
    Dev = 900/(c*c); % TDOA测量误差方差
    Q = 0.5*Dev*(eye(BSN -1)+ones(BSN -1)); % TDOA测量误差的协方差矩阵
    
    % S
    for i = 1: BSN-1,
        S(i, 1) = BS(1, i+1);
        S(i, 2) = BS(2, i+1);
    end
    
    % N
        % Z
        Z = ones(BSN-1);
        for i = 1: BSN-1,
            for j = 1: BSN-1,
               if j == i | j>i+1,
                   Z(i, j) = 0;
               end
            end
        end
        % D
        for i = 1: BSN,
            R(i) = sqrt((BS(1, i) - MS(1))^2 + (BS(2, i) - MS(2))^2);
        end
        a = sqrt(0.5*Dev)*randn(1);
        for i = 1: BSN-1,
            b = sqrt(0.5*Dev)*randn(1);
            Ri1(i) = R(i+1) - R(1) + a + b;
        end
        D = eye(BSN-1);
        for i = 1: BSN-1,
            D(i, i) = Ri1(i);
        end
        % I
        I = eye(BSN-1);
    N = (I-Z)*D;
    
    % u
        % K
        for i = 1: BSN-1,
            K(i) = BS(1, i+1)^2 + BS(2, i+1)^2;
        end
    for i = 1: BSN-1,
        u(i) = 0.5*(K(i) - Ri1(i)^2);
    end
    
    Za = inv(S'*N*Q*N'*S)*S'*N*Q*N'*u';
%    输出:
if nargout == 1,
    X = Za;
elseif nargout == 0,
    disp(Za);
end