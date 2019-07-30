function X = FWLSAlgorithm(BSN, MSP, R)
%
%FWLSALGORITHM ����������ʵ�����߶�λ�е�Friedlander�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - R    ΪС���뾶����λ(meter)��
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: FWLSAlgorithm.m


%   ������飺
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

%   ��ʼ������
    BSN = 4;
    MSP = [0.5, 0.5];
    R = 5000; % 
    BS = [0, sqrt(3)*R, 0.5*sqrt(3)*R, -0.5*sqrt(3)*R, -sqrt(3)*R, -0.5*sqrt(3)*R, 0.5*sqrt(3)*R;
          0, 0,         1.5*R,         1.5*R,          0,          -1.5*R,         -1.5*R];
    MS = R*MSP;
    c = 3*10^8; % ���ߵ粨�����ٶ�
    Dev = 900/(c*c); % TDOA��������
    Q = 0.5*Dev*(eye(BSN -1)+ones(BSN -1)); % TDOA��������Э�������
    
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
%    ���:
if nargout == 1,
    X = Za;
elseif nargout == 0,
    disp(Za);
end