function X = SXAlgorithm(BSN, MSP, Radius, Noise)
%SXALGORITHM ����������ʵ�����߶�λ�е�SX�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - R    ΪС���뾶����λ(meter)��
%               - Noise ������
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: SXAlgorithm.m


% ������飺
if  nargout>1,
    error('Too many output arguments.');
end
if nargin<2 | nargin>4,
    error('Wrong number of input arguments.');
end
if BSN < 3,
    error('The number of BSs must be larger than 3 for this program.');
end
flag = size(MSP);
if flag(1)~=1 | flag(2)~=2,
    error('Wrong position vector!');
end

% ��ʼ������
BS = Radius*NetworkTop(BSN);
MS = Radius*MSP;

% ��Ȩ����W��
% c = 3*10^8; % ���ߵ粨�����ٶ�
% Dev = Noise/(c*c); % TDOA��������
W = eye(BSN-1);
    
% Ri1:
R1 = sqrt(MS(1)^2 + MS(2)^2);
for i = 1: BSN-1,
    R(i) = sqrt((BS(1,i+1) - MS(1))^2 + (BS(2,i+1) - MS(2))^2);
end
for i = 1: BSN-1,
    Ri1(i) = R(i) - R1 + Noise*randn(1);
end
    
% delt:
for i = 1: BSN-1,
    K(i) = BS(1,i+1)^2 + BS(2,i+1)^2;
end
for i = 1: BSN-1,
    delt(i) = K(i) - Ri1(i)^2;
end
    
% S:
for i = 1: BSN-1,
    S(i, 1) = BS(1, i+1);
    S(i, 2) = BS(2, i+1);
end

% Sw:
Sw = inv(S'*W*S)*S'*W;

% a:
a = 4 - 4*Ri1*Sw'*Sw*Ri1';

% b:
b = 4*Ri1*Sw'*Sw*delt';

% c:
c = -delt*Sw'*Sw*delt';

% root:
root1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
root2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);

% Za:
% if root1 > 0,
%     Za = 0.5*Sw*(delt' - 2*root1^2);
% elseif root2 > 0,
%     Za = 0.5*Sw*(delt' - 2*root2^2);
% elseif abs(root1) < abs(root2),
%     Za = 0.5*Sw*(delt' - 2*abs(root1)^2);
% else
%     Za = 0.5*Sw*(delt' - 2*abs(root2)^2);
% end

if root1 > 0,
    Za = 0.5*Sw*(delt' - 2*root1^2);
else
    Za = 0.5*Sw*(delt' - 2*root2^2);
end

% if abs(root1) < abs(root2),
%     Za = 0.5*Sw*(delt' - 2*abs(root1)^2);
% else
%     Za = 0.5*Sw*(delt' - 2*abs(root2)^2);
% end


% ���:
if nargout == 1,
    X = Za;
elseif nargout == 0,
    disp(Za);
end