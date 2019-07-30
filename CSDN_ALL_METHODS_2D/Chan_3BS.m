function X = Chan_3BS(MSP,R,Noise)
%   Chan �㷨������3BS��MS���ж�λ��
%   CHAN_3BS:
%       ����˵����
%       Noise:   �������.
%           R��  С���뾶.
%  Also see: Chan_3BS.


%  �������:
if nargout ~=1,
    error('Too many output arguments!');
end 
if nargin ~= 3,
    error('input arguments error!');
end

%  �㷨��ʼ
MS = R*MSP;
BS = R*NetworkTop(3);

% A����:
X21 = BS(1,2) - BS(1,1);
X31 = BS(1,3) - BS(1,1);
Y21 = BS(2,2) - BS(2,1);
Y31 = BS(2,3) - BS(2,1);
A = inv([X21,Y21;X31,Y31]);

% B����:
R1 = sqrt((BS(1,1) - MS(1))^2 + (BS(2,1) - MS(2))^2);
R2 = sqrt((BS(1,2) - MS(1))^2 + (BS(2,2) - MS(2))^2);
R3 = sqrt((BS(1,3) - MS(1))^2 + (BS(2,3) - MS(2))^2);

R21 = R2 - R1 + MeaNoise(Noise);  % ��Ҫ������
R31 = R3 - R1 + MeaNoise(Noise);
B = [R21;R31];

% C����:
K1 = BS(1,1)^2 + BS(2,1)^2;
K2 = BS(1,2)^2 + BS(2,2)^2;
K3 = BS(1,3)^2 + BS(2,3)^2;
C = 0.5*[R21^2 - K2 + K1; R31^2 - K3 + K1];

% һԪ���η��̵�ϵ����
a = B'*A'*A*B - 1;
b = B'*A'*A*C + C'*A'*A*B;
c = C'*A'*A*C;

% ���̵���������
root1 = abs((-b + sqrt(b^2 - 4*a*c))/(2*a));
root2 = abs((-b - sqrt(b^2 - 4*a*c))/(2*a));

% ���鷽�̵ĸ���
if root1 < R,
    EMS = -A*(B*root1 + C);
else
    EMS = -A*(B*root2 + C);
end

% ��������
if nargout == 1,
    X = EMS;
else
    disp(EMS);
end