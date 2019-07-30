function X = TaylorAlgorithm(noi,noise)
%
%TAYLORALGORITHM ����������ʵ�����߶�λ�е�TAYLOR�㷨
%               - BSN  Ϊ��վ������3 < BSN <= 7��
%               - MSP  Ϊ�ƶ�̨�ĳ�ʼλ��, MSx, MSy��Ϊ[0,1]֮�������
%                      �ر�Ҫע�����С����MS֮��Ĺ�ϵ��MS��λ�ò���Խ�硣
%               - Noise ������
%               - R    ΪС���뾶����λ(meter)��
%               - X    Ϊ�ƶ�̨���㷨������λ��.
%See also: TaylorAlgorithm.m
rmse=0;

BSN=4;
BS=[0,100, 0 ,100;
    0, 0 ,100,100];
theta=1:1:100;
N=length(theta);
x=theta;
y=0*theta+20;
plot(x,y,'-r');hold on;

% TDOAЭ�������Q��
Q = 0.01*eye(BSN-1);
for k=1:N
    MS=[x(k),y(k)];    
% ��ʼ����λ�ã�
iEP = MS;

% h0:
for i = 1: BSN
    if 20<k&&k<30 && i==2
    MeaDist(i) = sqrt((MS(1) - BS(1,i))^2 + (MS(2) - BS(2,i))^2)+noise*randn(1);%+10;
  
    else
    MeaDist(i) = sqrt((MS(1) - BS(1,i))^2 + (MS(2) - BS(2,i))^2)+noise*randn(1);%;
    end
end
for i = 1: BSN-1,
    h0(i) = MeaDist(i+1) - MeaDist(1)+ noise*randn(1);%Noise %*randn(1);   %TDOA����ֵ
end

% �㷨��ʼ��
for n = 1: 10,
    % Rn:
    R1 = sqrt(iEP(1)^2 + iEP(2)^2);
    for i =1: BSN-1,
        R(i) = sqrt((iEP(1) - BS(1,i+1))^2 + (iEP(2) - BS(2,i+1))^2);        
    end
    
    % ht:
    for i = 1: BSN-1,
        h(i) = h0(i) - (R(i) - R1);
    end
    ht = h';
    
    % Gt:
    for i = 1: BSN-1,
        Gt(i, 1) = -iEP(1)/R1 - (BS(1, i+1) - iEP(2))/R(i);
        Gt(i, 2) = -iEP(2)/R1 - (BS(2, i+1) - iEP(2))/R(i);
    end
    
    % delt:
    delt = inv(Gt'*inv(Q)*Gt)*Gt'*inv(Q)*ht;
    
    EP = iEP + delt';

    iEP = EP;
end

out1(1,k)=EP(1);
out1(2,k)=EP(2);
rmse=rmse+sqrt((out1(1,k)-x(k))^2+(out1(2,k)-y(k))^2);
end
X=rmse/10;
k=1:N;
 plot(out1(1,k),out1(2,k),'-g');