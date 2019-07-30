function Estimation=Chan_TDOA_3D(BSNum,X,Xb,Real_ms,P)
c=3e8;  %Xb为基站坐标，bsnum为基站数量
Rb = sqrt((Real_ms(1) - Xb(1))^2+(Real_ms(2) - Xb(2))^2+(Real_ms(3)-Xb(3))^2);%移动台到基站的真实距离
%Q = eye(BSNum)*(delta0^2);
Q=(eye(BSNum)+ones(BSNum))/2;
%Q=chol(Q_a);
%P= normrnd(0,delta0,1,BSNum);%产生TDOA测量误差
Kb = sum(Xb.^2);
R = zeros(BSNum,1);
for i=1:BSNum                %产生TDOA测量值
    R(i) = -Rb+sqrt((Real_ms(1)- X(i,1))^2+(Real_ms(2) - X(i,2))^2+(Real_ms(3)-X(i,3))^2)+P(i);
end
Pbs = [Xb;X];
N = size(Pbs,1);
K = zeros(1,N);
K = Pbs(:,1).^2 + Pbs(:,2).^2+Pbs(:,3).^2;
ha = 0.5*(R.^2-K(2:N)+K(1));
Ga = -[Pbs(2:N,1)-Xb(1) Pbs(2:N,2)-Xb(2) Pbs(2:N,3)-Xb(3) R];

        %计算第一次WLS估计结果（远距算法）
Za=inv(Ga.'*inv(Q)*Ga)*Ga.'*inv(Q)*ha;
for i=1:BSNum
    Ba(i,i)=sqrt((Za(1)-X(i,1))^2+(Za(2)-X(i,2))^2+(Za(3)-X(i,3))^2);
end
Fa = c^2*Ba*Q*Ba;
Za1 = inv(Ga.'*inv(Fa)*Ga)*Ga.'*inv(Fa)*ha; %第一次WLS的估计结果,作为参考作最终定位结果的判决
if Za1(4)<0
   Za1(4)=abs(Za1(4));
end
Zacov = inv(Ga.'*inv(Fa)*Ga);
        %第一次ＷＬＳ计算（近距算法）
ha2=[(Za1(1)-Xb(1))^2;(Za1(2)-Xb(2))^2;(Za1(3)-Xb(3))^2;Za1(4)^2];
Ga2=[1 0 0;0 1 0;0 0 1;1 1 1];
Ba2=diag([Za1(1)-Xb(1),Za1(2)-Xb(2),Za1(3)-Xb(3),Za1(4)]);
Fa2=4*Ba2*Zacov*Ba2;

Za2=inv(Ga2.'*inv(Fa2)*Ga2)*Ga2.'*inv(Fa2)*ha2;

%Za2=inv(Ga2'*inv(Ba2)*Ga2*inv(Q)*Ga2'*inv(Ba2)*Ga2)*Ga2'*inv(Ba2)*Ga2*inv(Q)*Ga2'*inv(Ba2)*ha2;
%Zacov2 = inv(Ga2.'*inv(Fa2)*Ga2)/10;
Za2=abs(Za2);
out=diag(sign(Za1(1:3)-Pbs(1,:).'))*sqrt(Za2)+Pbs(1,:).';
Estimation=out.';
%-------------------------------------------------------以上chan算法