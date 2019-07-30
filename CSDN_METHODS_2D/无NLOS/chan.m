%**********************基于TDOA的Chan算法**********************************
function [MS]=chan(M)
%假设移动台坐标为
ms=[500,1000];
x=ms(1);y=ms(2);%移动台真实位置
c=3*10^8;%信号传播速度

X=[0,0,-4500,-4500,0,4500,4500,4500,-4500,-9000,-4500,4500,9000];
Y=[0,5196,2598,-2598,-5916,-2598,2598,7794,7794,0,-7794,-7794,0];
%假设小区半径是3000m.X,Y分别是基站位置横纵坐标,
basestx=X(1:M);
basesty=Y(1:M); %M是参与定位的基站数目,M的取值最大是13.
N=length(basestx); %参与定位的基站数目

Standarddeviation=[30,40,50,60,70,80,90,100,110,120]; %测量误差标准差 /m


ri1=[]; %第i（i>=2）个基站到移动台距离与第一个基站的（服务基站）到移动台距离的差值。
xi1=[]; %第i个基站与第一个基站位置横坐标的差值
yi1=[]; %第i个基站与第一个基站位置纵坐标的差值
k=[];
h=[];
Ga=[];

for i=2:N
     xi1(i-1)=basestx(i)-basestx(1);
     yi1(i-1)=basesty(i)-basesty(1);
end  %对xi1,yi1进行赋值
   
for i=1:N
    k(i)=(basestx(i))^2+(basesty(i))^2;
end  %对k进行赋值

rmse=[];
for j0=1:length(Standarddeviation)
 
   for i=2:N
       ri1(i-1)=sqrt((basestx(i)-x)^2+(basesty(i)-y)^2)- sqrt((basestx(1)-x)^2+(basesty(1)-y)^2)-Standarddeviation(j0);
     %在测量参数不知道的情况下为方便仿真，假设移动台位置已知，
     %则可以知道各个基站与移动台的实际距离差。
     %由于测量有误差，这里用实际距离差加上或减去测量误差标准差来表示测到的距离差
   end

   for i=2:N
        h(i-1)=0.5*((ri1(i-1))^2-k(i)+k(1));
   end  %对h进行赋值

   for i=1:3
       for j=2:N
           switch i,
                  case 1,
                        Ga(j-1,i)=-xi1(j-1);
                  case 2,
                        Ga(j-1,i)=-yi1(j-1);
                  case 3,
                        Ga(j-1,i)=-ri1(j-1);
           end
       end
   end  %对Ga进行复制运算
 
   Q=zeros(N-1,N-1);    
   for i=1:N-1
   Q(i,i)=(Standarddeviation(j0))^2; %非常重要此处,Q为测量误差的协方差矩阵
   end

   Za=inv(Ga'*inv(Q)*Ga)*Ga'*inv(Q)*h'; %第一次估计，假设移动台到每个基站距离均相等(移动台到基站距离较远)

   B1=[];
   for i=1:N-1
       B1(i,i)=sqrt((basestx(i+1)-Za(1))^2+(basesty(i+1)-Za(2))^2);
   end  %得到移动台估计位置，则可以到各个基站的距离

       P1=c^2*B1*Q*B1; %误差矢量的协方差矩阵
       Za1=inv(Ga'*inv(P1)*Ga)*Ga'*inv(P1)*h'; %第二次估计移动台位置
       C=inv(Ga'*inv(Q)*Ga);
       
       h1=[(Za1(1)-basestx(1))^2;(Za1(2)-basesty(1))^2;(Za1(3))^2];
       Ga1=[1,0;0,1;1,1];
       r1=sqrt((basestx(1)-Za1(1))^2+(basesty(1)-Za1(2))^2); %第一个基站与移动台间的距离
       B2=[Za1(1)-basestx(1),0,0;0,Za1(2)-basesty(1),0;0,0,r1];
       P2=4*B2*C*B2;
       Za2=inv(Ga1'*inv(P2)*Ga1)*Ga1'*inv(P2)*h1;
       
       ms0=sqrt(Za2)+[basestx(1);basesty(1)]; %利用第一个基站到移动台的距离与移动台位置的关系，改善估计精度
       rmse(j0)=sqrt((ms0(1)-x)^2+(ms0(2)-y)^2);
       MS(j0,:)=ms0';
end
MS;

rmse;
figure
plot(Standarddeviation,rmse,'^--r')
axis([30,120,0,100]);
legend('chan',2);
grid on;
hold on;
ylabel('定位误差均方根/m');
xlabel('TDOA误差标准差/m');
title('TDOA下M个基站参与定位')






















    
