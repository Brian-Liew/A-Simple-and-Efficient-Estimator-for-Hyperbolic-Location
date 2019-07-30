function compare()
theta=0:1:100;
N=length(theta);
x=theta;

y=0*theta+20;

noi=1;
noise=3;
a1(25,2)=0;
b1(25,2)=0;
c1(25,2)=0;
d1(25,2)=0;
e1(25,2)=0;
%for noi=1:40;
M=1;
for i=1:M
a=chan2(noise);
%figure;

%figure;
b=Taylorforline(noise);
%figure;
c=line1(noise);
d=drop(noise);
e=offset(noise);
a1=a1+a;
b1=b1+b;
c1=c1+c;
d1=d1+d;
e1=e1+e;
end
a1=a1/M;
b1+b1/M;
c1=c1/M;
d1=d1/M;
e1=e1/M;

%end
%figure
%i=1:40;
figure;
plot(x,y,'-k');hold on;

axis([0 100 0 100]);
plot(a1(:,1),a1(:,2),'+-b');
hold on;
plot(b1(:,1),b1(:,2),'x-r');
hold on;
legend('原始路径','chan','Taylor','location','northeast');
ylabel('y/m');
xlabel('x/m');
title('噪声标准差为3m的定位结果');

figure;
plot(x,y,'-k');hold on;
axis([0 100 0 100]);
plot(c1(:,1),c1(:,2),'+--b');

hold on;
plot(d1(:,1),d1(:,2),'x:c');

hold on;
plot(e1(:,1),e1(:,2),'-.m');
hold on;
legend('原始路径','标准卡尔曼算法','测量值丢弃法','整体偏移法','location','northeast');
ylabel('y/m');
xlabel('x/m');
title('噪声标准差为3m的定位结果');
