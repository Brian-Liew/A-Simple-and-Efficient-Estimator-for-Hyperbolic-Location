function compare()
theta=0:1:100;
N=length(theta);
x=theta;

y=0*theta+20;

noi=10;
noise=1;

%for noi=1:40;
a=chan2(noi,noise);
%figure;

%figure;
b=Taylorforline(noi,noise);
%figure;
c=line1(noi,noise);
d=drop1(noi,noise);
e=offset1(noi,noise,1.1);
%end
%figure
%i=1:40;
figure;
plot(x,y,'-k');hold on;

axis([0 100 0 100]);
plot(a(:,1),a(:,2),'+-b');
hold on;
plot(b(:,1),b(:,2),'x-r');
hold on;
legend('原始路径','chan','Taylor','location','northeast');
ylabel('y/m');
xlabel('x/m');
title('NLOS误差均值为10m的定位结果');

figure;
plot(x,y,'-k');hold on;
axis([0 100 0 100]);
plot(c(:,1),c(:,2),'+--b');

hold on;
plot(d(:,1),d(:,2),'x:c');

hold on;
plot(e(:,1),e(:,2),'-.m');
hold on;
legend('原始路径','标准卡尔曼算法','测量值丢弃法','整体偏移法','location','northeast');
ylabel('y/m');
xlabel('x/m');
title('NLOS误差均值为10m的定位结果');
