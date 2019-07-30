function compare()


noi=1;
noise=5;
a1=[0 0 0 0 0 ];
b1=[0 0 0 0 0 ];
c1=[0 0 0 0 0 ];
d1=[0 0 0 0 0 ];
e1=[0 0 0 0 0 ];
N=30;
%for noi=1:40;
for j=1:N
for i=5:5:25
%a(i/5)=chan2(i);
%figure;

%figure;
b(i/5)=dropoffset(i,noise,1.1);
%figure;
c(i/5)=line1(i,noise);
d(i/5)=drop1(i,noise);
e(i/5)=offset1(i,noise,1.1);
end
%a1=a1+a;
b1=b1+b;
c1=c1+c;
d1=d1+d;
e1=e1+e;
end
%a1=a1/N;
b1=b1/N;
c1=c1/N;
d1=d1/N;
e1=e1/N;
%figure
%i=1:40;
figure;
i=5:5:25;
%plot(i,a1(i*2),'+-g');
%hold on;
plot(i,b1(i/5),'x-r');
hold on;


plot(i,c1(i/5),'+--b');

hold on;
plot(i,d1(i/5),'x:k');

hold on;
plot(i,e1(i/5),'-.m');
hold on;
legend('Taylor','标准卡尔曼算法','测量值丢弃法','整体偏移法','location','northwest');
%axis([0 3 0 5]);
ylabel('测量误差均方根/m');
xlabel('NLOS噪声均值/m');
title('不同NLOS噪声均值的定位比较');
