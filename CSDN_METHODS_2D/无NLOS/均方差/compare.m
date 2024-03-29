function compare()


noi=1;
noise=3;
a1=[0 0 0 0 0 0 0 0];
b1=[0 0 0 0 0 0 0 0];
c1=[0 0 0 0 0 0 0 0];
d1=[0 0 0 0 0 0 0 0];
e1=[0 0 0 0 0 0 0 0];
N=20;
%for noi=1:40;
for j=1:N
for i=0.5:0.5:4
a(i*2)=chan2(i);
%figure;

%figure;
b(i*2)=Taylorforline(i);
%figure;
c(i*2)=line1(i);
d(i*2)=drop(i);
e(i*2)=offset(i);
end
a1=a1+a;
b1=b1+b;
c1=c1+c;
d1=d1+d;
e1=e1+e;
end
a1=a1/N;
b1=b1/N;
c1=c1/N;
d1=d1/N;
e1=e1/N;
%figure
%i=1:40;
figure;
i=0.5:0.5:4;
%plot(i,a1(i*2),'+-g');
%hold on;
plot(i,b1(i*2),'x-r');
hold on;


plot(i,c1(i*2),'+--b');

hold on;
plot(i,d1(i*2),'x:k');

hold on;
plot(i,e1(i*2),'-.m');
hold on;
legend('Taylor','标准卡尔曼算法','测量值丢弃法','整体偏移法','location','northwest');
%axis([0 3 0 5]);
ylabel('测量误差均方根/m');
xlabel('噪声误差均方根/m');
title('不同噪声均方根的定位比较');
