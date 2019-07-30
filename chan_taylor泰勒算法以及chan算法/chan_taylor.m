function [x,y,z]=chan_taylor(x0,y0,z0)
X=[0 -15 0 15];
Y=[0 15 15 15];
Z=[0 10 10 10];
BS=[X;Y;Z;];
x0=input('x0=');y0=input('y0=');z0=input('z0=');
test=[x0;y0;z0];
disp('基站坐标');disp(BS);
disp('假设目标坐标');disp(test);
figure(1)
plot3(X,Y,Z,'rp',x0,y0,z0,'g*');
text(X(1),Y(1),Z(1),'BS1');text(X(2),Y(2),Z(2),'BS2');text(X(3),Y(3),Z(3),'BS3');text(X(4),Y(4),Z(4),'BS4')
text(x0,y0,z0,'假设目标');
xlabel('x轴');ylabel('y轴');zlabel('z轴');
grid on;
axis([-30 30 -30 30 -30 30]);
box;
for i=1:4            
    r=sqrt((X-x0).^2+(Y-y0).^2+(Z-z0).^2);
    k=sqrt(X.^2+Y.^2+Z.^2);
end
disp('目标到基站距离r');disp(r);
disp('基站到原点距离k');disp(k);
for i=2:4      %有效代码
    H(i-1)=(k(i).^2-k(1).^2-(r(i)-r(1)).^2)/2;
    G((i-1),1)=X(i)-X(1);G((i-1),2)=Y(i)-Y(1);G((i-1),3)=Z(i)-Z(1);G((i-1),4)=r(i)-r(1);
end   %有效代码
disp(H');disp(G);
target=G\H';
disp(target);
hold on;
plot3(target(1),target(2),target(3),'ms');
