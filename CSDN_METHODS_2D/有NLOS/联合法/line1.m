function X=line(noi,noise)
rmse=0;
%target trajectory :circle,(x,y)
%theta = 0:0.1:2*pi+0.2;
%N=length(theta);
%x = 0.4*sin(theta)+0.5;
%y = 0.4*cos(theta)+0.5;
%plot(x,y,'-r');hold on;
theta=0:4:100;
N=length(theta);
x=theta;

y=0*theta+20;
plot(x,y,'-r');hold on;
%linearized mode
A = [1 1 0 0;...
    0 1 0 0;...
    0 0 1 1;...
    0 0 0 1];
C = cell(1,N);
%observation data
for k = 1 : N
    x(k) = x(k)+0*randn(1);%add system noise ,u(k)         先忽略
    y(k) = y(k)+0*randn(1);
    z1(k) = sqrt(x(k)^2+y(k)^2);
    z2(k) = sqrt((x(k)-100)^2+y(k)^2);
    z3(k) = sqrt(x(k)^2+(y(k)-100)^2);
    z4(k) = sqrt((x(k)-100)^2+(y(k)-100)^2);
end
z = cell(1,N);
%for k = 1 : N
%    z{k} = [z1(k);z2(k);z3(k);z4(k)] +[1;10;1;1];%add observation noise v(k)
%end
for k = 1 : N
   if k<=45/4||k>65/4
    z{k} = [z1(k);z2(k);z3(k);z4(k)] +noise*randn(4,1);%[1;1;1;1];%add observation noise v(k)
   else
    z{k} = [z1(k);z2(k);z3(k);z4(k)] +noise*randn(4,1)+[0;0;noi+0.1*randn(1);0];%[1;1;1;1];
   end
end

%estimated state:x_s
x_s = cell(1,N);
x_s{1} = [x(1);0;y(1);0];
P = cell(1,N);
P{1} = eye(4);
Q = 0.01*eye(4);       %设为很小值，视作不影响
R = 0.01*eye(4);
x_sTime = cell(1,N);
for k = 1 : N-1
    %C{k} = [x_s{k}(1)/sqrt(x_s{k}(1)^2+x_s{k}(3)^2) 0 x_s{k}(3)/sqrt(x_s{k}(1)^2+x_s{k}(3)^2) 0;...
     %       (2*x_s{k}(1)-2)/sqrt(2*((x_s{k}(1)-1)^2+x_s{k}(3)^2)) 0 x_s{k}(3)/sqrt((x_s{k}(1)-1)^2+x_s{k}(3)^2) 0;...
    %        x_s{k}(1)/sqrt((x_s{k}(3)-1)^2+x_s{k}(1)^2) 0 (2*x_s{k}(3)-2)/sqrt(2*((x_s{k}(3)-1)^2+x_s{k}(1)^2)) 0];
    C{k} = [ 1/(x_s{k}(1)^2+x_s{k}(3)^2)^(1/2)*x_s{k}(1),0, 1/(x_s{k}(1)^2+x_s{k}(3)^2)^(1/2)*x_s{k}(3),0;...
           1/2/(x_s{k}(1)^2-200*x_s{k}(1)+10000+x_s{k}(3)^2)^(1/2)*(2*x_s{k}(1)-200), 0, 1/(x_s{k}(1)^2-200*x_s{k}(1)+10000+x_s{k}(3)^2)^(1/2)*x_s{k}(3),0;...
           1/(x_s{k}(1)^2+x_s{k}(3)^2-200*x_s{k}(3)+10000)^(1/2)*x_s{k}(1), 0,1/2/(x_s{k}(1)^2+x_s{k}(3)^2-200*x_s{k}(3)+10000)^(1/2)*(2*x_s{k}(3)-200),0;...
           1/2/(x_s{k}(1)^2-200*x_s{k}(1)+20000+x_s{k}(3)^2-200*x_s{k}(3))^(1/2)*(2*x_s{k}(1)-200),0, 1/2/(x_s{k}(1)^2-200*x_s{k}(1)+20000+x_s{k}(3)^2-200*x_s{k}(3))^(1/2)*(2*x_s{k}(3)-200),0];
    K{k} = P{k}*C{k}'*inv(R+C{k}*P{k}*C{k}');
    z_e = [sqrt(x_s{k}(1)^2+x_s{k}(3)^2);sqrt((x_s{k}(1)-100)^2+x_s{k}(3)^2);sqrt(x_s{k}(1)^2+(x_s{k}(3)-100)^2);sqrt((x_s{k}(1)-100)^2+(x_s{k}(3)-100)^2)];
    x_sTime{k}=x_s{k}+ K{k}*(z{k}-z_e);%measurement update
    x_s{k+1} = A*x_sTime{k};%time update
    P{k+1} = A*P{k}*A' + Q - A*P{k}*C{k}'*inv(C{k}*P{k}*C{k}'+R)*C{k}*P{k}*A';
    
%     x_e(k) = x_sTime{k}(1);
%     y_e(k) = x_sTime{k}(3);
%     plot(x(k),y(k),'-^r');hold on;
%     plot(x_e(k),y_e(k),'*k');pause(0.1);hold on;
end
for k = 1 : N-1
    x_e(k) = x_sTime{k}(1);
    y_e(k) = x_sTime{k}(3);
     rmse=rmse+sqrt((x_e(k)-x(k))^2+( y_e(k)-y(k))^2);
end
X=rmse/(N-1);
%X=[x_e',y_e'];
%plot(x_e,y_e,'-k');
%axis([0 100 0 100]);
%xlabel('x');
%ylabel('y');
%title('Tracking Performance');
%legend('actual orbit','estimated orbit')
    
    
    
    
    
    
    
    
    
    
    
    
    
    