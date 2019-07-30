function [out, EE ] = Chan_3D_linear( num_of_sensor,sensor,location,sigma_d_square  )
%The function use the Chan to solve the mse of tdoa of 3D
% mes is the MSE; the function give the number of
% sensors and the location of the sensors and the source
%clear;
EE=0;
%num_of_sensor=8;
%sigma_d_square=0.0001./((3.0*10.^8).^2);
T=0.5.*ones(num_of_sensor-1);
for tt=1:num_of_sensor-1
    T(tt,tt)=1;
end
Q=sigma_d_square.*T;      %the Q matrix has been solved
for i=1:num_of_sensor
    R(i)=((sensor(1,i)-location(1)).^2+(sensor(2,i)-location(2)).^2+(sensor(3,i)-location(3)).^2).^(1./2);   %solve the ri
end
for tt=1:10000
for i=1:num_of_sensor
    r0_i1(i)=R(i)-R(1);    %generate the r0_i1
end
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,num_of_sensor);
r_i1=r0_i1+noise;              %make the noise and r_i1
for i=1:num_of_sensor
    K(i)=(sensor(1,i)).^2+sensor(2,i).^2+sensor(3,i).^2;    %make the K
end
for i=1:num_of_sensor-1
    h(i,1)=(1./2).*((r_i1(i+1)).^2-K(i+1)+K(1));  %make the h
end
for i=1:num_of_sensor-1
    Gl(i,1)=-(sensor(1,i+1)-sensor(1,1));
    Gl(i,2)=-(sensor(2,i+1)-sensor(2,1));
    Gl(i,3)=-(r_i1(i+1));          %make the Gl
end
Zl=((Gl'*(Q^(-1))*Gl)^(-1))*Gl'*(Q^(-1))*h;  %make the Zl
x1=Zl(1);
y1=Zl(2);
z1=sqrt(Zl(3)^2-x1^2-y1^2);  %solve the x1, y1 so that we will iterate them again
B=eye(num_of_sensor-1);
 for i=2:num_of_sensor
    r0(i-1)=sqrt((sensor(1,i)-x1)^2+(sensor(2,i)-y1)^2+(sensor(3,i)-z1)^2);
    B(i-1,i-1)=r0(i-1);
 end
Y=9.0*10.^16.*B*Q*B;   %Y
Zaa=((Gl'*(Y^-1)*Gl)^-1)*Gl'*(Y^-1)*h;  %zaa
x1=Zaa(1);
y1=Zaa(2);
z1=sqrt(Zaa(3)^2-x1^2-y1^2);   %solve them
EE=(x1-location(1)).^2+(y1-location(2)).^2+(z1-location(3)).^2+EE;
out=[x1;y1;z1];
end
EE=EE./10000;
end

