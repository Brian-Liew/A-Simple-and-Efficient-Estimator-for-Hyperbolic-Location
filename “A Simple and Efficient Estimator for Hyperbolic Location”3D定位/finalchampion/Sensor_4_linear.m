function [ out,EE ] = Sensor_4_linear( num_of_sensor,sensor,location,sigma_d_square )
%The function  solve the tdoa of 3D when the sensor only is 4
% out is the result and mes is the MSE; the function give the number of
% sensors and the location of the sensors and the source
%clear;
%location=[6,5,7];
%sensor=[0,1,2,3,4,5,6;0,1,2,3,4,5,6;0,1,2,3,4,5,6];
%num_of_sensor=4;
EE=0.0;
%sigma_d_square=0.0001./((3.0*10.^8).^2);
for i=1:num_of_sensor
    R(i)=sqrt((sensor(1,i)-location(1)).^2+(sensor(2,i)-location(2)).^2+(sensor(3,i)-location(3)).^2);   %solve the ri
end
for i=1:num_of_sensor
    r0_i1(i)=R(i)-R(1);    %generate the r0_i1
end
for i=1:num_of_sensor
    K(i)=(sensor(1,i)).^2+(sensor(2,i)).^2+(sensor(3,i)).^2;    %make the K
end
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,num_of_sensor);
r_i1=r0_i1+noise;
for i=1:num_of_sensor-1
    X_c(i,1)=sensor(1,i+1)-sensor(1,1);
    Y_c(i,1)=sensor(2,i+1)-sensor(2,1);
    Z_c(i,1)=sensor(3,i+1)-sensor(3,1);
end
B=[r_i1(2)^2-K(2)+K(1);r_i1(3)^2-K(3)+K(1);r_i1(4)^2-K(4)+K(1)];
syms r1;
syms x y z;
[x,y,r1]=solve('2*X_c(1)*x+2*Y_c(1)*y+2*r_i1(2)*r1+B(1)=0','2*X_c(2)*x+2*Y_c(2)*y+2*r_i1(3)*r1+B(2)=0','2*X_c(3)*x+2*Y_c(3)*y+2*r_i1(4)*r1+B(3)=0',x,y,r1);
f=x^2+y^2+z^2-r1^2;
z=solve(f,'z');
zz=eval(z);zz=abs(zz(1));
out=[eval(x);eval(y);zz];
EE=(out(1)-location(1)).^2+(out(2)-location(2)).^2+(out(3)-location(3)).^2+EE;

end

