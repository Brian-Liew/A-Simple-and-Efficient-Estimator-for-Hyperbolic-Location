function [ EE ] = CRLB_3D_linear( num_of_sensor,sensor,location,sigma_d_square  )
%The function use the CRLB to solve the mse of tdoa of 3D
% mes is the MSE; the function give the number of
% sensors and the location of the sensors and the source
%clear;
%clear;
%location   sensor
%location=[40,50,30];
%sensor=[0,10,20,30,40,50;0,10,20,30,40,50;0,10,20,30,40,50];
%num_of_sensors
%num_of_sensor=5;
EE=0;
%sigma_d_square=0.001./((3.0*10.^8).^2);
T=0.5.*ones(num_of_sensor-1);
for tt=1:num_of_sensor-1
    T(tt,tt)=1;
end
Q=sigma_d_square.*T;  %the Q matrix has been solved
Ga_=[1,0,0;0,1,0;0,0,1;1,1,1];       %make the Ga_ 
for i=1:num_of_sensor
    R(i)=sqrt((sensor(1,i)-location(1)).^2+(sensor(2,i)-location(2)).^2+(sensor(3,i)-location(3)).^2);   %solve the ri
end
B=eye(num_of_sensor-1);
 for i=1:num_of_sensor-1
    B(i,i)=R(i+1);   %the B
 end
for i=1:num_of_sensor
    r_i1(i)=R(i)-R(1);    %generate the r0_i1
end
B_=[1,0,0;0,1,0;location(1)/R(1),location(2)/R(1),location(3)/R(1)];
for i=1:num_of_sensor-1
    Ga(i,1)=-(sensor(1,i+1)-location(1));
    Ga(i,2)=-(sensor(2,i+1)-location(2));
    Ga(i,3)=-(r_i1(i+1));          %make the Ga
end
YY=(3.0*10^8)^2.*((B_'*Ga'*(B^-1)*(Q^-1)*(B^-1)*Ga*B_)^-1);
EE=sum(diag(YY));

end

