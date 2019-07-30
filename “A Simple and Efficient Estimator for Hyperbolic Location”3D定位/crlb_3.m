clear;
%location   sensor
location=[6,5,7];
sensor=[0,2,5,4,6,7;0,8,6,4,5,3;0,7,4,5,6,7];
%num_of_sensors
num_of_sensor=4;
EE=0;
sigma_d_square=0.001./((3.0*10.^8).^2);
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
B__=[(location(1)-sensor(1,1)), 0,0;0,(location(2)-sensor(2,1)),0;0,0,(location(3)-sensor(3,1))];     %the B__
B_=diag([location(1)-sensor(1,1),location(2)-sensor(2,1),location(3)-sensor(3,1),R(1)]);    %the B_
for i=1:num_of_sensor-1
    Ga(i,1)=-(sensor(1,i+1)-location(1));
    Ga(i,2)=-(sensor(2,i+1)-location(2));
    Ga(i,3)=-(sensor(3,i+1)-location(3));
    Ga(i,4)=-(r_i1(i+1));          %make the Ga
end
YY=(3.0*10^8)^2.*((B__*Ga_'*(B_^-1)*Ga'*(B^-1)*(Q^-1)*(B^-1)*Ga*(B_^-1)*Ga_*B__)^-1);
EE=sum(diag(YY));