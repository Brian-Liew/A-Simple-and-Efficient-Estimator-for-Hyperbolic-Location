function [ out,EE ] = Taylor_3D( num_of_sensor,sensor,location,sigma_d_square )
%The function use the Taylor to solve the tdoa of 3D
% out is the result and mes is the MSE; the function give the number of
% sensors and the location of the sensors and the source
%clear;
%location=[6,5,7];         %the sensors position
loc=location;
%out=loc;
%num_of_sensor=5;
%sensor=[0,2,5,4,6,7;0,8,6,4,5,3;0,7,4,5,6,7];
EE=0;         %initiate the MSE
%sigma_d_square=0.0001./((3.0*10.^8).^2);
T=0.5.*ones(-1);
for tt=1:num_of_sensor-1
    T(tt,tt)=1;
end
Q=sigma_d_square.*T;      %the Q matrix has been solved
for i=1:num_of_sensor
    R(i)=sqrt((sensor(1,i)-location(1)).^2+(sensor(2,i)-location(2)).^2+(sensor(3,i)-location(3)).^2);   %solve the ri
end
for pp=1:10000        %to iterate the process
for i=1:num_of_sensor
    r0_i1(i)=R(i)-R(1);    %generate the r0_i1
end
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,num_of_sensor);
r_i1=r0_i1+noise;              %make the noise and r_i1
nn=0;
deta=[1;1:1];        %the devitation of the x and y
while(norm(deta)>exp(-15)&&nn<10000)    %when the norm is small enough or the times is large enough
for i=1:num_of_sensor
    r(i)=sqrt((sensor(1,i)-location(1)).^2+(sensor(2,i)-location(2)).^2+(sensor(3,i)-location(3)).^2);   %solve the ri
end          
for i=1:num_of_sensor-1
    ht(i,1)=r_i1(i+1)-(r(i+1)-r(1));  %solve the ht
end
for i=1:num_of_sensor-1
    Gt(i,1)=(sensor(1,1)-location(1))./r(1)-(sensor(1,i+1)-location(1))./r(i+1);
    Gt(i,2)=(sensor(2,1)-location(2))./r(1)-(sensor(2,i+1)-location(2))./r(i+1);  %get the Gt
    Gt(i,3)=(sensor(3,1)-location(3))./r(1)-(sensor(3,i+1)-location(3))./r(i+1);
end
deta=(((Gt'*(Q^(-1))*Gt)^(-1))*Gt'*(Q^(-1))*ht);  %solve the deta
location(1)=location(1)+deta(1);
location(2)=location(2)+deta(2);    %iterate the x, y
location(3)=location(3)+deta(3);
nn=nn+1;
end
out=location;
EE=(location(1)-loc(1))^2+(location(2)-loc(2))^2+(location(3)-loc(3))^2+EE;   %to calculate the MSE
end
EE=EE./10000;

end

