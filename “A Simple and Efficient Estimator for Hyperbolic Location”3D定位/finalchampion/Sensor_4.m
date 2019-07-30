function [ pp,EE ] = Sensor_4( num_of_sensor,sensor,location,sigma_d_square )
%The function  solve the tdoa of 3D when the sensor only is 4
% out is the result and mes is the MSE; the function give the number of
% sensors and the location of the sensors and the source
%clear;
%location=[6,5,7];
%sensor=[0,2,5,4,6,7;0,8,6,4,5,3;0,7,4,5,6,7];
%num_of_sensor=4;
EE=0.0;
%sigma_d_square=0.001./((3.0*10.^8).^2);
for i=1:num_of_sensor
    R(i)=sqrt((sensor(1,i)-location(1)).^2+(sensor(2,i)-location(2)).^2+(sensor(3,i)-location(3)).^2);   %solve the ri
end
for i=1:num_of_sensor
    r0_i1(i)=R(i)-R(1);    %generate the r0_i1
end
for i=1:num_of_sensor
    K(i)=(sensor(1,i)).^2+(sensor(2,i)).^2+(sensor(3,i)).^2;    %make the K
end
%for tt=1:100
    syms r1_;
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,num_of_sensor);
r_i1=r0_i1+noise;
for i=1:3
    p(i,1)=sensor(1,i+1);
    p(i,2)=sensor(2,i+1);
    p(i,3)=sensor(3,i+1);
end
p=p^-1;
Q=-p*([r_i1(2);r_i1(3);r_i1(4)].*r1_+0.5*[r_i1(2)^2-K(2)-K(1);r_i1(3)^2-K(3)-K(1);r_i1(4)^2-K(4)-K(1)]);
r1=solve(Q(1)^2+Q(2)^2+Q(3)^2-r1_^2);   %just to solve the equation
r1=double(r1);
%r1=max(r1);  %to avoid the other solution
pp=-p*([r_i1(2);r_i1(3);r_i1(4)]*r1(1)+0.5*[r_i1(2)^2-K(2)-K(1);r_i1(3)^2-K(3)-K(1);r_i1(4)^2-K(4)-K(1)]);
EE=(pp(1)-location(1)).^2+(pp(2)-location(2)).^2+(pp(3)-location(3)).^2+EE;



end

