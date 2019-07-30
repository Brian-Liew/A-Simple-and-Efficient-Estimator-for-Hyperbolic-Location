function [ out, mse ] = SI_3D( num_of_sensor,sensor,loc,sigma_d_square )
%The function use the SI to solve the tdoa of 3D
% out is the result and mes is the MSE; the function give the number of
% sensors and the location of the sensors and the source
%clear;
%location   sensor
%loc=[6,5,7];
%sensor=[0,2,5,4,6,7;0,8,6,4,5,3;0,7,4,5,6,7];
%num_of_sensors
%num_of_sensor=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   K(i) of sensor
for i = 1: num_of_sensor-1,
    K(i) = sensor(1, i+1)^2 + sensor(2, i+1)^2+ sensor(3, i+1)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=0;
mse=0;
for pp=1:10000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%noise
%sigma_d_square=0.001./((3.0*10.^8).^2);
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,num_of_sensor-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1
location=loc;
R1 = sqrt(location(1)^2 + location(2)^2+ location(3)^2);
%R(i)
nn=0;
for i = 1: num_of_sensor-1,
    R(i) = sqrt((sensor(1, i+1) - location(1))^2 + (sensor(2, i+1) - location(2))^2+ (sensor(3, i+1) - location(3))^2);
end
%Ri1(i)
for i = 1: num_of_sensor-1,
    Ri1(i) = R(i) - R1 + noise(i);
end
for i = 1: num_of_sensor-1,
    delt(i) = K(i) - Ri1(i)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pd orthognol
eye_m = eye(num_of_sensor-1);
coef = Ri1*Ri1';
Pd_o = eye_m - (Ri1'*Ri1/coef);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W
W = eye(num_of_sensor-1);  
%%%%%%%%%%%%%%%%
% S
for i = 1: num_of_sensor-1,
    S(i, 1) = sensor(1, i+1);
    S(i, 2) = sensor(2, i+1);
    S(i, 3) = sensor(3, i+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Za
Za = 0.5*inv(S'*Pd_o*W*Pd_o*S)*S'*Pd_o*W*Pd_o*delt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
location(1)=Za(1);
location(2)=Za(2);
location(3)=Za(3);
out=Za;
 mse=mse+(location(1)-loc(1)).^2+(location(2)-loc(2)).^2+(location(3)-loc(3)).^2;
 n=n+1;
end
mse=mse./n;

end

