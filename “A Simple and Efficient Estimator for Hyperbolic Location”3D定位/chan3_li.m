clear;
%location   sensor
location=[6,5,7];
sensor=[1,2,3,4,5,6;1,2,3,4,5,6;1,2,3,4,5,6]';
    ini=[0,0,0];
%num_of_sensors
num_of_sensor=5;
c=3e8; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   K(i) of sensor
for i = 1: num_of_sensor
    K(i) = sensor(i,1).^2 + sensor(i,2).^2+sensor(i,3).^2;
end
% K1
K1=ini(1).^2 + ini(2).^2+ini(3).^2;
%   R(i)
for i = 1: num_of_sensor
    R(i) = sqrt((sensor(i,1) - location(1))^2 + (sensor(i,2) - location(2))^2+ (sensor(i,3) - location(3))^2);
end
%R1
R1 = sqrt((location(1) - ini(1))^2+(location(2) - ini(2))^2+(location(3)-ini(3))^2);
%   Q
sigma_d_square=0.0001./((3.0*10.^8).^2);  %%%%%%
T=0.5.*ones(num_of_sensor);
for tt=1:num_of_sensor     
    T(tt,tt)=1;
end
Q=sigma_d_square.*T; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mse=0;
for pp=1:10000
%noise
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,num_of_sensor);
%Ri1(i)
for i = 1: num_of_sensor
    Ri1(i) = R(i) - R1 + noise(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h
h = 0.5*(Ri1.^2-K(1:num_of_sensor)+K1)'; 
% Ga
Ga = -[sensor(1:num_of_sensor,1)-ini(1) sensor(1:num_of_sensor,2)-ini(2)  Ri1(1:num_of_sensor)'];
% Za to get the first result
Za=inv(Ga.'*inv(Q)*Ga)*Ga.'*inv(Q)*h;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  B
for i=1:num_of_sensor
    B(i,i)=sqrt((Za(1)-ini(1))^2+(Za(2)-ini(2))^2+(sqrt(Za(3)^2-Za(1)^2-Za(2)^2)-ini(3))^2);
end
% Fa
Fa = c^2*B*Q*B;
% Za1  to get the first wls
Za1 = inv(Ga.'*inv(Fa)*Ga)*Ga.'*inv(Fa)*h;
if Za1(3)<0
  Za1(3)=abs(Za1(3));
end
out=[Za1(1);Za1(2);sqrt(Za1(3)^2-Za1(1)^2-Za1(2)^2)];
mse=mse+sqrt((out(1)-location(1)).^2+(out(2)-location(2)).^2+(out(3)-location(3)).^2);
end
mse=mse./10000;
