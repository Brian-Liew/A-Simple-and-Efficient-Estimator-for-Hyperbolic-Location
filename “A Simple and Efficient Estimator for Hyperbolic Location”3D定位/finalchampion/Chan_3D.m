function [ out,mse ] = Chan_3D( num_of_sensor,sensor,location,sigma_d_square,ini)
%The function is to use the chan (proposed_2) to solve the 3D tdoa
%   out is the result and mes is the MSE; the function give the number of
%   sensors and the location of the sensors and the source while ini is the
%   initial location we give
%   sensor
%location   sensor
%location=[6,5,7];
%sensor=[2,5,4,6,7;8,6,4,5,3;7,4,5,6,7]';
%ini=[0,0,0];
%num_of_sensors
%num_of_sensor=4;
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
%sigma_d_square=0.0001./((3.0*10.^8).^2);  %%%%%%
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
Ga = -[sensor(1:num_of_sensor,1)-ini(1) sensor(1:num_of_sensor,2)-ini(2) sensor(1:num_of_sensor,3)-ini(3) Ri1(1:num_of_sensor)'];
% Za to get the first result
Za=inv(Ga.'*inv(Q)*Ga)*Ga.'*inv(Q)*h;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  B
for i=1:num_of_sensor
    B(i,i)=sqrt((Za(1)-ini(1))^2+(Za(2)-ini(2))^2+(Za(3)-ini(3))^2);
end
% Fa
Fa = c^2*B*Q*B;
% Za1  to get the first wls
Za1 = inv(Ga.'*inv(Fa)*Ga)*Ga.'*inv(Fa)*h;
if Za1(4)<0
   Za1(4)=abs(Za1(4));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cov
cov = inv(Ga.'*inv(Fa)*Ga);
% h2
h2=[(Za1(1)-ini(1))^2;(Za1(2)-ini(2))^2;(Za1(3)-ini(3))^2;Za1(4)^2];
% Ga2
Ga2=[1 0 0;0 1 0;0 0 1;1 1 1];
% B2
B2=diag([Za1(1),Za1(2),Za1(3),Za1(4)]);
% Fa2
Fa2=4*B2*cov*B2;
% Za2
Za2=inv(Ga2.'*inv(Fa2)*Ga2)*Ga2.'*inv(Fa2)*h2;
Za2=abs(Za2);
out=diag(sign(Za1(1:3)-ini'))*sqrt(Za2)+ini';
mse=mse+sqrt((out(1)-location(1)).^2+(out(2)-location(2)).^2+(out(3)-location(3)).^2);
end
mse=mse./10000;

end

