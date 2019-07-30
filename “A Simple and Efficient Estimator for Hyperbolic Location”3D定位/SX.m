clear;
%location   sensor
loc=[8,22];
sensor=[0,-5,4,-2,7,-7,2,-4,3,1;
        0, 8,6, 4,3, 5,5, 2,3,8];
%num_of_sensors
num_of_sensor=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   K(i) of sensor
for i = 1: num_of_sensor-1,
    K(i) = sensor(1, i+1)^2 + sensor(2, i+1)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%noise
sigma_d_square=0.001./((3.0*10.^8).^2);
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,num_of_sensor-1);
mse=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1
location=loc;
R1 = sqrt(location(1)^2 + location(2)^2);
%R(i)
for i = 1: num_of_sensor-1,
    R(i) = sqrt((sensor(1, i+1) - location(1))^2 + (sensor(2, i+1) - location(2))^2);
end
for pp=1:100000
%Ri1(i)
for i = 1: num_of_sensor-1,
    Ri1(i) = R(i) - R1 + noise(i);
end
% delt:
for i = 1: num_of_sensor-1,
    K(i) = sensor(1,i+1)^2 + sensor(2,i+1)^2;
end
for i = 1: num_of_sensor-1,
    delt(i) = K(i) - Ri1(i)^2;
end
  W = eye(num_of_sensor-1);    
% S:
for i = 1: num_of_sensor-1,
    S(i, 1) = sensor(1, i+1);
    S(i, 2) = sensor(2, i+1);
end

% Sw:
Sw = inv(S'*W*S)*S'*W;

% a:
a = 4 - 4*Ri1*Sw'*Sw*Ri1';

% b:
b = 4*Ri1*Sw'*Sw*delt';

% c:
c = -delt*Sw'*Sw*delt';

% root:
root1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
root2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
out=[abs(root2),abs(root1)];
mse=mse+(out(1)-loc(1)).^2+(out(2)-loc(2)).^2;
end
mse=mse./100000;