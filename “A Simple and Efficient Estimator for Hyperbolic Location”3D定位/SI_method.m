function mse = SI_method( num_of_sensor,sensor,location,sigma_d_square )
%The function is the SI method to solve the mse
%   sensor is the location of sensors and the location is the initial
%   location while the num_of_sensor is the number of the sensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   K(i) of sensor
for i = 1: num_of_sensor-1,
    K(i) = sensor(1, i+1)^2 + sensor(2, i+1)^2;
end
%  R1
R1 = sqrt(location(1)^2 + location(2)^2);
%R(i)
for i = 1: num_of_sensor-1,
    R(i) = sqrt((sensor(1, i+1) - location(1))^2 + (sensor(2, i+1) - location(2))^2);
end
% W
W = eye(num_of_sensor-1);  
% S
for i = 1: num_of_sensor-1,
    S(i, 1) = sensor(1, i+1);
    S(i, 2) = sensor(2, i+1);
end
loc=location;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=0;
mse=0;
for pp=1:100000
    noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,num_of_sensor-1);
    %Ri1(i)
    for i = 1: num_of_sensor-1
        Ri1(i) = R(i) - R1 + noise(i);
    end
    %delt
    for i = 1: num_of_sensor-1
        delt(i) = K(i) - Ri1(i)^2;
    end
    % Pd orthognol
    eye_m = eye(num_of_sensor-1);
    coef = Ri1*Ri1';
    Pd_o = eye_m - (Ri1'*Ri1/coef);
    % Za
    Za = 0.5*inv(S'*Pd_o*W*Pd_o*S)*S'*Pd_o*W*Pd_o*delt';
    location(1)=Za(1);
    location(2)=Za(2);
    mse=mse+(location(1)-loc(1)).^2+(location(2)-loc(2)).^2;
    n=n+1;
end
   mse=mse./n; 
end

