function EE = Three_sensor_linear(  )
%the function use the three sensors to solve the near  linear sources 
%   M is the number of the sensors,sigma_d_square influence the Q,it output
%   the MSE of the method
x(1)=0.0;  y(1)=0.0;
x(2)=2.0; y(2)=0.0;
x(3)=-2.0;  y(3)=0.0; 
X=8.0;Y=22.0;
M=3;
EE=0.0;
L1=2;L2=2;
sigma_d_square=0.0001./((3.0*10.^8).^2);
for ii=1:M
    R(ii)=sqrt((x(ii)-X).^2+(y(ii)-Y).^2);   %solve the ri
end
for ii=1:M
    r0_i1(ii)=R(ii)-R(1);    %generate the r0_i1
end
for tt=1:10000
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,M);
r_i1=r0_i1+noise;
rr=(L1.*(1-(r_i1(2)./L1).^2)+L2*(1-(r_i1(3)./L2).^2))./(2*(r_i1(3)./L2+r_i1(2)./L1));   %the r
x1=-(L2^2*r_i1(2)-L1^2*r_i1(3)-r_i1(2)*r_i1(3)*(r_i1(3)-r_i1(2)))./(2*(r_i1(2)*L2+L1*r_i1(3)));%the x1
y1=sqrt(rr^2-x1^2); %the y
EE=(x1-X).^2+(y1-Y).^2+EE;
end
EE=(EE)./10000;

end

