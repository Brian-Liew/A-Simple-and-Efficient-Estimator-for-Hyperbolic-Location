function EE = Taylor_series_linear( M, sigma_d_square)
%the function use the Taylor-series to solve the near  linear sources 
%   M is the number of the sensors,sigma_d_square influence the Q,it output
%   the MSE of the method
x(1)=0.0;  y(1)=0.0;
x(2)=2.0; y(2)=0.0;
x(3)=-2.0;  y(3)=0.0;
x(4)=4.0; y(4)=0.0;
x(5)=-4.0;  y(5)=0.0;
x(6)=6.0; y(6)=0.0;
x(7)=-6.0;  y(7)=0.0;
x(8)=8.0; y(8)=0.0;
x(9)=-8.0;  y(9)=0.0;
x(10)=10.0; y(10)=0.0;
x0=8.0;    y0=22.0;
x_=x0;y_=y0;
EE=0;
%sigma_d_square=0.001./((3.0*10.^8).^2);
T=0.5.*ones(M-1);
for tt=1:M-1
    T(tt,tt)=1;
end
Q=sigma_d_square.*T;      %the Q matrix has been solved
for ii=1:M
    R(ii)=sqrt((x(ii)-x_).^2+(y(ii)-y_).^2);   %solve the ri
end
for pp=1:10000
for ii=1:M
    r0_i1(ii)=R(ii)-R(1);    %generate the r0_i1
end
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,M);
r_i1=r0_i1+noise;              %make the noise and r_i1
nn=0;
deta=[1;1];
while(norm(deta)>exp(-15)&&nn<10000)
for ii=1:M
    r(ii)=sqrt((x(ii)-x_).^2+(y(ii)-y_).^2);   %solve the ri
end          
for ii=1:M-1
    ht(ii,1)=r_i1(ii+1)-(r(ii+1)-r(1));  %solve the ht
end
for ii=1:M-1
    Gt(ii,1)=(x(1)-x_)./r(1)-(x(ii+1)-x_)./r(ii+1);
    Gt(ii,2)=(y(1)-y_)./r(1)-(y(ii+1)-y_)./r(ii+1);  %get the Gt
end
deta=(((Gt'*(Q^(-1))*Gt)^(-1))*Gt'*(Q^(-1))*ht);  %solve the deta
x_=x_+deta(1);
y_=y_+deta(2);
nn=nn+1;
end
EE=(x_-x0)^2+(y_-y0)^2+EE;
end
EE=EE./10000;

end

