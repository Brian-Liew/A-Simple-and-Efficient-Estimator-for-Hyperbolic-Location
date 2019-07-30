function EE = CRLB_distant_linear(  M,sigma_d_square  )
%the function use the CRLB to solve the distant linear sources 
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
x0=-50.0;    y0=250.0;
x_=x0; y_=y0;
EE=0;
%sigma_d_square=0.0001./((3.0*10.^8).^2);
T=0.5.*ones(M-1);
for tt=1:M-1
    T(tt,tt)=1;
end
Q=sigma_d_square.*T;  %the Q matrix has been solved
for ii=1:M
    R(ii)=sqrt((x(ii)-x0).^2+(y(ii)-y0).^2);   %solve the ri
end
B=eye(M-1);
 for ii=1:M-1
    B(ii,ii)=R(ii+1);
 end
for ii=1:M
    r_i1(ii)=R(ii)-R(1);    %generate the r0_i1
end
for ii=1:M-1
    Gl(ii,1)=-(x(ii+1)-x(1));
    Gl(ii,2)=-(r_i1(ii+1));          %make the Gl
end
T_=[1,0;(x0-x(1))/R(1),(y0-y(1))/R(1)];  %make the T_
YY=(3.0*10^8)^2.*((T_'*Gl'*(B^-1)*(Q^-1)*(B^-1)*Gl*T_)^-1);
EE=sum(diag(YY));
end

