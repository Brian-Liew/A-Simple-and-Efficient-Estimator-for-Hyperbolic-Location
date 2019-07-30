function EE= CRLB( M, sigma_d_square)
%the function use the CRLB to solve the MSE
% M is the number of the sensors,sigma_d_square influence the Q,it output
% the MSE of the method
x(1)=0.0;  y(1)=0.0;
x(2)=-5.0; y(2)=8.0;
x(3)=4.0;  y(3)=6.0;
x(4)=-2.0; y(4)=4.0;
x(5)=7.0;  y(5)=3.0;
x(6)=-7.0; y(6)=5.0;
x(7)=2.0;  y(7)=5.0;
x(8)=-4.0; y(8)=2.0;
x(9)=3.0;  y(9)=3.0;
x(10)=1.0; y(10)=8.0;
x0=8.0;    y0=22.0;
x_=x0; y_=y0;
EE=0;
%sigma_d_square=0.001./((3.0*10.^8).^2);
T=0.5.*ones(M-1);
for tt=1:M-1
    T(tt,tt)=1;
end
Q=sigma_d_square.*T;  %the Q matrix has been solved
Ga_=[1,0;0,1;1,1];       %make the Ga_ 
for ii=1:M
    R(ii)=sqrt((x(ii)-x0).^2+(y(ii)-y0).^2);   %solve the ri
end
B=eye(M-1);
 for ii=1:M-1
    B(ii,ii)=R(ii+1);   %the B
 end
for ii=1:M
    r_i1(ii)=R(ii)-R(1);    %generate the r0_i1
end
B__=[(x0-x(1)), 0;0,(y0-y(1))];     %the B__
B_=diag([x0-x(1),y0-y(1),R(1)]);    %the B_
for ii=1:M-1
    Ga(ii,1)=-(x(ii+1)-x(1));
    Ga(ii,2)=-(y(ii+1)-y(1));
    Ga(ii,3)=-(r_i1(ii+1));          %make the Ga
end
YY=(3.0*10^8)^2.*((B__*Ga_'*(B_^-1)*Ga'*(B^-1)*(Q^-1)*(B^-1)*Ga*(B_^-1)*Ga_*B__)^-1);
EE=sum(diag(YY));


end

