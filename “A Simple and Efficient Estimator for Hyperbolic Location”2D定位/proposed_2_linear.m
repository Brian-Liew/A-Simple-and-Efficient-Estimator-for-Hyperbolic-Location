function EE = proposed_2_linear( M,sigma_d_square )
%the function use the proposed_2 to solve the near  linear sources 
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
EE=0;
%sigma_d_square=0.0001./((3.0*10.^8).^2);
T=0.5.*ones(M-1);
for tt=1:M-1
    T(tt,tt)=1;
end
Q=sigma_d_square.*T;      %the Q matrix has been solved
for ii=1:M
    R(ii)=((x(ii)-x0).^2+(y(ii)-y0).^2).^(1./2);   %solve the ri
end
for tt=1:10000
for ii=1:M
    r0_i1(ii)=R(ii)-R(1);    %generate the r0_i1
end
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,M);
r_i1=r0_i1+noise;              %make the noise and r_i1
for ii=1:M
    K(ii)=(x(ii)).^2+(y(ii)).^2;    %make the K
end
for ii=1:M-1
    h(ii,1)=(1./2).*((r_i1(ii+1)).^2-K(ii+1)+K(1));  %make the h
end
for ii=1:M-1
    Gl(ii,1)=-(x(ii+1)-x(1));
    Gl(ii,2)=-(r_i1(ii+1));          %make the Gl
end
Zl=((Gl'*(Q^(-1))*Gl)^(-1))*Gl'*(Q^(-1))*h;  %make the Zl
x1=Zl(1);
y1=sqrt(Zl(2)^2-x1^2);  %solve the x1, y1 so that we will iterate them again
B=eye(M-1);
 for i=2:M
    r0(i-1)=sqrt((x(1)-x1)^2+(y(1)-y1)^2);
    B(i-1,i-1)=r0(i-1);
 end
Y=9.0*10.^16.*B*Q*B;   %Y
Zaa=((Gl'*(Y^-1)*Gl)^-1)*Gl'*(Y^-1)*h;  %zaa
x1=Zaa(1);
y1=sqrt(Zaa(2)^2-x1^2);    %solve them
EE=(x1-x0).^2+(y1-y0).^2+EE;
end
EE=EE./10000;

end

