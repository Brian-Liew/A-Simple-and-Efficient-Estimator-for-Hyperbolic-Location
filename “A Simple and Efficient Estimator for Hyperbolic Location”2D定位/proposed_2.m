function EE = proposed_2( M, sigma_d_square)
%the function use the proposed_2 to solve the MSE
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
EE=0;
%sigma_d_square=0.001./((3.0*10.^8).^2);
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
    Ga(ii,1)=-(x(ii+1)-x(1));
    Ga(ii,2)=-(y(ii+1)-y(1));
    Ga(ii,3)=-(r_i1(ii+1));          %make the Ga
end
Za=((Ga'*(Q^(-1))*Ga)^(-1))*Ga'*(Q^(-1))*h;  %make the Za
x1=Za(1);
y1=Za(2);   %iterate the x and y
B=eye(M-1);
 for i=2:M
    r0(i-1)=sqrt((x(1)-x1)^2+(y(1)-y1)^2);
    B(i-1,i-1)=r0(i-1);
 end
Y=9.0*10.^16.*B*Q*B;   %Y
Zaa=((Ga'*(Y^-1)*Ga)^-1)*Ga'*(Y^-1)*h;
x1=Zaa(1);
y1=Zaa(2);    %iterate the x and y again
Ga_=[1,0;0,1;1,1];       %make the Ga_ 
B_=diag([x1-x(1),y1-y(1),Za(3)]);  %make the B_
h_=[(Zaa(1)-x(1))^2;(Zaa(2)-y(1)).^2;(Zaa(3)).^2];  %h_
cov_Za=(Ga'*(Y^(-1))*Ga)^(-1);    %conv_Za
Y_=4.*B_*cov_Za*B_;             % Y_
Za_=((Ga_'*(Y_^(-1))*Ga_)^(-1))*Ga_'*(Y_^(-1))*h_; %Za_
Zp=sqrt(Za_)+[x(1);y(1)];   %make the Zp
EE=(Zp(1)-x0).^2+(Zp(2)-y0).^2+EE;
end
EE=EE./10000;


end

