function EE = Taylor_three( )
%the function is used to solve the three sensors 
% to solve the MSE
x(1)=0.0;  y(1)=0.0;
x(2)=-5.0; y(2)=8.0;
x(3)=4.0;  y(3)=6.0; 
X=8.0;Y=22.0;
M=3;
EE=0.0;
sigma_d_square=0.001./((3.0*10.^8).^2);
for ii=1:M
    R(ii)=sqrt((x(ii)-X).^2+(y(ii)-Y).^2);   %solve the ri
end
for ii=1:M
    r0_i1(ii)=R(ii)-R(1);    %generate the r0_i1
end
for ii=1:M
    K(ii)=(x(ii)).^2+(y(ii)).^2;    %make the K
end
for tt=1:1000
    syms r1_;
noise=3*10.^8*normrnd(0,sqrt(sigma_d_square./2),1,M);
r_i1=r0_i1+noise;
Q=-[(x(2)-x(1)),(y(2)-y(1));(x(3)-x(1)),(y(3)-y(1))]^-1*([r_i1(2);r_i1(3)].*r1_+0.5*[r_i1(2)^2-K(2)-K(1);r_i1(3)^2-K(3)-K(1)]);
r1=solve(Q(1)^2+Q(2)^2-r1_^2);   %just to solve the equation
r1=double(r1);
r1=max(r1);  %to avoid the other solution
pp=-[(x(2)-x(1)),(y(2)-y(1));(x(3)-x(1)),(y(3)-y(1))]^-1*([r_i1(2);r_i1(3)]*r1(1)+0.5*[r_i1(2)^2-K(2)-K(1);r_i1(3)^2-K(3)-K(1)]);
EE=(pp(1)-X).^2+(pp(2)-Y).^2+EE;
end
EE=(EE)./1000;

end

