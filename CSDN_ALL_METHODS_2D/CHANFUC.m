 % BSN表示基站的个数，STATEA的前2列是基站的横纵坐标，
   %第三列是测得的代表移动台与和服务基站（基站1）和其它各基站的距离差的TDOA值.
function F=CHANFUC(BSN,STATEA)
  
           for n=1:BSN
               K(n)=STATEA(n,1)*STATEA(n,1)+STATEA(n,2)*STATEA(n,2);
           end
           %计算X21,Y21....
           for n=1:BSN-1
               Xi1(n)=STATEA(n+1,1)-STATEA(1,1);
               Yi1(n)=STATEA(n+1,2)-STATEA(1,2);
           end    
           %计算到各基站的距离
           for n=2:BSN 
               ri1(n)=STATEA(n,3);
           end
           
           for n=1:BSN-1
               h(n,1)=0.5*(ri1(n+1)^2-K(n+1)+K(1));
           end    
           
           for n=1:BSN-1
               Ga(n,1)=-Xi1(n);
               Ga(n,2)=-Yi1(n);
               Ga(n,3)=-ri1(n+1);
           end 
         
            Q=eye(BSN-1);
                    
            Za=pinv(Ga'*pinv(Q)*Ga)*Ga'*pinv(Q)*h;
            x1=Za(1);
            y1=Za(2);
            r1=Za(3);
            
            for n=1:BSN-1
                rno(n)=sqrt((x1-STATEA(n+1,1))*(x1-STATEA(n+1,1))+(y1-STATEA(n+1,2))*(y1-STATEA(n+1,2)));
            end
            
            B=diag(rno);
           
             V=B*Q*B;
             Za1=pinv(Ga'*pinv(V)*Ga)*Ga'*pinv(V)*h;
    
            x10=Za1(1);
            y10=Za1(2);
            r10=Za1(3);
            if r10<0 
                r10=0;
            end
            %第一次加权ls定位结果
            %第二次加权定位
            covza=pinv(Ga'*pinv(V)*Ga);
             B2=diag([x10 y10 r10]);
             V2=4*B2*covza*B2;
             Ga1=[1 0;0 1;1 1];
             h1=[x10*x10;y10*y10;r10*r10];
            Zae=pinv(Ga1'*pinv(V2)*Ga1)*Ga1'*pinv(V2)*h1;
            Zp=sqrt(Zae);
            if x10<0 
             Zp1=-Zp(1);
            else
             Zp1=Zp(1);
            end
            if y10<0
              Zp2=-Zp(2);
            else
              Zp2=Zp(2);
            end%第二次加权ls定位结果
%             Zp1 = Zp(1);
%             Zp2 = Zp(2);
            F=[Zp1;Zp2];
return