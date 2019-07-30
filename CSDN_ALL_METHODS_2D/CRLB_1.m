function X = CRLB( BSN, MSP, R, Noise )
%   Cramer-Rao Lower Bound ����ƫ���Ƶ������½�
%   CRML����˵��:
%       BSN: ��վ�ĸ�����
%       MS:  �ƶ�̨��λ�ã�����MSx/MSy����[0��1]��
%       R:   С���뾶��
%
%       X:   CRLB��һ�������ʽ�����Խ����ϵ�һ��Ԫ�ر�ʾ�������CRLB��
%            �ڶ���Ԫ�ر�ʾ��λ���ȵ�CRLB��
%   See also CRLB.m


%   ������飺
if  nargout>1,
    error('Too many output arguments.');
end
if nargin ~= 4,
    error('Wrong number of input arguments.');
end
% if BSN <= 3,
%     error('The number of BSs must be larger than 3 for this program.');
% end
flag = size(MSP);
if flag(1)~=1 | flag(2)~=2,
    error('Wrong position vector!');
end

%   ��ʼ��������:
    BS  = NetworkTop(BSN);
    c   = 3*10^8;
    Dev = Noise^2/c^2;
    Q  = 0.5*Dev*(eye(BSN-1) + ones(BSN-1));
%     Q = eye(BSN-1);
    %Q = Dev*eye(BSN-1);
%   �㷨���̣�
    MS = R*MSP;
    % MS��BS֮�����ʵ����
    for i = 1: BSN,
        BSR(i) = sqrt((BS(1, i) - MS(1))^2 + (BS(2, i ) - MS(2))^2);
    end
    
    %Gt
    for i = 1: BSN-1,
        Gt(i, 1) = (BS(1, 1) - MS(1))/BSR(1) - (BS(1, i+1) - MS(1))/BSR(i+1);
        Gt(i, 2) = (BS(2, 1) - MS(2))/BSR(1) - (BS(2, i+1) - MS(2))/BSR(i+1);
    end
    
    mCRLB = c^2*inv(Gt'*inv(Q)*Gt);
    Crlb = sqrt(mCRLB(1,1) + mCRLB(2,2));
    
%     BSR_1 = sqrt(MS(1)^2 + MS(2)^2); % R(1)
%     
%     % R(2) --- R(BSN)
%     for i = 1: BSN-1,
%         BSR(i) = sqrt((BS(1, i+1) - MS(1))^2 + (BS(2, i+1) - MS(2))^2);
%     end
%     % Ga0
%     for i = 1:BSN -1,
%         Ga(i,1) = -BS(1, i+1);
%         Ga(i,2) = -BS(2, i+1);
%         Ga(i,3) = -BSR(i);
%     end
%     % Ga'
%     mGa = [1, 0; 0, 1; 1, 1];
%     % B
%     B = zeros(BSN-1, BSN-1);
%     for i = 1: BSN-1,
%         B(i, i) = BSR(i);
%     end
%     % B'
%     mB = [MS(1), 0, 0; 0, MS(2), 0; 0, 0, BSR_1];
%     % B''
%     mmB = [MS(1), 0; 0, MS(2);];
%     
% %   ���
%     mCrlb = c*c*inv(mmB*mGa'*inv(mB)*Ga'*inv(B)*inv(Q)*inv(B)*Ga*inv(mB)*mGa*mmB);
%     Crlb = sqrt(mCrlb(1,1) + mCrlb(2,2));
    
    if nargout == 1, 
        X = Crlb;
    elseif nargout == 0,
        disp( Crlb );
    end;