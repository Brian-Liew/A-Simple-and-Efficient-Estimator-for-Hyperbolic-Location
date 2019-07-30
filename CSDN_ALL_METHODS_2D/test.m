% test_Taylor_SS 本函数用于实现场强定位算法。
%               - BSN  为基站个数，3 < BSN <= 7；
%               - MSP  为移动台的初始位置, MSx, MSy均为[0,1]之间的数；
%                      特别要注意服务小区与MS之间的关系，MS的位置不能越界。
%               - ACP  为定位站的初始位置, ACx, ACy均已知，并作为移动台的初始估计位置
%               - R    为小区半径，单位(meter)；
%               - X    移动台的估计位置；
%               - zone 讨论小区环境
% See also: test_Taylor_SS.m

BSN  = 7;
%Noise = [10^2, 20^2, 30^2, 50^2, 100^2, 200^2, 300^2];
% Noise = [30, 50, 100, 120, 150, 180, 210];
Noise = [30, 60, 90, 150, 210, 300];
R = 3000;
%R = 1000;
Cal_Num = 1000;

% 加角度信息：
for i = 1: 6,
    for m = 1: Cal_Num,
%         MSP = [i*sqrt(3)/12,i/12];
%         ACP1 = [0.2,0.3];
        MSP = [0.5, 0.5];
        ACP = [0,0];
        MS = MSP*R;
        MSr(m, 1) = MS(1);
        MSr(m, 2) = MS(2);

        MSe1 = ChanAlgorithm(BSN, MSP, R, Noise(i));
        MSer1(m, 1) = MSe1(1);
        MSer1(m, 2) = MSe1(2);
        
%         MSe2 = ChanAlgorithm_A(BSN, MSP, ACP, R, Noise(i), 0.1);
%         MSer2(m, 1) = MSe2(1);
%         MSer2(m, 2) = MSe2(2);
        MSe2 = TaylorAlgorithm(BSN, MSP, R, Noise(i));
        MSer2(m, 1) = MSe2(1);
        MSer2(m, 2) = MSe2(2);
    end
    RMSE1(i) = TDOA_RMSE(MSr, MSer1);
    RMSE2(i) = TDOA_RMSE(MSr, MSer2);
%     RMSE3(i) = TDOA_RMSE(MSr, MSer3);
%     RMSE4(i) = TDOA_RMSE(MSr, MSer4);
%     RMSE5(i) = TDOA_RMSE(MSr, MSer5);
%     RMSE3(i) = TDOA_RMSE(MSr, MSer3);
end
