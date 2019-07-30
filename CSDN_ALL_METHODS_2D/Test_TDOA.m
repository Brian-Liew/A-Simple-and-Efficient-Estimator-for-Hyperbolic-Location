%  本程序实现测试、比较无线定位算法

%  TESTALGORITHM
%    参数说明：
%       移动台初始位置由MSPosition给出
%       输出移动台的估计位置。
%       比较各个算法的RMSE、CRLB、GDOP
%       比较各个算法的时间复杂度。


BSN = [7,6,5,4];
% BSN = [7];
%BSN = 3;
R = 3000;
CalNum = 10000;
%Noise = [50^2,100^2,150^2,200^2,300^2,400^2,500^2,700^2,1000^2];
%Noise = 30^2;
%Noise = [40^2,80^2,90^2,120^2,150^2,200^2,500^2,700^2,1000^2];

% 几种算法的比较  
% for m = 1:4,
%     for n = 1:7,
%         for i = 1: 1000,
%             MSP = MSPosition();
%             MS(i ,1) = MSP(1)*R;
%             MS(i, 2) = MSP(2)*R;
%             
%             EMSMI =  EvLChanAlgorithm(BSN(m),MSP,R,Noise(n));
%             EMSM(i ,1) = EMSMI(1);
%             EMSM(i, 2) = EMSMI(2);
% 
%             EMSvCI = EvChanAlgorithm(BSN(m),MSP,R,Noise(n));
%             EMSvC(i ,1) = EMSvCI(1);
%             EMSvC(i, 2) = EMSvCI(2);
% 
%             EMSCI = ChanAlgorithm(BSN(m),MSP,R,Noise(n));
%             EMSC(i ,1) = EMSCI(1);
%             EMSC(i, 2) = EMSCI(2);
% 
%             EMSTI = TaylorAlgorithm(BSN(m),MSP,R,Noise(n));
%             EMST(i ,1) = EMSTI(1);
%             EMST(i, 2) = EMSTI(2);
%             
%             EMSvLCI = EvTaylorAlgorithm(BSN(m),MSP,R,Noise(n));
%             EMSvLC(i ,1) = EMSvLCI(1);
%             EMSvLC(i, 2) = EMSvLCI(2);
%         end
% 
%         rmse_vc(m,n) = TDOA_RMSE(MS,EMSvC);
%         rmse_vt(m,n) = TDOA_RMSE(MS,EMSvLC);
%         rmse_c(m,n) = TDOA_RMSE(MS,EMSC);
%         rmse_m(m,n) = TDOA_RMSE(MS,EMSM);
%         rmse_t(m,n) = TDOA_RMSE(MS,EMST);
%     end
% end

%几种传统定位算法的比较
% Noise = [30^2,60^2,90^2,120^2,150^2,180^2];
% Noise = [30];
Noise = [30, 60, 90, 150, 210, 300];
for m = 1:4,
    for n = 1:6,
        for i = 1: CalNum,
            MSP = MSPosition();
            MS(i,1) = MSP(1)*R;
            MS(i,2) = MSP(2)*R;

            EMSCI = ChanAlgorithm(BSN(m), MSP, R, Noise(n));
            EMSC(i, 1) = EMSCI(1);
            EMSC(i, 2) = EMSCI(2);

            EMSTI = TaylorAlgorithm(BSN(m), MSP, R, Noise(n));
            EMST(i ,1) = EMSTI(1);
            EMST(i, 2) = EMSTI(2);

            EMSSII = SIAlgorithm(BSN(m), MSP, R, Noise(n));
            EMSSI(i ,1) = EMSSII(1);
            EMSSI(i, 2) = EMSSII(2);

%             EMSSXI = SXAlgorithm(BSN(m), MSP, R, Noise(n));
%             EMSSX(i ,1) = EMSSXI(1);
%             EMSSX(i, 2) = EMSSXI(2);
        end

        rmse_c(m,n) = TDOA_RMSE(MS,EMSC);
        rmse_t(m,n) = TDOA_RMSE(MS,EMST);
        rmse_si(m,n) = TDOA_RMSE(MS,EMSSI);
        crlb(m,n) = CRLB(BSN(m), MSP, R, Noise(n));
%         rmse_sx(m,n) = TDOA_RMSE(MS,EMSSX);
    end
end

% % % 三种方法的比较：
%  Noise = [30^2,60^2,120^2,240^2,210^2,600^2,1200^2,2400^2];
% for m = 1,
%     for n = 1:8,
%         for i = 1: 10000,
%             MSP = MSPosition();
%             MS(i ,1) = MSP(1)*R;
%             MS(i, 2) = MSP(2)*R;
% 
%             EMSvCI = EvChanAlgorithm(BSN(m),MSP,R,Noise(n));
%             EMSvC(i ,1) = EMSvCI(1);
%             EMSvC(i, 2) = EMSvCI(2);
% 
%             EMSCI = ChanAlgorithm(BSN(m),MSP,R,Noise(n));
%             EMSC(i ,1) = EMSCI(1);
%             EMSC(i, 2) = EMSCI(2);
%  
%             EMSvLCI = EvTaylorAlgorithm(BSN(m),MSP,R,Noise(n));
%             EMSvLC(i ,1) = EMSvLCI(1);
%             EMSvLC(i, 2) = EMSvLCI(2);
%             
%             EMSTI = TaylorAlgorithm(BSN(m), MSP, R, Noise(n));
%             EMST(i ,1) = EMSTI(1);
%             EMST(i, 2) = EMSTI(2);           
%          end
% 
%         rmse_vc(m,n) = TDOA_RMSE(MS,EMSvC);
%         rmse_c(m,n) = TDOA_RMSE(MS,EMSC);
%         rmse_vl(m,n) = TDOA_RMSE(MS,EMSvLC);
%         rmse_t(m,n) = TDOA_RMSE(MS,EMST);
%     end
% end

% % CHAN算法、Fang算法以及CRLB的比较
%     for n = 1:5,
%         for i = 1: 10000,
%             MSP = MSPosition();
%             MS(i, 1) = R*MSP(1);
%             MS(i, 2) = R*MSP(2);
%             
%             EMS_Chan3i = Chan_3BS(MSP,R,Noise(n));
%             EMS_Chan3(i ,1) = EMS_Chan3i(1);
%             EMS_Chan3(i, 2) = EMS_Chan3i(2);
%             
%             EMS_Fangi = FangAlgorithm(MSP,R,Noise(n));
%             EMS_Fang(i ,1) = EMS_Fangi(1);
%             EMS_Fang(i, 2) = EMS_Fangi(2);
%         end
%         
%         rmse_c(n) = TDOA_RMSE(MS, EMS_Chan3);
%         rmse_f(n) = TDOA_RMSE(MS, EMS_Fang);
% %         crlb(n) = CRLB(BSN, MSP, R, Noise(n))
%     end

% % Chan算法性能分析
% Noise = [60^2,90^2,120^2,150^2,180^2,210^2,600^2];
% for m = 1:1,
%     for n = 1:7,
%         for i = 1: 1000,
%             MSP = MSPosition();
%             MS(i,1) = MSP(1)*R;
%             MS(i,2) = MSP(2)*R;
% 
%             EMSCI = ChanAlgorithm(BSN(m), MSP, R, Noise(n));
%             EMSC(i, 1) = EMSCI(1);
%             EMSC(i, 2) = EMSCI(2);
% 
%             EMSTI = TaylorAlgorithm(BSN(m), MSP, R, Noise(n));
%             EMST(i ,1) = EMSTI(1);
%             EMST(i, 2) = EMSTI(2);
%             
%         end
% 
%         rmse_c(m,n) = TDOA_RMSE(MS,EMSC);
%         rmse_t(m,n) = TDOA_RMSE(MS,EMST);
%     end
% end

% % Taylor算法分析:
% Noise = [30^2,60^2,90^2,120^2,150^2,180^2,210^2];
% for m = 1:1,
%     for n = 1:7,
%         for i = 1: 10000,
%             MSP = MSPosition();
%             MS(i,1) = MSP(1)*R;
%             MS(i,2) = MSP(2)*R;
% 
%             EMSTI_1 = TaylorAlgorithm(BSN(m), MSP, R, Noise(n));
%             EMST_1(i ,1) = EMSTI_1(1);
%             EMST_1(i, 2) = EMSTI_1(2);
%             
%             MSP_1 = MSP + [0.05 0.05];
%             EMSTI_2 = TaylorAlgorithm(BSN(m),MSP_1, R, Noise(n));
%             EMST_2(i ,1) = EMSTI_2(1);
%             EMST_2(i, 2) = EMSTI_2(2);
%             
%             MSP_2 = [sqrt(3)/4, sqrt(3)*tand(15)/4];
%             EMSTI_3 = TaylorAlgorithm(BSN(m),MSP_2, R, Noise(n));
%             EMST_3(i ,1) = EMSTI_3(1);
%             EMST_3(i, 2) = EMSTI_3(2);
%             
%             MSP_3 = MSP + [0.01 0.01];
%             EMSTI_4 = TaylorAlgorithm(BSN(m),MSP_3, R, Noise(n));
%             EMST_4(i ,1) = EMSTI_4(1);
%             EMST_4(i, 2) = EMSTI_4(2);            
%         end
% 
%         rmse_t1(m,n) = TDOA_RMSE(MS,EMST_1);
%         rmse_t2(m,n) = TDOA_RMSE(MS,EMST_2);
%         rmse_t3(m,n) = TDOA_RMSE(MS,EMST_3);
%         rmse_t4(m,n) = TDOA_RMSE(MS,EMST_4);
%     end
% end
