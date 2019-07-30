clear;
bj=-30:1:-12;
sigm=((10.^(bj/10))./(3*10^8)).^2;
for  mp=1:19
mse_1(mp)=10.*log(CRLB_distant(9,sigm(mp)));
end
plot(bj,mse_1,'r--');
%scatter(bj,mse_1,'r');
hold on;
for  mp=1:19
mse_1(mp)=10.*log(proposed_1_distant(9,sigm(mp)));
end
plot(bj,mse_1,'g');
%scatter(bj,mse_1,'g');
hold on;
for  mp=1:19
mse_1(mp)=10.*log(Taylor_series_distant_(9,sigm(mp)));
end
plot(bj,mse_1,'k');
%scatter(bj,mse_1,'k');
hold on;for  mp=1:19
mse_1(mp)=10.*log(proposed_2_distant(9,sigm(mp)));
end
plot(bj,mse_1,'y');
%scatter(bj,mse_1,'y');
hold on;
legend('CRLB','proposed_1','Taylor','proposed_2');
xlabel('10log(cd)');
ylabel('10log(MSE)');
title('Comparison of 4 methods between CRLB and new method and Taylor');
