% num_of_sensor,sensor,location,sigma_d_square
clear;
num_of_sensor=6;
location=[6,5,7];
sensor=[0,-15,5,-10,10,2;0,20,2,16,10,3;0,10,20,9,-15,7];
sigma_d_square=0.001./((3.0*10.^8).^2);

figure(1)
plot3(sensor(1,:),sensor(2,:),sensor(3,:),'rp',6,5,7,'g*');
axis([-30 30 -30 30 -30 30]);
box;
text(6,5,7,'信号源');grid on;
xlabel('x轴');ylabel('y轴');zlabel('z轴');hold on;
out=SI_3D(num_of_sensor,sensor,location,sigma_d_square);
plot3(out(1),out(2),out(3),'ms');text(out(1),out(2),out(3),'估计源');grid on;
legend('sensor','source','estimated point');
title('SI')
figure(2)
plot3(sensor(1,:),sensor(2,:),sensor(3,:),'rp',6,5,7,'g*');
axis([-30 30 -30 30 -30 30]);
box;
text(6,5,7,'信号源');grid on;
xlabel('x轴');ylabel('y轴');zlabel('z轴');hold on;
out=Taylor_3D(num_of_sensor,sensor,location,sigma_d_square);
plot3(out(1),out(2),out(3),'ms');text(out(1),out(2),out(3),'估计源');grid on;
legend('sensor','source','estimated point');
title('Taylor')

figure(3)
plot3(sensor(1,:),sensor(2,:),sensor(3,:),'rp',6,5,7,'g*');
axis([-30 30 -30 30 -30 30]);
box;
text(6,5,7,'信号源');grid on;
xlabel('x轴');ylabel('y轴');zlabel('z轴');hold on;
out=Sensor_4(num_of_sensor,sensor,location,sigma_d_square);
plot3(out(1),out(2),out(3),'ms');text(out(1),out(2),out(3),'估计源');grid on;
legend('sensor','source','estimated point');
title('4_sensor')

figure(4)
sensor=[-15,5,-10,10,2;20,2,16,10,3;10,20,9,-15,7];
plot3(sensor(1,:),sensor(2,:),sensor(3,:),'rp',6,5,7,'g*');
axis([-30 30 -30 30 -30 30]);
box;
text(6,5,7,'信号源');grid on;
xlabel('x轴');ylabel('y轴');zlabel('z轴');hold on;
out=Chan_3D(num_of_sensor-1,sensor',location,sigma_d_square,[0;0;0]');
plot3(out(1),out(2),out(3),'ms');text(out(1),out(2),out(3),'估计源');grid on;
legend('sensor','source','estimated point');
title('Chan')