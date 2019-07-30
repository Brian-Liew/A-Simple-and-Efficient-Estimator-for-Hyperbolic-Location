clear;
location=[6,5,7];
sensor=[0,2,5,4,6,7,8,9;0,8,6,4,12,5,3,9;0,7,4,5,6,14,7,14];
sensor1=[2,5,4,6,7,15,8;8,6,4,5,16,3,9;7,4,14,5,6,7,14];
sigma_d_square=0.0001./((3.0*10.^8).^2);
fprintf('                     Table 1                                          ');
fprintf('\n');
fprintf('---------------------------------------------------------');
fprintf('\n');
fprintf('MSE           M=4       M=5       M=6       M=7       M=8   ');
fprintf('\n');
fprintf('--------------------------------------------------------');
fprintf('\n');
fprintf('SI ');
fprintf('               ');
 for num_=5:6
    [a,b]=SI_3D(num_,sensor,location,sigma_d_square);
    fprintf('%10.5f',b);
 end
 
 fprintf('\n');fprintf('Chan_3D ');
 [a,b]=Sensor_4(4,sensor,location,sigma_d_square);
fprintf('%10.5f',b);
  for num_=5:8
    [a,b]=Chan_3D(num_-1,sensor1',location,sigma_d_square,[0;0;0]');
    fprintf('%10.5f',b);
  end
  
  fprintf('\n');fprintf('Taylor  ');
  [a,b]=Sensor_4(4,sensor,location,sigma_d_square);
fprintf('%10.5f',b);
  for num_=5:8
    [a,b]=Taylor_3D(num_,sensor,location,sigma_d_square);
    fprintf('%10.5f',b);
  end
  fprintf('\n');fprintf('CRLB    ');
  for num_=4:8
    [b]=CRLB_3D(num_,sensor,location,sigma_d_square*10);
    fprintf('%10.5f',b);
  end
 fprintf('\n');     %set the form
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 clear;
location=[118,222,40];
sensor=[0 0 0;120 106 0;200 398 0;46 320 0;340 354 0;136 420 0;86 35 0;38 465 0;208 147 0;143 369 0]';
sensor1=[120 106 0;200 398 0;46 320 0;340 354 0;136 420 0;86 35 0;38 465 0;208 147 0;143 369 0]';
sigma_d_square=0.0001./((3.0*10.^8).^2);
fprintf('                     Table 2 (linear)                            ');
fprintf('\n');
fprintf('---------------------------------------------------------');
fprintf('\n');
fprintf('MSE           M=4       M=5       M=6       M=7       M=8   ');
fprintf('\n');
fprintf('--------------------------------------------------------');
 fprintf('\n');fprintf('Chan_3D ');
 [a,b]=Sensor_4_linear(4,sensor,location,sigma_d_square);
fprintf('%10.5f',b);
  for num_=5:8
    [a,b]=Chan_3D_linear(num_,sensor,location,sigma_d_square);
    fprintf('%10.5f',b);
  end
  
  fprintf('\n');fprintf('Taylor  ');
  [a,b]=Sensor_4_linear(4,sensor,location,sigma_d_square);
fprintf('%10.5f',b);
  for num_=5:8
    [a,b]=Taylor_3D(num_,sensor,location,sigma_d_square);
    fprintf('%10.5f',b);
  end
  fprintf('\n');fprintf('CRLB    ');
  for num_=4:8
    [b]=CRLB_3D_linear(num_,sensor,location,sigma_d_square*10);
    fprintf('%10.5f',b);
  end
 fprintf('\n');     %set the form