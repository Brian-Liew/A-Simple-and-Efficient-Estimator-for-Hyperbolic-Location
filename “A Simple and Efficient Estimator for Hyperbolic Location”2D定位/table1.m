s1=0.001./((3.0*10.^8).^2);
s2=0.0001./((3.0*10.^8).^2);
s3=0.00001./((3.0*10.^8).^2);
fprintf('                                  Table 1                                          ');
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------');
fprintf('\n');
fprintf('MSE      M=3      M=4       M=5       M=6       M=7       M=8       M=9       M=10   ');
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------');
fprintf('\n');
fprintf('B  ');
fprintf('%10.5f',Three_sensor_linear());
 for num_=4:10
    fprintf('%10.5f',Taylor_series(num_,s1));
 end
 fprintf('\n');fprintf('C  ');
 fprintf('%10.5f',Three_sensor_linear());
  for num_=4:10
    fprintf('%10.5f',proposed_1(num_,s1));
  end
  fprintf('\n');fprintf('D  ');
  fprintf('%10.5f',Three_sensor_linear());
  for num_=4:10
    fprintf('%10.5f',proposed_2(num_,s1));
  end
  fprintf('\n');fprintf('E  ');
  for num_=3:10
    fprintf('%10.5f',CRLB(num_,s1));
  end
 fprintf('\n');     %set the form