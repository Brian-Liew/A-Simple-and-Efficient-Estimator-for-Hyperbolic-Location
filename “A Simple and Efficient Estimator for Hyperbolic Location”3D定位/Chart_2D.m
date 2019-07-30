location=[6,5,7];
sensor=[0,-15,5,-10,10,2;0,20,2,16,10,3;0,10,20,9,-15,7];
sigma_d_square=0.001./((3.0*10.^8).^2);
fprintf('                     Table 1                                          ');
fprintf('\n');
fprintf('-----------------------------------------');
fprintf('\n');
fprintf('MSE      M=4      M=5       M=6   ');
fprintf('\n');
fprintf('------------------------------------------');
fprintf('\n');
fprintf('SI ');
 for num_=4:6
    fprintf('%10.5f',SI_method(num_,sensor,location,sigma_d_square));
 end
 fprintf('\n');fprintf('SX ');
  for num_=4:6
    fprintf('%10.5f',SX_method(num_,sensor,location,sigma_d_square));
  end
 
 fprintf('\n');     %set the form