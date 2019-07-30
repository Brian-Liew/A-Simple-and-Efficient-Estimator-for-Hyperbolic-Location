s1=0.001./((3.0*10.^8).^2);
s2=0.0001./((3.0*10.^8).^2);
s3=0.00001./((3.0*10.^8).^2);
fprintf('                                  Table 4                                          ');
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------');
fprintf('\n');
fprintf('MSE     M=4        M=5       M=6       M=7       M=8       M=9       M=10   ');
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------');
fprintf('\n');
fprintf('B  ');
 for num_=4:10
    fprintf('%10.2f',Taylor_distant_linear(num_,s3));
 end
 fprintf('\n');fprintf('C  ');
  for num_=4:10
    fprintf('%10.2f',proposed_1_distant_linear(num_,s3));
  end
  fprintf('\n');fprintf('E  ');
  for num_=4:10
    fprintf('%10.2f',CRLB_distant_linear(num_,s3));
  end
 fprintf('\n');