function X = MeaNoise(Noise)
%  ���������ɲ�����������Ӹ�˹�ֲ�
%  MeaNoise
%    ����˵����
%        Noise:   ��˹�ֲ�������Թ���ƽ���Ľ��
%  Also see: MeaNoise.


%  �������:
if nargout>1,
    error('Too many output arguments!');
end 
if nargin ~= 1
    error('input arguments error!');
end

% ������
Dev = Noise;

% �����
X = sqrt(Dev)*randn(1);

% ��������
if nargout == 1,
    X;
else
    disp(X);
end
