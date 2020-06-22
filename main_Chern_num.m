parameters();
exchange_tensors();
bond_info();
return_value = local_to_global();

if return_value == 0
    return
end

global dk

N = 100;
dk = 1.0/N;
k1 = -0.5:dk:0.5-dk;
k2 = -0.5:dk:0.5-dk;

%******* calculation of Berry curvature and Chern number *****%


BerryCurvature = zeros(8, N, N);

for i=1:(N)
    for j=1:(N)
      q = [k1(i), k2(j)];
      BerryCurvature(:, i, j) = berrycurvature(q);
    end  
end


berryC1 = reshape(BerryCurvature(1,:,:), N, N);
berryC2 = reshape(BerryCurvature(2,:,:), N, N);
berryC3 = reshape(BerryCurvature(3,:,:), N, N);
berryC4 = reshape(BerryCurvature(4,:,:), N, N);
% berryC5 = reshape(BerryCurvature(5,:,:),N,N);
% berryC6 = reshape(BerryCurvature(6,:,:),N,N);
% berryC7 = reshape(BerryCurvature(7,:,:),N,N);
% berryC8 = reshape(BerryCurvature(8,:,:),N,N);

pos1=find(abs(berryC1(:))>10^4);
berryC1(pos1)=0;
Chern_num1 = (dk)^2*sum(berryC1(:))*(2*pi);

pos2=find(abs(berryC2(:))>10^4);
berryC2(pos2)=0;
Chern_num2 = (dk)^2*sum(berryC2(:))*(2*pi);

pos3=find(abs(berryC3(:))>10^4);
berryC3(pos3)=0;
Chern_num3 = (dk)^2*sum(berryC3(:))*(2*pi);

pos4=find(abs(berryC4(:))>10^4);
berryC4(pos4)=0;
Chern_num4 = (dk)^2*sum(berryC4(:))*(2*pi);

% pos5=find(abs(berryC5(:))>10^4);
% berryC5(pos3)=0;
% Chern_num5 = (dk)^2*sum(berryC5(:))*(2*pi);
% 
% pos6=find(abs(berryC6(:))>10^4);
% berryC6(pos6)=0;
% Chern_num6 = (dk)^2*sum(berryC6(:))*(2*pi);
% 
% pos7=find(abs(berryC7(:))>10^4);
% berryC7(pos7)=0;
% Chern_num7 = (dk)^2*sum(berryC7(:))*(2*pi);
% 
% pos8=find(abs(berryC8(:))>10^4);
% berryC8(pos8)=0;
% Chern_num8 = (dk)^2*sum(berryC8(:))*(2*pi);

% **********plots***********

[K1, K2]=meshgrid(k1, k2);

figure 
subplot(2,2,1)
pcolor(K1,K2,berryC1)
colormap(redblue)
shading interp
colorbar
caxis([-10,10])
title(['band4, $\frac{1}{2\pi}\int d^2k \Omega_{n,k}=$', num2str(Chern_num1)], ...
    'Interpreter', 'latex','FontSize', 16, 'Color', 'k')

subplot(2,2,2)
pcolor(K1,K2,berryC2)
colormap(redblue)
shading interp
colorbar
caxis([-10,10])
title(['band3, $\frac{1}{2\pi}\int d^2k \Omega_{n,k}=$', num2str(Chern_num2)], ...
    'Interpreter', 'latex','FontSize', 16, 'Color', 'k')

subplot(2,2,3)
pcolor(K1,K2,berryC3)
colormap(redblue)
shading interp
colorbar
caxis([-10,10])
title(['band2, $\frac{1}{2\pi}\int d^2k \Omega_{n,k}=$', num2str(Chern_num3)], ...
    'Interpreter', 'latex','FontSize', 16, 'Color', 'k')

subplot(2,2,4)
pcolor(K1,K2,berryC4)
colormap(redblue)
shading interp
colorbar
caxis([-10,10])
title(['band1, $\frac{1}{2\pi}\int d^2k \Omega_{n,k}=$', num2str(Chern_num4)], ...
    'Interpreter', 'latex','FontSize', 16, 'Color', 'k')
