parameters();
exchange_tensors();
bond_info();
return_value = local_to_global();

if return_value == 0
    return
end

Nw = 100;
w_min = 0.0;
w_max = 1.0;
dw = (w_max-w_min)/Nw;
omega = w_min+dw:dw:w_max;

BZpath = importdata('BZpath.dat');
kk = BZpath(:, 1);
Nk = length(kk);
%Nk = 200;
%dk = (1.8-0.2)/Nk;
%kk = 0.2:dk:1.8;     % module of input k
%kk = kk*6.77;        % units: convert in (1/a) for calculation a=6.77A is
                      % the lattice constant

intensity_sc = zeros(Nw, Nk);
sqw_mat = zeros(Nw, Nk, 3, 3);
kz = 0.0;


for flag = 1:Nk
    q1 = BZpath(flag, 2);
    q2 = BZpath(flag, 3);
    [kx, ky] = q12toqxy(q1, q2);
    [intensity_sc(:, flag), sqw_mat(:, flag, :, :)] ... 
    = intensity_v2(omega, kx, ky, kz);
end

Sxx = reshape(sqw_mat(:, :, 1, 1), Nw, Nk);
Syy = reshape(sqw_mat(:, :, 2, 2), Nw, Nk);
Szz = reshape(sqw_mat(:, :, 3, 3), Nw, Nk);
Sxy = reshape(sqw_mat(:, :, 1, 2), Nw, Nk);
Syx = reshape(sqw_mat(:, :, 2, 1), Nw, Nk);

figure


[Q, O] = meshgrid(kk, omega);
subplot(2, 3, 6);
pcolor(Q, O, intensity_sc)
ylim([0,0.9])
caxis([0, 20])
shading interp
title('total intensity')
colorbar

subplot(2, 3, 1)
pcolor(Q, O, Sxx)
ylim([0,0.9])
caxis([0, 20])
shading interp
title('Sxx')
colorbar

subplot(2, 3, 2)
pcolor(Q, O, Syy)
ylim([0,0.9])
caxis([0, 20])
shading interp
title('Syy')
colorbar

subplot(2, 3, 3)
pcolor(Q, O, Szz)
ylim([0,0.9])
caxis([0, 20])
shading interp
title('Szz')
colorbar

subplot(2, 3, 4)
pcolor(Q, O, Sxy)
ylim([0,0.9])
caxis([-20, 20])
shading interp
title('Sxy')
colorbar

subplot(2, 3, 5)
pcolor(Q, O, Syx)
ylim([0,0.9])
caxis([-20, 20])
shading interp
title('Syx')
colorbar

% figure 
% band=1; Cut1= reshape(real(ek_q(band,:)),1,Nk1);
% band=2; Cut2= reshape(real(ek_q(band,:)),1,Nk1);
% band=3; Cut3= reshape(real(ek_q(band,:)),1,Nk1);
% band=4; Cut4= reshape(real(ek_q(band,:)),1,Nk1);
% plot(k1,Cut1,'r-');
% hold on
% plot(k1,Cut2,'r-');
% hold on
% plot(k1,Cut3,'r-');
% hold on
% plot(k1,Cut4,'r-');
% hold on
% band=1; Cut1= reshape(real(ek_qmQ(band,:)),1,Nk1);
% band=2; Cut2= reshape(real(ek_qmQ(band,:)),1,Nk1);
% band=3; Cut3= reshape(real(ek_qmQ(band,:)),1,Nk1);
% band=4; Cut4= reshape(real(ek_qmQ(band,:)),1,Nk1);
% plot(k1,Cut1,'b-');
% hold on
% plot(k1,Cut2,'b-');
% hold on
% plot(k1,Cut3,'b-');
% hold on
% plot(k1,Cut4,'b-');
% hold on
% ylim([0,1])
% xlabel('k1')
% ylabel('meV')
% hold off


