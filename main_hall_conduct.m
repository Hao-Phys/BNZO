kBT = 0.001:0.01:2.0; 

% epsilon1 = 0.6977;
% epsilon2 = 0.7511;
% epsilon3 = 1.0978;
% epsilon4 = 1.2470;

epsilon1 = 0.6480;
epsilon2 = 0.7773;
epsilon3 = 1.1208;
epsilon4 = 1.2580;

n1 = bose(epsilon1, kBT);
n2 = bose(epsilon2, kBT);
n3 = bose(epsilon3, kBT);
n4 = bose(epsilon4, kBT);

C1 = ctwo(n1) - pi^2/3;
C2 = ctwo(n2) - pi^2/3;
C3 = ctwo(n3) - pi^2/3;
C4 = ctwo(n4) - pi^2/3;

berry1 = -1;
berry2 = 0;
berry3 = 1.1532;
berry4 = -0.15333;

kappa = -kBT .* 2.*pi .* (C1.*berry1+C2.*berry2+C3.*berry3+C4.*berry4);

figure;
plot(kBT,kappa, 'r-');
xlabel('$k_BT$', 'Interpreter', 'latex');
ylabel('$\kappa_{xy}\hbar V/k_B$', 'Interpreter', 'latex');
ylim([0,3])