

% This file tests the role of DM and the anisotropic exchange interaction

global S
global J_intra Delta_intra_OP Delta_intra_IP
global J_inter Delta_inter J_inter_xpy J_inter_xmy DM
global J_diag Delta_diag J_diag_xpy J_diag_xmy
global h_ext
global A_mat broadening

S = 0.5;
broadening = 0.005;
A_mat = zeros(8,8);

for i=1:4
    A_mat(i,i) = 1;
end
for i=5:8
    A_mat(i,i) = -1;
end

J_intra = -1.2;
Delta_intra_OP = 0.5;
Delta_intra_IP = 0.5;

J_inter = 0.3;
Delta_inter = 0.5;
J_inter_xpy = 0.0;
J_inter_xmy = 0.0; % zero in the experimental order
DM = 0.0;


J_diag = 0.3;
Delta_diag = 0.5;
J_diag_xpy = 0.0;    % does not affect the LSW energy at zero field
J_diag_xmy = 0.0;

h_ext = 0.2;

exchange_tensors();
bond_info();
return_value = local_to_global();

q = [0, 0.5];
hsw1 = sw_hamiltonian(q);
[ubov1] = eigensystem_bov(q);

%disp(ubov1'*hsw1*ubov1);

J_inter_xpy = 0.2;
DM = 0.3;

exchange_tensors();
bond_info();
return_value1 = local_to_global();

q = [0, 0.5];
[hsw2] = sw_hamiltonian(q);
hdm = hsw2-hsw1;
hdm_bogo=ubov1'*hdm*ubov1;
disp(hdm_bogo);

h_proj(1:2, 1:2) = hdm_bogo(3:4, 3:4);
h_proj(1:2, 3:4) = hdm_bogo(3:4, 7:8);
h_proj(3:4, 1:2) = hdm_bogo(7:8, 3:4);
h_proj(3:4, 3:4) = hdm_bogo(7:8, 7:8);

Asp = diag([1, 1, -1, -1]);
ek = eig(Asp*h_proj);
disp(ek);


