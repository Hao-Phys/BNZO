
function parameters()

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


J_intra = -0.74513;       %(-1.0,-1.4)
Delta_intra_OP = 0.081608; %(0,0.2)
Delta_intra_IP = 0.16742; %(0,0.2)

J_inter = 0.001374;    %(-0.1,0.1)
Delta_inter = 0.082688; %(0,1)
J_inter_xpy = 0.048482; %0.1; %(-0.4,0.4)
J_inter_xmy = 0.16067; %(-0.4,0.4)             % zero in the experimental order
DM = 0.10662;          %0.2;%-0.04;        %(-0.2,0.2)


J_diag = 0.01669;     %(0.01,0.2)
Delta_diag = 0.3058;  %(0, 1)
J_diag_xpy = 0.045295; %(-0.5, 0.5)   % does not affect the LSW energy at zero field
J_diag_xmy = 0.006486; %(-0.4, 0.4)

g_factor = 2.0;
mu_B = 5.788381E-2;
field = 1.0;  %Telsa
h_ext = g_factor*mu_B*field;  


end
