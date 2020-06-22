J_intra = -1.2;
Delta_intra_OP = 0.2;
Delta_intra_IP = 0.1;

J_inter = 0.03;
Delta_inter = 0.5;
J_inter_xpy = 0.1;
J_inter_xmy = 0.2; % zero in the experimental order
DM = -0.04;


J_diag = 0.06;
Delta_diag = 0.5;
J_diag_xpy = -0.4;    % does not affect the LSW energy at zero field
J_diag_xmy = 0.02;

h_ext = 0.0;

phi = 0.9592816472;
t2phi = tan(2*phi);

expr = (8*J_inter_xmy-4*J_diag_xpy)/(J_intra*(1-Delta_intra_IP)-4*J_diag_xmy);

disp(expr)
disp(t2phi)