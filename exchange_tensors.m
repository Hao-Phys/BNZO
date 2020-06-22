function exchange_tensors()

global J_intra Delta_intra_OP Delta_intra_IP
global J_inter Delta_inter J_inter_xpy J_inter_xmy DM
global J_diag Delta_diag J_diag_xpy J_diag_xmy

global J_ex DM_ex

% the 10 exchange tensors

J_intra_13 = J_intra * [1, 0, 0; ...
                        0, Delta_intra_IP, 0; ...
                        0, 0, Delta_intra_OP];
                       
J_intra_24 = J_intra * [Delta_intra_IP, 0, 0; ...
                        0, 1, 0; ...
                        0, 0, Delta_intra_OP];
                       
J_inter_a = [J_inter+J_inter_xmy, J_inter_xpy, 0; ...
             J_inter_xpy, J_inter-J_inter_xmy, 0; ...
             0, 0, J_inter*Delta_inter];
             
J_inter_b = [J_inter+J_inter_xmy, -J_inter_xpy, 0; ...
             -J_inter_xpy, J_inter-J_inter_xmy, 0; ...
             0, 0, J_inter*Delta_inter];
   
J_inter_c = [J_inter-J_inter_xmy, -J_inter_xpy, 0; ...
             -J_inter_xpy, J_inter+J_inter_xmy, 0; ...
             0, 0, J_inter*Delta_inter];
         
J_inter_d = [J_inter-J_inter_xmy, J_inter_xpy, 0; ...
             J_inter_xpy, J_inter+J_inter_xmy, 0; ...
             0, 0, J_inter*Delta_inter];
         
J_diag_a = [J_diag+J_diag_xmy, J_diag_xpy, 0; ...
             J_diag_xpy, J_diag-J_diag_xmy, 0; ...
             0, 0, J_diag*Delta_diag];
             
J_diag_b = [J_diag+J_diag_xmy, -J_diag_xpy, 0; ...
             -J_diag_xpy, J_diag-J_diag_xmy, 0; ...
             0, 0, J_diag*Delta_diag];
   
J_diag_c = [J_diag-J_diag_xmy, -J_diag_xpy, 0; ...
             -J_diag_xpy, J_diag+J_diag_xmy, 0; ...
             0, 0, J_diag*Delta_diag];
         
J_diag_d = [J_diag-J_diag_xmy, J_diag_xpy, 0; ...
             J_diag_xpy, J_diag+J_diag_xmy, 0; ...
             0, 0, J_diag*Delta_diag];
         
         
                
J_ex = zeros(14, 3, 3);

J_ex(1, :, :) = J_intra_13;                       
J_ex(2, :, :) = J_intra_24;
J_ex(3, :, :) = J_inter_a;            
J_ex(4, :, :) = J_inter_c;
J_ex(5, :, :) = J_inter_a;
J_ex(6, :, :) = J_inter_c;
J_ex(7, :, :) = J_inter_b;
J_ex(8, :, :) = J_inter_d;
J_ex(9, :, :) = J_inter_b;
J_ex(10, :, :) = J_inter_d;
J_ex(11, :, :) = J_diag_a;
J_ex(12, :, :) = J_diag_b;
J_ex(13, :, :) = J_diag_c;
J_ex(14, :, :) = J_diag_d;

DM_ex = DM * [0, 1, 0; -1, 0, 0; 0, 0, 0];
             