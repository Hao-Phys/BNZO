function [hlsw] = sw_hamiltonian(q)

global unit_idx sub_idx
global delta_ij
global S
global J_ex DM_ex
global R_l2g
global h_ext cant_ang

hlsw = zeros(8,8);
hlsw_11 = zeros(4,4);
hlsw_22 = zeros(4,4);
hlsw_12 = zeros(4,4);
hlsw_21 = zeros(4,4);

for bond = 1:14
    
    u1 = unit_idx(bond,1);
    u2 = unit_idx(bond,2);
    alpha = sub_idx(bond,1);
    beta  = sub_idx(bond,2);
 
    bond_vec = reshape(delta_ij(:,bond),2,1);
    phase = exp(1i*2*pi*q*bond_vec);
    cphase = conj(phase);
    
    J_bond = reshape(J_ex(bond, :, :), 3, 3);
    R1 = reshape(R_l2g((u1-1)*4 + alpha, :, :), 3, 3);
    R2 = reshape(R_l2g((u2-1)*4 + beta, :, :), 3, 3);
    
    J_local = transpose(R1)*J_bond*R2;
    
    J_a = -S * J_local(3,3);
    
    J_b = 0.5 * S * (J_local(1,1) + 1i*J_local(2,1) ...
                   - 1i*J_local(1,2) + J_local(2,2));
    cJ_b = conj(J_b);
    
    J_c = 0.5 * S * (J_local(1,1) - 1i*J_local(2,1) ...
                   - 1i*J_local(1,2) - J_local(2,2));
    cJ_c = conj(J_c);
               
    hlsw_11(alpha, alpha) = hlsw_11(alpha, alpha) + 0.5 * J_a;
    hlsw_11(beta, beta) = hlsw_11(beta, beta) + 0.5 * J_a;
    hlsw_22(alpha, alpha) = hlsw_22(alpha, alpha) + 0.5 * J_a;
    hlsw_22(beta, beta) = hlsw_22(beta, beta) + 0.5 * J_a;

    hlsw_11(alpha, beta) = hlsw_11(alpha, beta) ...
                         + 0.5 * J_b * phase;
    hlsw_11(beta, alpha) = hlsw_11(beta, alpha) ...
                         + 0.5 * cJ_b * cphase; 
    hlsw_22(alpha, beta) = hlsw_22(alpha, beta) ...
                         + 0.5 * cJ_b * phase;
    hlsw_22(beta, alpha) = hlsw_22(beta, alpha) ...
                         + 0.5 * J_b * cphase;

    hlsw_21(alpha, beta) = hlsw_21(alpha, beta) ...
                         + 0.5 * J_c * phase;
    hlsw_21(beta, alpha) = hlsw_21(beta, alpha) ...
                         + 0.5 * J_c * cphase;
    hlsw_12(alpha, beta) = hlsw_12(alpha, beta) ...
                         + 0.5 * cJ_c * phase;
    hlsw_12(beta, alpha) = hlsw_12(beta, alpha) ...
                         + 0.5 * cJ_c * cphase; 
                     
     if ((bond >= 3) && (bond <=10))
          J_local = R1'*DM_ex*R2;
          J_a = -S * J_local(3,3);
          J_b = 0.5 * S * (J_local(1,1) + 1i*J_local(2,1) ...
                    - 1i*J_local(1,2) + J_local(2,2));
          cJ_b = conj(J_b);
          J_c = 0.5 * S * (J_local(1,1) - 1i*J_local(2,1) ...
                    - 1i*J_local(1,2) - J_local(2,2));
          cJ_c = conj(J_c);
                
          hlsw_11(alpha, alpha) = hlsw_11(alpha, alpha) + 0.5 * J_a;
          hlsw_11(beta, beta) = hlsw_11(beta, beta) + 0.5 * J_a;
          hlsw_22(alpha, alpha) = hlsw_22(alpha, alpha) + 0.5 * J_a;
          hlsw_22(beta, beta) = hlsw_22(beta, beta) + 0.5 * J_a;
 
          hlsw_11(alpha, beta) = hlsw_11(alpha, beta) ...
                          + 0.5 * J_b * phase;
          hlsw_11(beta, alpha) = hlsw_11(beta, alpha) ...
                          + 0.5 * cJ_b * cphase; 
          hlsw_22(alpha, beta) = hlsw_22(alpha, beta) ...
                          + 0.5 * cJ_b * phase;
          hlsw_22(beta, alpha) = hlsw_22(beta, alpha) ...
                          + 0.5 * J_b * cphase;
 
          hlsw_21(alpha, beta) = hlsw_21(alpha, beta) ...
                          + 0.5 * J_c * phase;
          hlsw_21(beta, alpha) = hlsw_21(beta, alpha) ...
                          + 0.5 * J_c * cphase;
          hlsw_12(alpha, beta) = hlsw_12(alpha, beta) ...
                          + 0.5 * cJ_c * phase;
          hlsw_12(beta, alpha) = hlsw_12(beta, alpha) ...
                          + 0.5 * cJ_c * cphase; 
     end
        
       
end

hlsw(1:4, 1:4) = hlsw_11;
hlsw(1:4, 5:8) = hlsw_12;
hlsw(5:8, 1:4) = hlsw_21;
hlsw(5:8, 5:8) = hlsw_22;

hlsw = hlsw + 0.5 * h_ext * diag(cos(cant_ang));


end
