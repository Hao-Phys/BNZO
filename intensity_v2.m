function [sc_inten, sqw_mat] = intensity_v2(omega, qx, qy, qz)

global S Q_vec 

[q1, q2] = qxytoq12(qx, qy);
q = [q1, q2];
len_omega = length(omega);

chi_mat = zeros(len_omega, 3, 3);

% green function with momentum q and q-Q 
gf_q = greenfunction(omega, q);
gf_qmQ = greenfunction(omega, q-Q_vec);

% green function blocks
gf_q_11 = gf_q(:, 1:4, 1:4);
gf_q_12 = gf_q(:, 1:4, 5:8);
gf_q_21 = gf_q(:, 5:8, 1:4);
gf_q_22 = gf_q(:, 5:8, 5:8);

gf_qmQ_11 = gf_qmQ(:, 1:4, 1:4);
gf_qmQ_12 = gf_qmQ(:, 1:4, 5:8);
gf_qmQ_21 = gf_qmQ(:, 5:8, 1:4);
gf_qmQ_22 = gf_qmQ(:, 5:8, 5:8);

%*** calculation of the in-plane components of chi_{mu_0nu_0}(q,omega)***%

% the relative phase in the in-plane components
phase = sus_phase();

for mu0 = 1:2
  for nu0 = 1:2
      [C_a, C_b, C_c, C_d] = sus_coefficients(mu0, nu0);
      
      for alpha = 1:4
         for beta = 1:4
            tmp = -(S/8.0) .* phase(alpha, beta) ...
               .*(C_a(alpha, beta) .* gf_qmQ_12(:, alpha, beta) ...
                + C_b(alpha, beta) .* gf_qmQ_21(:, alpha, beta) ...
                + C_c(alpha, beta) .* gf_qmQ_22(:, alpha, beta) ...
                + C_d(alpha, beta) .* gf_qmQ_11(:, alpha, beta) );

            chi_mat(:, mu0, nu0) = chi_mat(:, mu0, nu0) + tmp;
          
          end
      end
      
   end
end

%*** calculation of the out-of-plane component chi_{z_0z_0}(q,omega)***%
[C_a, C_b, C_c, C_d] = sus_coefficients(3, 3);

 for alpha = 1:4
    for beta = 1:4
          tmp = -(S/8.0) ...
                *(C_a(alpha, beta) .* gf_q_12(:, alpha, beta) ...
                + C_b(alpha, beta) .* gf_q_21(:, alpha, beta) ...
                + C_c(alpha, beta) .* gf_q_22(:, alpha, beta) ...
                + C_d(alpha, beta) .* gf_q_11(:, alpha, beta) );

          chi_mat(:, 3, 3) = chi_mat(:, 3, 3) + tmp;
     end
 end
 
 % -2 should be included
 sqw_mat = -2*imag(chi_mat);
 
 sc_inten = projector(qx, qy, qz, sqw_mat); 

 
end 

