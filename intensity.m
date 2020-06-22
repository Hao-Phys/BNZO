function sc_inten = intensity(omega, qx, qy, qz)

% This function calculates the single-crystal intensity
% (qx, qy, qz) is the momentum transfer of the neutron 
% in the global reference frame
% omega is the energy range of the neutron

global S Q_vec broadening

[q1, q2] = qxytoq12(qx, qy);
q = [q1, q2];
len_omega = length(omega);

chi_mat = zeros(len_omega, 3, 3);

phase = sus_phase();

[eq, ubov] = eigensystem(q);
[emq, ubov_m] = eigensystem(-q);
[eqmQ, ubov_3] = eigensystem(q-Q_vec);
[eQmq, ubov_4] = eigensystem(Q_vec-q);

ubov_m = conj(ubov_m);
ubov_4 = conj(ubov_4);

u11 = reshape(ubov(1:4, 1:4), 4, 4);
u21 = reshape(ubov(5:8, 1:4), 4, 4);
u11_m = reshape(ubov_m(1:4, 1:4), 4, 4);
u21_m = reshape(ubov_m(5:8, 1:4), 4, 4);
u11_3 = reshape(ubov_3(1:4, 1:4), 4, 4);
u21_3 = reshape(ubov_3(5:8, 1:4), 4, 4);
u11_4 = reshape(ubov_4(1:4, 1:4), 4, 4);
u21_4 = reshape(ubov_4(5:8, 1:4), 4, 4);


for band = 1:4

  minus_mat_11 = u11(:, band) * (u11(:, band))';
  minus_mat_12 = u11(:, band) * (u21(:, band))';
  minus_mat_21 = u21(:, band) * (u11(:, band))';
  minus_mat_22 = u21(:, band) * (u21(:, band))';

  plus_mat_11 = u21_m(:, band) * (u21_m(:, band))';
  plus_mat_12 = u21_m(:, band) * (u11_m(:, band))';
  plus_mat_21 = u11_m(:, band) * (u21_m(:, band))';
  plus_mat_22 = u11_m(:, band) * (u11_m(:, band))';

  minus_mat_Q_11 = u11_3(:, band) * (u11_3(:, band))';
  minus_mat_Q_12 = u11_3(:, band) * (u21_3(:, band))';
  minus_mat_Q_21 = u21_3(:, band) * (u11_3(:, band))';
  minus_mat_Q_22 = u21_3(:, band) * (u21_3(:, band))';

  plus_mat_Q_11 = u21_4(:, band) * (u21_4(:, band))';
  plus_mat_Q_12 = u21_4(:, band) * (u11_4(:, band))';
  plus_mat_Q_21 = u11_4(:, band) * (u21_4(:, band))';
  plus_mat_Q_22 = u11_4(:, band) * (u11_4(:, band))';

  for mu0 = 1:2
    for nu0 = 1:2

      [C_a, C_b, C_c, C_d] = sus_coefficients(mu0, nu0);

      tmp1 = (-2.0*S) .* phase .* (C_a .* minus_mat_Q_12 + C_b .* minus_mat_Q_21 ...
           + C_c .* minus_mat_Q_22 + C_d .* minus_mat_Q_11);
      tmp2 = (-2.0*S) .* phase .* (C_a .* plus_mat_Q_12 + C_b .* plus_mat_Q_21 ...
           + C_c .* plus_mat_Q_22 + C_d .* plus_mat_Q_11);
      
      temp = reshape(- sum(tmp1(:)) ./ (omega - eqmQ(band) + 1i*broadening) ...
               + sum(tmp2(:)) ./ (omega - eQmq(band) + 1i*broadening), ...
               len_omega, 1);
           
      chi_mat(:, mu0, nu0) = chi_mat(:, mu0, nu0) + temp;
                       
    end
  end 

  [C_a, C_b, C_c, C_d] = sus_coefficients(3, 3);

  tmp1 = (-2.0*S) .* (C_a .* minus_mat_12 + C_b .* minus_mat_21 ...
      + C_c .* minus_mat_22 + C_d .* minus_mat_11);
  tmp2 = (-2.0*S) .* (C_a .* plus_mat_12 + C_b .* plus_mat_21 ...
      + C_c .* plus_mat_22 + C_d .* plus_mat_11);
  
  temp = reshape(- sum(tmp1(:)) ./ (omega - eq(band) + 1i*broadening) ...
               + sum(tmp2(:)) ./ (omega - emq(band) + 1i*broadening), ...
               len_omega, 1);
           
  chi_mat(:, 3, 3) = chi_mat(:, 3, 3) + temp;

end
  
sqw_mat = -2 * imag(chi_mat);
% pro_mat = projector(qx, qy, qz);
% 
% tmp = reshape(pro_mat .* sqw_mat;
sc_inten = sqw_mat;

end

