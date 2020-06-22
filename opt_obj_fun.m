function fun = opt_obj_fun(v)

global unit_idx sub_idx
global S
global J_ex DM_ex h_ext

fun = 0.0;

for bond = 1:14

  u1 = unit_idx(bond,1);
  u2 = unit_idx(bond,2);

  alpha = sub_idx(bond,1);
  beta  = sub_idx(bond,2);

  % pass the variables
  theta1 = v((u1-1)*4+alpha);
  phi1   = v((u1-1)*4+alpha+8);

  theta2 = v((u2-1)*4+beta);
  phi2   = v((u2-1)*4+beta+8);
  
  u1_p = mod(u1, 2) + 1;
  u2_p = mod(u2, 2) + 1;

  theta1_p = v((u1_p-1)*4+alpha);
  phi1_p   = v((u1_p-1)*4+alpha+8);

  theta2_p = v((u2_p-1)*4+beta);
  phi2_p   = v((u2_p-1)*4+beta+8);

  S1 = S * [sin(theta1)*cos(phi1), sin(theta1)*sin(phi1), cos(theta1)];
  S2 = S * [sin(theta2)*cos(phi2); sin(theta2)*sin(phi2); cos(theta2)];

  S1_p = S * [sin(theta1_p)*cos(phi1_p), sin(theta1_p)*sin(phi1_p), cos(theta1_p)];
  S2_p = S * [sin(theta2_p)*cos(phi2_p); sin(theta2_p)*sin(phi2_p); cos(theta2_p)];

  J_bond = reshape(J_ex(bond, :, :), 3, 3);

  fun = fun + S1 * J_bond * S2 ...
      + S1_p * J_bond * S2_p;

  if ((bond >=3) && (bond <=10))
    J_bond = DM_ex;
    fun = fun + S1 * J_bond * S2 ...
        + S1_p * J_bond * S2_p;
  end
  
end

for sub_lat = 1:8
    fun = fun - h_ext * S * cos(v(sub_lat));
end

end

