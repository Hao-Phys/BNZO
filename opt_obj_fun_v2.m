function fun = opt_obj_fun_v2(v)

global unit_idx sub_idx
global S
global J_ex DM_ex h_ext

fun = 0.0;

theta = v(9);

for bond = 1:14

  u1 = unit_idx(bond,1);
  u2 = unit_idx(bond,2);

  alpha = sub_idx(bond,1);
  beta  = sub_idx(bond,2);
  
  J_bond = reshape(J_ex(bond, :, :), 3, 3);

  % pass the variables
  phi1   = v((u1-1)*4+alpha);
  phi2   = v((u2-1)*4+beta);

  S1 = S * [sin(theta)*cos(phi1), sin(theta)*sin(phi1), cos(theta)];
  S2 = S * [sin(theta)*cos(phi2); sin(theta)*sin(phi2); cos(theta)];

  u1_p = mod(u1, 2) + 1;
  u2_p = mod(u2, 2) + 1;

  phi1_p   = v((u1_p-1)*4+alpha);
  phi2_p   = v((u2_p-1)*4+beta);

  S1_p = S * [sin(theta)*cos(phi1_p), sin(theta)*sin(phi1_p), cos(theta)];
  S2_p = S * [sin(theta)*cos(phi2_p); sin(theta)*sin(phi2_p); cos(theta)];

  fun = fun + S1 * J_bond * S2 ...
      + S1_p * J_bond * S2_p;

  if ((bond >=3) && (bond <=10))
     J_bond = DM_ex;
     fun = fun + S1 * J_bond * S2 ...
         + S1_p * J_bond * S2_p;
     
  end

end

fun = fun - 8*h_ext*S*cos(theta); 

end
