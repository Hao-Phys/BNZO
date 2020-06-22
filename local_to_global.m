function return_val=local_to_global()

global R_l2g h_ext cant_ang
global J_intra Delta_intra_IP 
global J_inter_xmy
global  J_diag_xmy
%J_diag_xpy

cant_ang = zeros(8, 1);

% if external field zero, in-plane order
% if J_inter_xmy = 0, the moments are perpendicular to the dimer
% if J_inter_xmy != 0, the moments gain an uniform canting angle
% phi_cant defined later

% if external field non-zero, use optimization to find the order!

if h_ext == 0
   
   if J_inter_xmy == 0
       
      return_val = 1;
      theta1 = pi/2;
      phi1 = 0;

      theta2 = pi/2;
      phi2 =  3*pi/2;
   
      theta3 = pi/2;
      phi3 = 0;

      theta4 = pi/2;
      phi4 = 3*pi/2;
   
      theta1_p = pi/2;
      phi1_p = pi;

      theta2_p = pi/2;
      phi2_p = pi/2;

      theta3_p = pi/2;
      phi3_p = pi;

      theta4_p = pi/2;
      phi4_p = pi/2;
   
      theta = [theta1, theta2, theta3, theta4, ...
             theta1_p, theta2_p, theta3_p, theta4_p];

      phi = [phi1, phi2, phi3, phi4, ...
             phi1_p, phi2_p, phi3_p, phi4_p];
         
   else 
       
      return_val = 1;
      expr = (8*J_inter_xmy)... 
          /(J_intra*(1-Delta_intra_IP)-4*J_diag_xmy);
          % delete 4*J_diag_xpy Dec 4, 2019 by Hao 
            
      phi_cant = 0.5 * atan(expr);
      
      theta1 = pi/2;
      phi1 = phi_cant;

      theta2 = pi/2;
      phi2 =  3*pi/2 + phi_cant;
   
      theta3 = pi/2;
      phi3 = phi_cant;

      theta4 = pi/2;
      phi4 = 3*pi/2 + phi_cant;
   
      theta1_p = pi/2;
      phi1_p = pi + phi_cant;

      theta2_p = pi/2;
      phi2_p = pi/2 + phi_cant;

      theta3_p = pi/2;
      phi3_p = pi + phi_cant;

      theta4_p = pi/2;
      phi4_p = pi/2 + phi_cant;
   
      theta = [theta1, theta2, theta3, theta4, ...
             theta1_p, theta2_p, theta3_p, theta4_p];

      phi = [phi1, phi2, phi3, phi4, ...
             phi1_p, phi2_p, phi3_p, phi4_p];
         
   end
       
   
else
    reason = isfile('opt_order.txt');
    return_val = reason;
    if reason == 0
        disp('optimal order file does not exist!');
        disp('please run main_opt or main_opt_v2 to generate optimal order');
        return
    end
    
    dat = importdata('opt_order.txt');
    dat = dat.data;
    theta = zeros(8, 1);
    phi = zeros(8, 1);
    
    for flag = 1:8
        theta(flag) = dat(flag, 1);
        phi(flag) = dat(flag, 2);
    end
    
end


% transformation matrix from local to global frame

R_l2g = zeros(8, 3, 3);

for i = 1:8
    R_l2g(i, :, :) = [-sin(phi(i)), -cos(phi(i))*cos(theta(i)), cos(phi(i))*sin(theta(i)); ... 
                      cos(phi(i)), -sin(phi(i))*cos(theta(i)), sin(phi(i))*sin(theta(i)); ...
                      0, sin(theta(i)), cos(theta(i))];
    cant_ang(i) = theta(i);
end


end