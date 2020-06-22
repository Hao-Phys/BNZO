function bond_info()

global theta
global unit_idx sub_idx delta_ij 
global a1 a2 h1 h2 h3 h4
global Q_vec

theta = 0.9842;    % the lattice angle defined in the note
Q_vec = [0.5, 0.5];

a1 = [1; 0];
a2 = [0; 1];
h1 = [0; 0];
h2 = 0.5 * [1; (sin(theta/2)-cos(theta/2))./(cos(theta/2)+sin(theta/2))];
h3 = sin(theta/2)/(cos(theta/2)+sin(theta/2)) * [1; 1];
h4 = 0.5 * [(sin(theta/2)-cos(theta/2))./(cos(theta/2)+sin(theta/2)); 1];

unit_idx = zeros(14, 2);
sub_idx = zeros(14, 2);
delta_ij = zeros(2, 14);

unit_idx(1, 1) = 1;
unit_idx(1, 2) = 1;
sub_idx(1, 1) = 1;
sub_idx(1, 2) = 3;
delta_ij(:, 1) = h3;

unit_idx(2, 1) = 1;
unit_idx(2, 2) = 1;
sub_idx(2, 1) = 2;
sub_idx(2, 2) = 4;
delta_ij(:, 2) = a1-a2-(h2-h4);

unit_idx(3, 1) = 1;
unit_idx(3, 2) = 1;
sub_idx(3, 1) = 1;
sub_idx(3, 2) = 2;
delta_ij(:, 3) = h2;

unit_idx(4, 1) = 1;
unit_idx(4, 2) = 2;
sub_idx(4, 1) = 2;
sub_idx(4, 2) = 3;
delta_ij(:, 4) = h4-a2;

unit_idx(5, 1) = 2;
unit_idx(5, 2) = 2;
sub_idx(5, 1) = 3;
sub_idx(5, 2) = 4;
delta_ij(:, 5) = -h2;

unit_idx(6, 1) = 2;
unit_idx(6, 2) = 1;
sub_idx(6, 1) = 4;
sub_idx(6, 2) = 1;
delta_ij(:, 6) = a2-h4;

unit_idx(7, 1) = 1;
unit_idx(7, 2) = 1;
sub_idx(7, 1) = 4;
sub_idx(7, 2) = 1;
delta_ij(:, 7) = -h4;

unit_idx(8, 1) = 1;
unit_idx(8, 2) = 2;
sub_idx(8, 1) = 1;
sub_idx(8, 2) = 2;
delta_ij(:, 8) = h2-a1;

unit_idx(9, 1) = 2;
unit_idx(9, 2) = 2;
sub_idx(9, 1) = 2;
sub_idx(9, 2) = 3;
delta_ij(:, 9) = h4;

unit_idx(10, 1) = 2;
unit_idx(10, 2) = 1;
sub_idx(10, 1) = 3;
sub_idx(10, 2) = 4;
delta_ij(:, 10) = a1-h2;

unit_idx(11, 1) = 1;
unit_idx(11, 2) = 2;
sub_idx(11, 1) = 1;
sub_idx(11, 2) = 3;
delta_ij(:, 11) = h3-a2;

unit_idx(12, 1) = 2;
unit_idx(12, 2) = 1;
sub_idx(12, 1) = 3;
sub_idx(12, 2) = 1;
delta_ij(:, 12) = a1-h3;

unit_idx(13, 1) = 2;
unit_idx(13, 2) = 1;
sub_idx(13, 1) = 4;
sub_idx(13, 2) = 2;
delta_ij(:, 13) = a2+h2-h4;

unit_idx(14, 1) = 2;
unit_idx(14, 2) = 1;
sub_idx(14, 1) = 2;
sub_idx(14, 2) = 4;
delta_ij(:, 14) = a1+h4-h2;


end



