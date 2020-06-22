function gf = greenfunction(omega, q)

% This function calculates the green function of H-P bosons
% with momentum q and energy range omega
% omega can be either a single energy or a discrete energy range.

global broadening

len_omega = length(omega);

% for each omega, the green function is a 8*8 matrix

gf = zeros(len_omega, 8, 8);

[ek, ubov] = eigensystem(q);
[ek_m, ubov_m] = eigensystem(-q);

% see the note, use the properties of the Bogoliubov 
% coefficients

ubov_m = conj(ubov_m);

minus_mat = zeros(8, 8);
plus_mat = zeros(8, 8);

u11 = reshape(ubov(1:4, 1:4), 4, 4);
u21 = reshape(ubov(5:8, 1:4), 4, 4);
u11_m = reshape(ubov_m(1:4, 1:4), 4, 4);
u21_m = reshape(ubov_m(5:8, 1:4), 4, 4);

for band = 1:4

  minus_mat(1:4, 1:4) = u11(:, band) * (u11(:, band))';
  minus_mat(1:4, 5:8) = u11(:, band) * (u21(:, band))';
  minus_mat(5:8, 1:4) = u21(:, band) * (u11(:, band))';
  minus_mat(5:8, 5:8) = u21(:, band) * (u21(:, band))';

  plus_mat(1:4, 1:4) = u21_m(:, band) * (u21_m(:, band))';
  plus_mat(1:4, 5:8) = u21_m(:, band) * (u11_m(:, band))';
  plus_mat(5:8, 1:4) = u11_m(:, band) * (u21_m(:, band))';
  plus_mat(5:8, 5:8) = u11_m(:, band) * (u11_m(:, band))';
  
  
  for i = 1:8
      for j = 1:8
          
          tmp1 = reshape((omega - ek(band) + 1i * broadening),len_omega,1);
          tmp2 = reshape((omega + ek_m(band) + 1i * broadening),len_omega,1); 
          
          gf(:, i, j) = gf(:,i,j) ...
                       - minus_mat(i, j)./tmp1 + plus_mat(i, j)./tmp2;
         
      end
  end

end
  


end
