function [ek, ubov, hlsw] = eigensystem(q)

global A_mat

hlsw = sw_hamiltonian(q);
[ubov, ek] = eigenshuffle(A_mat*hlsw);

% para-renormalization of the Bogoliubov vectors

tmp = ubov' * A_mat * ubov;

for k=1:8
  ubov(:, k) = ubov(:, k)/sqrt(abs(tmp(k, k)));
end

ek = 2 * ek;

end
