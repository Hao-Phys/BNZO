function [q1, q2] = qxytoq12(qx, qy)
% this function is a transformation of change of basis
% from the global frame to the basis of S-S lattice
q1 = 0.5 * sqrt(2) * (qx + qy);
q2 = 0.5 * sqrt(2) * (-qx + qy);
end