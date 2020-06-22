function [qx, qy] = q12toqxy(q1, q2)
% this function is a transformation of change of basis
% from the basis of S-S lattice to the global reference frame

qx = 0.5 * sqrt(2) * (q1-q2);
qy = 0.5 * sqrt(2) * (q1+q2);
end