lattice_angle = 0.4921;
d = sqrt(2) * (sin(lattice_angle)+cos(lattice_angle));
phi = 3*pi/4 - lattice_angle;

num = sqrt(2)*cos(lattice_angle)/d;
lattice_vec1 = [num, num];

figure
plot(lattice_vec1(1), lattice_vec1(2), 'r*')

