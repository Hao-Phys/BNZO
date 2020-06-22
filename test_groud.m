parameters();
exchange_tensors();
bond_info();

data1 = importdata('opt_order.txt');
data1 = data1.data;
data2 = importdata('../src_field/optimal_mean_field.dat');


v1(1:8) = data1(1:8, 2);
v1(9) = data1(1, 1);

v2(1:8) = data2(1:8, 2);
v2(9) = data2(1, 1);

ecl1 = classical_energy(v1);
ecl2 = classical_energy(v2);

disp(ecl1)
disp(ecl2)