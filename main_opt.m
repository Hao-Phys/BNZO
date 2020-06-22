parameters();
exchange_tensors();
bond_info();

x_global_opt = zeros(16, 1);
e_global_opt = 1000.0;

x0 = zeros(16, 1);
lb0 = zeros(16, 1);
ub0 = zeros(16, 1);

lb0(1:16) = -0.001;
ub0(1:8) =  pi + 0.001;
ub0(9:16) = 2*pi + 0.001;

% number of random processes
N_rand = 200;

for flag = 1:N_rand
  rng('shuffle');
  x0(1:8) = pi * rand(8, 1);
  x0(9:16) = 2.0 * pi * rand(8,1);

  options_opt = optimoptions('fmincon','MaxIteration', 5000, 'Display', 'iter');
  problem_opt.objective = @opt_obj_fun;
  problem_opt.x0 = x0;
  problem_opt.lb = lb0;
  problem_opt.ub = ub0;
  problem_opt.solver = 'fmincon';
  problem_opt.options = options_opt;
  [x_opt, e_opt] = fmincon(problem_opt);

  if (e_opt <= e_global_opt)
    x_global_opt = x_opt;
    e_global_opt = e_opt;
    disp('find new local minimum!')
  end

end

disp('write optimal order to the opt_order.txt!')
disp(x_global_opt)
disp('write optimal order to the opt_order.txt!')
disp(x_global_opt)
disp('the optimal ground state energy per site is')
out = [num2str(e_global_opt), ' meV'];
disp(out)

fname = 'opt_order.txt';

fileID = fopen(fname, 'w');
fprintf(fileID, 'theta');
fprintf(fileID, 'phi\n');
for i = 1:8
    fprintf(fileID, '%12.8f', x_global_opt(i));
    fprintf(fileID, '%12.8f\n', x_global_opt(i+8));
end
   


fclose(fileID);
    

