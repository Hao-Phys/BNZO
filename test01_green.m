parameters();
exchange_tensors();
bond_info();
return_value = local_to_global();
global broadening A_mat


if return_value == 0
    return
end

q = [0.45, 0.67];
omega1 = 0.13;
omega2 = 0.33;
omega3 = 0.12;

hsw = sw_hamiltonian(q);
gf1 = -(omega1 + 1i*broadening) * A_mat + 2*hsw;
inv_gf1 = inv(gf1);
gf2 = -(omega2 + 1i*broadening) * A_mat + 2*hsw;
gf3 = -(omega3 + 1i*broadening) * A_mat + 2*hsw;

omega = [omega1, omega2, omega3];
gf4 = greenfunction(omega, q);

re1 = gf1 * reshape(gf4(1, :, :), 8, 8);
re2 = gf2 * reshape(gf4(2, :, :), 8, 8);
re3 = gf3 * reshape(gf4(3, :, :), 8, 8);

  disp(re1)
%  disp(re2)
%  disp(re3)

disp('the hamiltonian is')
disp(2*real(hsw))
disp(2*imag(hsw))

disp('the green function of omega1=0.13 is')
disp('real part')
disp(real(reshape(gf4(1, :, :), 8, 8)))
disp('imag part')
disp(imag(reshape(gf4(1, :, :), 8, 8)))

disp('the result of direct inverse is')
disp('real part')
disp(real(inv_gf1))
disp('imag part')
disp(imag(inv_gf1))


