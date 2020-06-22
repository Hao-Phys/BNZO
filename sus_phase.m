function phase = sus_phase()
global h1 h2 h3 h4
global Q_vec

phase = zeros(4, 4);
h_vec = zeros(2, 4);
h_vec(:, 1) = h1;
h_vec(:, 2) = h2;
h_vec(:, 3) = h3;
h_vec(:, 4) = h4;

for alpha = 1:4
  for beta = 1:4
    phase(alpha, beta) = exp(-1i*2*pi*Q_vec*(h_vec(:, alpha)-h_vec(:, beta)));
  end
end 

end

