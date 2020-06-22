parameters();
exchange_tensors();
bond_info();
return_value = local_to_global();

if return_value == 0
    return
end

% number of sites
N = 50;
dk = 1.0/N;
k1 = -0.5:dk:0.5-dk;
k2 = -0.5:dk:0.5-dk;

band_energy = zeros(8,N,N);
epsilon = zeros(8,N,N);

% dispersion in the whole BZ

for x1 = 1:N
    for x2 = 1:N        
        q= [k1(x1), k2(x2)];
        [epsilon(:, x1, x2), ubov] = eigensystem(q);
    end
end

[K1, K2] = meshgrid(k1, k2);

band=1; suf1= reshape(real(epsilon(band, : , :)), N, N);
band=2; suf2= reshape(real(epsilon(band, : , :)), N, N);
band=3; suf3= reshape(real(epsilon(band, : , :)), N, N);
band=4; suf4= reshape(real(epsilon(band, : , :)), N, N);
  
figure
surf(K1, K2, suf1)
hold on
surf(K1, K2, suf2)
hold on 
surf(K1, K2, suf3)
hold on
surf(K1, K2, suf4)
hold on
xlabel('q1')
ylabel('q2')
zlabel('meV')
zlim([0.0,2.0])
hold off


% dispersion along q2=0
ecut1 = zeros(8,N);
bogo_vec = zeros(8,8);

for x = 1:N
    q = [k1(x), 0.5];
    [ecut1(:, x), ubov] = eigensystem(q);
     
    criterion = abs(imag(ecut1(:, x)));
    criterion = sum(criterion(:));

     if criterion > 10E-8
      disp('UNSTABLE SYSTEM')
      disp('at the point')
      disp(q)
     end
end

figure 
subplot(2,1,1)
band=1; Cut1= reshape(real(ecut1(band,:)),1,N);
band=2; Cut2= reshape(real(ecut1(band,:)),1,N);
band=3; Cut3= reshape(real(ecut1(band,:)),1,N);
band=4; Cut4= reshape(real(ecut1(band,:)),1,N);
plot(2*pi*k1,Cut1);
hold on
plot(2*pi*k1,Cut2);
hold on
plot(2*pi*k1,Cut3);
hold on
plot(2*pi*k1,Cut4);
hold on
xlim([-pi,pi])
ylim([0.0,2.5])
xlabel('q1')
ylabel('meV')
text(0.0, 2.0, 'q_2=\pi','FontSize', 18)

hold off

for x = 1:N
    q = [0.5, k2(x)];
    [ecut1(:, x), ubov] = eigensystem(q);
     
    criterion = abs(imag(ecut1(:, x)));
    criterion = sum(criterion(:));

     if criterion > 10E-8
      disp('UNSTABLE SYSTEM')
      disp('at the point')
      disp(q)
     end
end

subplot(2,1,2)
band=1; Cut1= reshape(real(ecut1(band,:)),1,N);
band=2; Cut2= reshape(real(ecut1(band,:)),1,N);
band=3; Cut3= reshape(real(ecut1(band,:)),1,N);
band=4; Cut4= reshape(real(ecut1(band,:)),1,N);
plot(2*pi*k2,Cut1);
hold on
plot(2*pi*k2,Cut2);
hold on
plot(2*pi*k2,Cut3);
hold on
plot(2*pi*k2,Cut4);
hold on
xlim([-pi,pi])
ylim([0.0,2.5])
xlabel('q2')
ylabel('meV')
text(0, 2.0, 'q_1=\pi', 'FontSize', 18)
hold off

