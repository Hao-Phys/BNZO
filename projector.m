function inten = projector(qx, qy, qz, sqw_mat)

q = [qx, qy, qz];
norm_sq = qx.^2 + qy.^2 + qz.^2;
inten = 0.0;

if norm_sq == 0
    disp('momentum transfer zero, calculation terminates!');
    return
end

% project out the parallel components

for mu0 = 1:3
    for nu0 = 1:3
     
        if mu0 == nu0
           inten = inten ...
                 + (1-q(mu0)*q(nu0)/norm_sq) .* sqw_mat(:, mu0, nu0);
        else 
            inten = inten ...
                  - (q(mu0)*q(nu0)/norm_sq) .* sqw_mat(:, mu0, nu0);
        end
        
    end
end


end
      