function [ubov] = eigensystem_bov(q)

global A_mat

[tmp, tmp1] = eigensystem(q);
[tmp, tmp2] = eigensystem(-q);

tmp2 = conj(tmp2);

ubov(1:4, 1:4) = reshape(tmp1(1:4, 1:4), 4, 4);
ubov(5:8, 1:4) = reshape(tmp1(5:8, 1:4), 4, 4);
ubov(1:4, 5:8) = reshape(tmp2(5:8, 1:4), 4, 4);
ubov(5:8, 5:8) = reshape(tmp2(1:4, 1:4), 4, 4);

tmp3 = ubov' * A_mat * ubov;

for x = 1:8
    ubov(:, x) = ubov(:, x)/sqrt(abs(tmp3(x, x)));
end

end


