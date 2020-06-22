function c2 = ctwo(x)
c2 = (1+x) .* (log((1+x)./x)).^2 - (log(x)).^2 - 2 .* polylog(2,-x);
end