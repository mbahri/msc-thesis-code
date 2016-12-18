function M = weight_cauchy(X,B,sigma)

    M = 1 ./ ( 1 + (X-B).^2/sigma^2 );


return;