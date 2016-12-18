function M = weight_welsch(X,B,sigma)

    M = exp( - (X-B).^2/sigma^2 );


return;