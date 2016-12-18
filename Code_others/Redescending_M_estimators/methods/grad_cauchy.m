function grad = grad_cauchy(H,X,B,c)

      tmp = H.*(X-B);
      tmp_exp = 1 + ( (tmp.^2)/c^2);
      
      grad =  tmp ./ tmp_exp;

return;