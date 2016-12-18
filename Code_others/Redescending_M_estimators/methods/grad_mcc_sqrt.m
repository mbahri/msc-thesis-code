function grad = grad_mcc_sqrt(H,X,B,sigma)

      tmp = H.*(X-B);
      tmp_exp = exp(-(tmp.^2)/sigma^2);
      
      grad = 1/(2*sigma)* tmp_exp ./ sqrt(1- tmp_exp) .* tmp;

      grad(isnan(grad)==1) = 0;
      
      grad
return;