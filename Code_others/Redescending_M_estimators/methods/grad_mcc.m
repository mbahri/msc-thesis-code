function grad = grad_mcc(H,X,B,sigma)

      tmp = H.*(X-B);
      tmp_exp = exp(-(tmp.^2)/sigma^2);
      
      grad = tmp_exp.*tmp;

return;