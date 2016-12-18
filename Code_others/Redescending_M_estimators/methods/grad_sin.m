function grad = grad_sin(H,X,B,c)

      tmp = H.*(X-B);
      
      grad = c* sin( tmp/c);
      
      grad(abs(tmp) > c*pi) = 0;
       

return;