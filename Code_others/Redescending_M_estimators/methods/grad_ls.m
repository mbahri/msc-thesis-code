function grad = grad_ls(H,X,B,sigma)
 
      
      grad = H.*(X-B);

return;