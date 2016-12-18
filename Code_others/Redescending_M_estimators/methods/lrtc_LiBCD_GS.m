
function [tmp_X, iter, ranks] = lrtc_LiBCD_GS(B, H, X_cell, lambda,  sigma, alpha, grad_f, thresholding_operator,eps,max_iter)
% This function solves the tensor completion problem
% min_{X_{(1)},\ldots,X_{(N)}} J( H.* ( \sum X_{(i)}  ) -B ) + \lambda sum
% r(X_{(i)}
% via the block coordinate descent type methods.
% The implementation of tensor foldings and unfoldings relies on the matlab
% package Tensorlab, which can be downloaded from
% http://www.tensorlab.net/.
% low rank tensor completion
% B: the observed date
% H: the projection operator, represented as a Omega \circ B, where
% \circ is the Hardamard operator
% X_cell: cell array contains the initial guess of X
% lambda: regularization or constrained parameters
% sigma: parameter of the loss functions
% grad_f: functional handle of the gradient of a specified loss
% thresholding_operator: functional handle of a specified thresholding
% operator
% eps: stopping criterion
% max_iter: maximum iterations

% Yang Y., Feng Y., Suykens J.A.K., ``Robust Low Rank Tensor Recovery with Regularized Redescending M-Estimator'', Internal Report 14-97, ESAT-SISTA, KU Leuven (Leuven, Belgium), 2014.

size_tens = size(B);

len_lambda = length(lambda);
len_X = length(X_cell);

if len_lambda ~= len_X
    error('number of parameters not equal!');
end


idx = 1: ndims(B);



X0 = cell(size(X_cell));


err  = 1;
iter = 0;

% auxiliary variable used to compute the gradient
% it also acts as the output tensor X
tmp_X = X_cell{1};
for i = 2: len_X
    tmp_X = tmp_X + X_cell{i};
end


while (err   > eps ) && (iter <=max_iter )
    iter = iter+1;
    X0 = X_cell;
    % Added by Mehdi Bahri
    fprintf('[%s][%s] - It %d - err = %f\n', func2str(grad_f), func2str(thresholding_operator), iter, err);
    
    % update coordinatively along each mode of the tensor
    for i = 1: len_X
         
        Yi =  X0{i}   - (1/alpha) *   grad_f(H,tmp_X,B,sigma);
        
        
        % unfolding the tnesor Yi along its i-th mode
        mB = tens2mat(Yi, i, idx(idx~=i));
        % perform the thresholding operator to the unfolding matrix
        mD = thresholding_operator(mB,lambda(i)/alpha);
        
        % fold the matrix to tensor
        X_cell{i} = mat2tens(mD,size_tens, i,idx(idx~=i));
        
        
        % update the auxiliary variable from X^{k+1}_1 + ... + X^{k+1}_{i-1} +
        % X^k_i + ... + X^k_N to
        % X^{k+1}_1 + ... + X^{k+1}_{i-1} +
        % X^{k+1}_i + ... + X^k_N
        tmp_X = tmp_X - X0{i} + X_cell{i};
        
        
        
    end
    

    err = frob(X_cell{1}-X0{1});
       
end

ranks = zeros(1,len_X);
 for i = 1: len_X
     ranks(i) = rank( tens2mat(X_cell{i},i,idx(idx~=i) ));
 end
% X = X_cell{1};
% for i = 2: len_X
%     X = X + X_cell{i};
% end


return;


