function [E] = solve_l1l2(W,lambda)
% SOLVE_L1L2 Vectorized implementation of the prox operator of the l1l2 norm
% 
% Mehdi Bahri
% Imperial College London
% August, 2016

E = W*diag(max(1 - lambda ./ diag(sqrt(W'*W)), 0));
end

% function [E] = solve_l1l2(W,lambda)
% n = size(W,2);
% E = W;
% for i=1:n
%     E(:,i) = solve_l2(W(:,i),lambda);
% end
% end
% 
% function [x] = solve_l2(w,lambda)
% % min lambda |x|_2 + |x-w|_2^2
% 
% % x = max(||w||_2 - lambda, 0) * w / ||w||_2
% nw = norm(w);
% if nw>lambda
%     x = (nw-lambda)*w/nw;
% else
%     x = zeros(length(w),1);
% end
% end